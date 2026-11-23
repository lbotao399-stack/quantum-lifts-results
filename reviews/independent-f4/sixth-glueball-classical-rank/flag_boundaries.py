"""Classical commutator boundaries in Taylor type (4,4,3,3), evaluated on scalar flags."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from functools import cache
import gzip
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
HEIGHTS = (4, 4, 3, 3)


@cache
def trace_words():
    signed, representatives = {}, {}
    first = 1
    for n in range(2, 9):
        mask = 4 ** n - 1
        codes = np.arange(mask + 1, dtype=np.int32)
        reverse = np.zeros_like(codes)
        for j in range(n):
            reverse = 4 * reverse + ((codes >> (2*j)) & 3)

        def least_rotation(a):
            low = a.copy()
            for _ in range(n-1):
                a = ((a << 2) & mask) + (a >> (2*(n-1)))
                low = np.minimum(low, a)
            return low

        forward = least_rotation(codes)
        backward = least_rotation(reverse)
        canonical = np.minimum(forward, backward)
        sign = np.ones_like(codes) if n % 2 == 0 else np.sign(backward-forward)
        representatives[n] = np.flatnonzero((codes == canonical) & (sign != 0))
        order = np.zeros(len(codes), dtype=np.int32)
        order[representatives[n]] = np.arange(first, first+len(representatives[n]))
        signed[n] = sign * order[canonical]
        first += len(representatives[n])
    assert first < 1 << 14
    return signed, representatives, first


@cache
def adjoint():
    data = json.loads(gzip.decompress((HERE/"inputs.json.gz").read_bytes()))
    matrices = np.zeros((52,52,52), dtype=np.int64)
    for i,j,terms in data["brackets"]:
        for k,v in terms:
            matrices[i,k,j] = v
            matrices[j,k,i] = -v
    return matrices


def sample(seed, prime, integral=False):
    rng = np.random.Generator(np.random.PCG64(seed))
    coefficients = rng.integers(-3, 4, size=(4, 52))
    x = (coefficients @ adjoint().reshape(52,-1)).reshape(4,52,52) % prime
    # Products and flattened trace pairings have integer magnitude at most
    # 52^2 (prime-1)^2. All their binary64 sums are exact under this bound.
    assert 52**2 * (prime-1)**2 < 2**53
    if not integral:
        x = x.astype(np.float64)
    powers = {1:x}
    for d in range(2,5):
        powers[d] = ((powers[d-1][:,None] @ x[None]) % prime).reshape(-1,52,52)
    _, reps, size = trace_words()
    result = np.zeros(size, dtype=np.uint16)
    start = 1
    for n, codes in reps.items():
        left, right = powers[n//2], powers[n-n//2]
        if integral:
            # Use separate indexed multiplications for the independent check.
            values = [np.sum(left[k//len(right)] * right[k%len(right)].T) for k in codes]
        else:
            full = left.reshape(len(left),-1) @ right.swapaxes(1,2).reshape(len(right),-1).T
            values = full.ravel()[codes]
        result[start:start+len(codes)] = np.asarray(values, dtype=np.int64) % prime
        start += len(codes)
    return result


@cache
def flags():
    options = []
    for h in HEIGHTS:
        group = []
        for perm in permutations(range(h)):
            sg = (-1)**sum(x > y for i,x in enumerate(perm) for y in perm[i+1:])
            group.append((perm, sg))
        options.append(group)
    rows, signs = [], []
    for a,b,c,d in product(*options):
        rows.append(a[0]+b[0]+c[0]+d[0])
        signs.append(a[1]*b[1]*c[1]*d[1])
    assert len(rows) == 20736
    return np.asarray(rows,dtype=np.int32).T.copy(), np.asarray(signs,dtype=np.int32)


@cache
def patterns():
    def split(total, ceiling, remaining):
        if not total:
            yield ()
        elif remaining:
            for k in range(min(ceiling, total),1,-1):
                for tail in split(total-k,k,remaining-1):
                    yield (k,)+tail
    return tuple(shape for shape in split(14,8,4) if max(shape) >= 3)


def recipe(seed):
    sizes = patterns()[seed % len(patterns())]
    rng = np.random.Generator(np.random.PCG64(6364433+seed))
    slots = rng.permutation(14)
    groups = tuple(np.split(slots,np.cumsum(sizes)[:-1]))
    eligible = [i for i, a in enumerate(groups) if len(a) >= 3]
    mark = eligible[rng.integers(len(eligible))]
    return groups, mark


def polynomial(seed):
    digits, signs = flags()
    signed, _, _ = trace_words()
    factors, marked = recipe(seed)

    def trace(slots):
        key = digits[slots[0]]
        for slot in slots[1:]:
            key = 4*key + digits[slot]
        return signed[len(slots)][key]

    def monomials(factor_ids):
        sg = signs.copy()
        for t in factor_ids:
            sg *= np.sign(t)
        keep = sg != 0
        ids = np.sort(np.abs(np.stack(factor_ids)[:,keep]),axis=0)
        keys = np.zeros(ids.shape[1],dtype=np.int64)
        for t in ids:
            keys = (keys << 14) + t
        return keys, sg[keep]

    ids = [trace(a) for a in factors]
    left, a = monomials(ids)
    swapped = factors[marked].copy()
    swapped[0], swapped[1] = swapped[1], swapped[0]
    ids[marked] = trace(swapped)
    right, b = monomials(ids)
    keys, inverse = np.unique(np.concatenate((left,right)),return_inverse=True)
    coefficients = np.bincount(inverse,weights=np.concatenate((a,-b))).astype(np.int64)
    keep = coefficients != 0
    keys = keys[keep]
    letters = np.stack([(keys >> (14*i)) & ((1<<14)-1) for i in range(len(factors))])
    return coefficients[keep], letters


def evaluate(seed, frame, prime):
    coeff, words = polynomial(seed)
    row = np.zeros(frame.shape[1],dtype=np.int64)
    for j in range(0,len(coeff),64):
        part = words[:,j:j+64]
        values = frame[part[0]].astype(np.int64)
        for t in part[1:]:
            values *= frame[t]
        row += coeff[j:j+64] @ (values % prime)
    return row % prime


def direct_control(seed, sample_seed, prime):
    # Re-evaluate a young symmetrization without cyclic or reversal reduction.
    digits, signs = flags()
    coeff = np.random.Generator(np.random.PCG64(sample_seed)).integers(-3,4,size=(4,52))
    x = (coeff @ adjoint().reshape(52,-1)).reshape(4,52,52) % prime
    groups, mark = recipe(seed)
    @cache
    def raw_trace(word):
        if len(word) == 1:
            return x[word[0]]
        return (raw_trace(word[:-1]) @ x[word[-1]]) % prime
    @cache
    def tr(word):
        return int(raw_trace(word).trace()) % prime
    total = 0
    for vals, sign in zip(digits.T, signs):
        score = int(sign)
        for k,slots in enumerate(groups):
            word = tuple(vals[slots])
            value = tr(word)
            if k == mark:
                value -= tr((word[1],word[0],*word[2:]))
            score = score * value % prime
        total += score
    return total % prime


def calculate(count, probes, prime, workers, replay=None, check=False, checkpoint=None):
    assert prime in (1009,1013) and count > 0 and 0 < probes <= 704
    start = 144336
    trace_words()
    adjoint()
    flags()
    with ThreadPoolExecutor(workers) as pool:
        frame = np.stack(list(pool.map(lambda i: sample(start+i,prime),range(probes))),axis=1)
        basis, selected, pivots = {}, [], []
        determinant = 1
        seeds = list(range(count)) if replay is None else replay
        for lo in range(0,len(seeds),32):
            batch = seeds[lo:lo+32]
            for seed, row in zip(batch,pool.map(lambda i:evaluate(i,frame,prime),batch)):
                for k in sorted(basis):
                    if row[k]:
                        row = (row-row[k]*basis[k]) % prime
                nonzero = np.flatnonzero(row)
                if len(nonzero):
                    k = int(nonzero[0])
                    y = int(row[k])
                    determinant = determinant*y % prime
                    selected.append(seed)
                    pivots.append(k)
                    basis[k] = row*pow(y,-1,prime) % prime
            if checkpoint is not None and (lo+len(batch)) % 128 == 0:
                partial = dict(prime=prime, points=probes, first_point_seed=start, rank=len(selected),
                               recipes=selected, pivot_points=pivots, determinant=determinant,
                               recipes_processed=lo+len(batch))
                temporary = checkpoint.with_suffix(".part")
                temporary.write_text(json.dumps(partial)+"\n")
                temporary.replace(checkpoint)
                print(lo+len(batch), len(selected), flush=True)
            if len(basis) == probes:
                break
    # With columns in discovery order, every later reduced row has zeros at
    # all earlier pivots. The product of the discovery pivots is the determinant.
    record = dict(prime=prime, points=probes, first_point_seed=start, rank=len(selected),
                  recipes=selected, pivot_points=pivots, determinant=determinant)
    if check:
        assert np.array_equal(frame[:,0],sample(start,prime,integral=True))
        witness = next(i for i in selected if evaluate(i,frame[:,:1],prime)[0])
        direct = direct_control(witness,start,prime)
        assert direct == evaluate(witness,frame[:,:1],prime)[0]
        record["direct_integer_control"] = dict(recipe=witness, point=0, value=direct)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=2000)
    parser.add_argument("--points", type=int, default=160)
    parser.add_argument("--prime", type=int, default=1009)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    ref = json.loads(args.replay.read_bytes()) if args.replay else None
    result = calculate(args.count,ref["points"] if ref else args.points,args.prime,args.workers,
                       ref["recipes"] if ref else None,args.check,args.checkpoint)
    if ref and ref["prime"] == args.prime:
        for key, value in result.items():
            if key in ref:
                assert value == ref[key], key
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["prime"],result["rank"],result["determinant"],flush=True)
