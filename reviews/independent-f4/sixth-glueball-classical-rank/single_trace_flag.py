"""Span the degree-fourteen dihedral single traces in the F4 source's classical image."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import cache
import gzip
import json
from math import gcd, lcm
from pathlib import Path

import numpy as np

from flag_boundaries import HERE, adjoint, flags, polynomial, trace_words


CONTENT = (4,4,4,2)
HALF = 4**7


def canonical(codes):
    mask = 4**14-1
    rev = np.zeros_like(codes)
    for i in range(14):
        rev = (rev << 2) + ((codes >> (2*i)) & 3)
    low = np.minimum(codes,rev)
    for _ in range(13):
        codes = ((codes << 2) & mask) + (codes >> 26)
        rev = ((rev << 2) & mask) + (rev >> 26)
        low = np.minimum(low,np.minimum(codes,rev))
    return low


@cache
def catalogue():
    words = np.arange(HALF,dtype=np.int32)
    digits = np.stack([(words >> (2*j)) & 3 for j in range(7)])
    weights = np.stack([np.sum(digits == j,axis=0) for j in range(4)]).T
    groups = {}
    for w in np.unique(weights,axis=0):
        if np.all(w <= CONTENT):
            groups[tuple(w)] = words[np.all(weights == w,axis=1)]
    reps = []
    total = 0
    for weight, left in groups.items():
        right = groups[tuple(np.subtract(CONTENT,weight))]
        for i in range(0,len(left),32):
            codes = (left[i:i+32,None]*HALF + right[None,:]).ravel()
            reps.append(codes[codes == canonical(codes)])
            total += len(codes)
    assert total == 3153150
    return np.sort(np.concatenate(reps))


def single(seed):
    digits, signs = flags()
    slots = np.random.Generator(np.random.PCG64(144433+seed)).permutation(14)
    codes = digits[slots[0]]
    for slot in slots[1:]:
        codes = (codes << 2) + digits[slot]
    keys, inverse = np.unique(canonical(codes),return_inverse=True)
    coefficients = np.bincount(inverse,weights=signs).astype(np.int32)
    keep = coefficients != 0
    indices = np.searchsorted(catalogue(),keys[keep])
    assert np.array_equal(catalogue()[indices], keys[keep])
    return indices.astype(np.int32), coefficients[keep]


def character_dimension():
    from source_symmetry import character
    shape = (4,4,4,2)
    rotations = [character(shape,(14//gcd(14,i),)*gcd(14,i)) for i in range(14)]
    mirrors = [character(shape,(1,1)+(2,)*6),character(shape,(2,)*7)]
    dim = (sum(rotations)+7*sum(mirrors))//28
    assert dim == 474 and sum(rotations)+7*sum(mirrors) == 28*dim
    return dict(rotations=rotations, reflections=mirrors, dimension=dim)


def extend(basis,row,prime):
    row = row.astype(np.int64,copy=True) % prime
    for j in sorted(basis):
        if row[j]:
            row = (row-row[j]*basis[j]) % prime
    active = np.flatnonzero(row)
    if not len(active):
        return None
    j = int(active[0])
    value = int(row[j])
    basis[j] = row*pow(value,-1,prime) % prime
    return j,value


def formal(check=None):
    group = character_dimension()
    prime = 1009
    number = len(catalogue())
    slots = np.sort(np.random.Generator(np.random.PCG64(714443)).choice(number,800,replace=False))
    restriction = np.full(number,-1,dtype=np.int32)
    restriction[slots] = np.arange(len(slots))
    basis = {}
    seeds, columns = [], []
    determinant = 1
    for seed in (range(10000) if check is None else check):
        coords, coeff = single(seed)
        pos = restriction[coords]
        keep = pos >= 0
        row = np.zeros(len(slots), dtype=np.int64)
        row[pos[keep]] = coeff[keep]
        pivot = extend(basis,row,prime)
        if pivot:
            j,y = pivot
            seeds.append(seed)
            columns.append(int(slots[j]))
            determinant = determinant*y % prime
            if len(seeds) == group["dimension"]:
                break
    assert len(seeds) == group["dimension"]
    return dict(character=group, necklaces=number, prime=prime, rank=len(seeds), recipes=seeds,
                necklace_columns=columns, determinant=determinant)


@cache
def native26():
    data = json.loads(gzip.decompress((HERE/"representation26.json.gz").read_bytes()))
    t = np.zeros((52,26,26), dtype=np.int64)
    for k,row in enumerate(data["matrices"]):
        for i,j,v in row:
            t[k,i,j] = v
    h = np.zeros((26,26),dtype=np.int64)
    den = lcm(*(Fraction(v).denominator for _,_,v in data["form"]))
    for i,j,v in data["form"]:
        h[i,j] = int(den*Fraction(v))
    assert np.array_equal(h,h.T) and np.count_nonzero(h@t+t.swapaxes(1,2)@h) == 0
    pivots = {}
    for row in h:
        extend(pivots,row,1009)
    assert len(pivots) == 26
    table = json.loads(gzip.decompress((HERE/"inputs.json.gz").read_bytes()))
    metric = np.zeros((52,52),dtype=np.int64)
    for i,j,v in table["metric"]:
        metric[i,j] = v
    assert np.array_equal(np.einsum("aij,bji->ab",t,t),6*metric)
    brackets = adjoint()
    for i in range(52):
        rhs = np.einsum("kj,kab->jab",brackets[i],t)
        assert np.array_equal(t[i]@t-t@t[i],rhs)
    return t


def power_tables(x,prime):
    dim = x.shape[1]
    out = {1:x}
    for n in range(2,8):
        prev = out[n-1].reshape(-1,dim)
        next_power = np.empty((4**n,dim,dim))
        for k in range(4):
            value = prev @ x[k]
            # All entries are nonnegative integers below 2^53 and the
            # quotient is less than 2^53/prime. Multiples divide exactly.
            value -= np.floor(value/prime)*prime
            next_power[k::4] = value.reshape(-1,dim,dim)
        out[n] = next_power
        if n >= 6:
            del out[n-1]
    return out


def trace_values(powers,codes,n,prime):
    if n == 14:
        a = b = powers[7]
    else:
        a,b = powers[n//2], powers[n-n//2]
    result = np.empty(len(codes),dtype=np.int32)
    for j in range(0,len(codes),192):
        keys = codes[j:j+192]
        entries = np.einsum("nij,nji->n",a[keys//len(b)],b[keys%len(b)])
        result[j:j+192] = entries.astype(np.int64) % prime
    return result


def seed_matrix(seeds):
    cols, entries, offsets = [], [], [0]
    for seed in seeds:
        idx, val = single(seed)
        cols.append(idx)
        entries.append(val)
        offsets.append(offsets[-1]+len(idx))
    assert all(b > a for a,b in zip(offsets,offsets[1:]))
    return np.concatenate(entries).astype(np.int64), np.concatenate(cols), np.array(offsets[:-1])


def short_values(powers,terms,prime):
    _,reps,size = trace_words()
    frame = np.zeros(size,dtype=np.int64)
    start = 1
    for n,words in reps.items():
        frame[start:start+len(words)] = trace_values(powers,words,n,prime)
        start += len(words)
    values = []
    for coeff,indices in terms:
        value = np.prod(frame[indices],axis=0)
        values.append(int(coeff @ (value % prime)) % prime)
    return values


def point(number,prime,native,matrix,short_terms,short_recipes,rows):
    coefficients = np.random.Generator(np.random.PCG64(144336+number)).integers(-3,4,size=(4,52))
    output = {}
    for name in native:
        rep = native26() if name == "26" else adjoint()
        size = rep.shape[1]
        x = ((coefficients @ rep.reshape(52,-1)).reshape(4,size,size) % prime).astype(np.float64)
        tables = power_tables(x,prime)
        vector = trace_values(tables,catalogue(),14,prime)
        entries, columns, starts = matrix
        values = np.add.reduceat(entries*vector[columns],starts) % prime
        output[name] = values
        if name == "52" and short_recipes:
            output["short"] = short_values(tables,short_terms,prime)
        del tables
    return np.array([output[family][index] for family,index in rows],dtype=np.int64)


def physical(args):
    prime = args.prime
    shape = json.loads((HERE/"single_trace_formal.json").read_bytes())
    short = json.loads((HERE/"flag_135.json").read_bytes())
    default = ([["short",s] for s in short["recipes"]]
               + [[name,s] for name in ("26","52") for s in shape["recipes"]])
    ref = json.loads(args.replay.read_bytes()) if args.replay else None
    physical_rows = ref["rows"] if ref else default
    points = ref["points"] if ref else list(range(args.points))
    native = [k for k in ("26","52") if any(f == k or (k == "52" and f == "short") for f,_ in physical_rows)]
    seeds = sorted({s for family,s in physical_rows if family != "short"})
    short_recipes = [s for f,s in physical_rows if f == "short"]
    indices = {s:i for i,s in enumerate(seeds)}
    shorter = {s:i for i,s in enumerate(short_recipes)}
    rows = [(f, (shorter if f == "short" else indices)[s]) for f,s in physical_rows]
    if "26" in native:
        native26()
    matrix = seed_matrix(seeds)
    terms = [polynomial(s) for s in short_recipes]
    basis, selected, kept = {}, [], []
    det = 1
    def record():
        return dict(prime=prime, points=kept, rows=selected, rank=len(kept), determinant=det)
    with ThreadPoolExecutor(args.workers) as pool:
        for lo in range(0,len(points),args.workers):
            batch = points[lo:lo+args.workers]
            def evaluate(j):
                return point(j,prime,native,matrix,terms,short_recipes,rows)
            for j, values in zip(batch,pool.map(evaluate,batch)):
                pivot = extend(basis,values,prime)
                if pivot:
                    i,y = pivot
                    selected.append(physical_rows[i])
                    kept.append(j)
                    det = det*y % prime
            if args.checkpoint:
                tmp = args.checkpoint.with_suffix(".part")
                tmp.write_text(json.dumps(record())+"\n")
                tmp.replace(args.checkpoint)
            if lo//args.workers % 4 == 0:
                print("points",lo+len(batch),"rank",len(kept),flush=True)
            assert len(kept) <= 703
    result = record()
    if ref and ref["prime"] == prime:
        assert result == ref
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("formal","physical"))
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--points",type=int,default=704)
    parser.add_argument("--prime",type=int,choices=(1009,1013),default=1009)
    parser.add_argument("--checkpoint",type=Path)
    parser.add_argument("--replay",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    if args.stage == "formal":
        ref = json.loads(args.replay.read_bytes()) if args.replay else None
        result = formal(ref["recipes"] if ref else None)
        assert ref is None or result == ref
    else:
        result = physical(args)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["rank"],result["determinant"],flush=True)
