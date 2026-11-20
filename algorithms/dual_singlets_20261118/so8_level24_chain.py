"""Construct and replay rational corrections for the level-24 defining-trace BMN class."""
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat

import mixed_incoming as m
from so8_level24_generator import even_class, small_rational
from trace_descent import action, evaluate, settings


def product(coefficients, matrix):
    """Exact product of an integer row with integer evaluation coordinates."""
    terms = [(i, int(n)) for i, n in enumerate(coefficients) if n]
    if not terms:
        return [0] * matrix.shape[1]
    indices, numbers = zip(*terms)
    selected = matrix[np.asarray(indices)]
    limit = max(abs(int(selected.min())), abs(int(selected.max())))
    if sum(map(abs, numbers)) * max(1, limit) < 2**63:
        return (np.asarray(numbers, dtype=np.int64) @ selected).tolist()
    return list(np.asarray(numbers, dtype=object) @ selected.astype(object))


def clear_denominators(values):
    values = [Fraction(v) for v in values]
    denominator = math.lcm(*(v.denominator for v in values))
    return [int(v * denominator) for v in values], denominator


def solve_image(matrix, target, prime=m.P):
    """Use modular pivots and accept a solution only after every integer equality."""
    wanted, scale = clear_denominators(target)
    if not any(wanted):
        return {}, dict(pivot_rank=None, reconstruction_digits=0)
    field = nmod_mat((matrix % prime).tolist(), prime)
    columns = m.pivotcols(field)
    rows = m.pivotcols(nmod_mat((matrix[:, columns].T % prime).tolist(), prime))
    assert len(rows) == len(columns) and columns
    sub = matrix[np.ix_(rows, columns)]
    inverse = nmod_mat((sub.T % prime).tolist(), prime).inv()
    pivot_target = [wanted[j] for j in columns]
    residual = pivot_target[:]
    modulus = 1
    residues = [0] * len(rows)
    def verified(nums, denominator):
        if product(nums, matrix[rows]) != [n * denominator for n in wanted]:
            return None
        return {i: Fraction(n, denominator * scale) for i, n in zip(rows, nums) if n}
    for step in range(1, 13):
        digit = inverse * nmod_mat(len(rows), 1, [v % prime for v in residual], prime)
        values = list(map(int, digit.entries()))
        residues = [n + modulus * v for n, v in zip(residues, values)]
        modulus *= prime
        difference = [n - v for n, v in zip(residual, product(values, sub))]
        assert all(v % prime == 0 for v in difference)
        residual = [v // prime for v in difference]
        fractions = [small_rational(v, modulus) for v in residues]
        if all(v is not None for v in fractions):
            nums, denominator = clear_denominators(fractions)
            answer = verified(nums, denominator)
            if answer is not None:
                return answer, dict(pivot_rank=len(rows), reconstruction_digits=step)
    exact = fmpz_mat(sub.T.tolist()).solve(fmpz_mat([[v] for v in pivot_target]))
    nums, denominator = clear_denominators([Fraction(str(v)) for v in exact.entries()])
    answer = verified(nums, denominator)
    if answer is None:
        raise ValueError("The selected classical rows do not solve the full rational equation.")
    return answer, dict(pivot_rank=len(rows), reconstruction_digits=None)


def first_image(operator, frame):
    op = even_class(operator)
    quantum = action(op, 1, 8, "SO")
    assert len(op) == 13449 and quantum
    # At deficit two, three, four, or more inputs cannot supply three
    # different theta labels for each fan vertex.
    assert all(sum(sorted((len(l[2]) for tr in w for l in tr),reverse=True)[:3]) < 6 for w in op)
    denominators = math.lcm(*(Fraction(v).denominator for v in quantum.values()))
    row = {w: int(v * denominators) for w, v in quantum.items()}
    coordinate_frames = {}
    for degree, data in frame["target"].items():
        alphabet = [(i,j,tuple(t)) for i,j,t in data["alphabet"]]
        coordinate_frames[int(degree)] = dict(data, alphabet=alphabet)
    values = evaluate([row], coordinate_frames, 8, "SO")[0]
    return [Fraction(int(v), denominators) for v in values], len(quantum)


def run(operator, band, grades, out, replay=None):
    settings()
    data = json.loads(band.read_bytes())
    assert (data["group"], data["colours"], data["flavours"], data["spin"],
            data["target_deficit"], data["maximum_direct_order"]) == ("SO", 8, [4,4,4], 0, 1, 1)
    out.mkdir(exist_ok=True, parents=True)
    previous = {}
    receipt = dict(operator_sha256=hashlib.sha256(operator.read_bytes()).hexdigest(),
                   band_sha256=hashlib.sha256(band.read_bytes()).hexdigest(), corrections=[])
    for grade in range(1,5):
        with np.load(grades/f"grade{grade}.npz", allow_pickle=False) as archive:
            blocks = {}
            for source in (grade-1, grade):
                block = archive[f"{source}_{grade}"]
                assert block.dtype == np.int64 and block.shape == (
                    data["source_dimensions"][source], data["target_dimensions"][grade])
                assert hashlib.sha256(m.compact(block.tolist())).hexdigest() == data["block_hashes"][f"{source}-{grade}"]
                blocks[source] = block
        if grade == 1:
            rhs, receipt["first_action_terms"] = first_image(operator, data["frames"][grade])
        else:
            coefficients = [previous.get(i, 0) for i in range(data["source_dimensions"][grade-1])]
            nums, den = clear_denominators(coefficients)
            rhs = [Fraction(v, den) for v in product(nums, blocks[grade-1])]
        if replay is None:
            previous, info = solve_image(blocks[grade], [-v for v in rhs])
        else:
            saved = replay["corrections"][grade-1]
            assert saved["grade"] == grade
            previous = {i:Fraction(n) for i,n in saved["classical_rows"]}
            info = {k:saved[k] for k in ("pivot_rank", "reconstruction_digits")}
        nums, den = clear_denominators([previous.get(i,0) for i in range(blocks[grade].shape[0])])
        check = product(nums, blocks[grade])
        assert all(Fraction(v, den) + x == 0 for v, x in zip(check, rhs, strict=True))
        record = dict(grade=grade, rhs_nonzero=sum(bool(x) for x in rhs),
                      coordinates_checked=len(check), classical_rows=[[i,str(n)] for i,n in sorted(previous.items())],
                      **info)
        receipt["corrections"].append(record)
        m.write(out/"chain.json", receipt)
        print("grade", grade, "equations", len(check), "nonzero parents", len(previous), flush=True)
    assert len(receipt["corrections"]) == 4
    if replay is not None:
        assert replay == receipt
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("operator","band","grades","out"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--replay", type=Path)
    args = parser.parse_args()
    assert __debug__
    run(args.operator,args.band,args.grades,args.out,
        json.loads(args.replay.read_bytes()) if args.replay else None)
