"""Pair native F4 trace identities, after glueball multiplication, with the corrected row."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import json
from pathlib import Path

from flint import fmpq_mat, fmpz_mat

from flag_boundaries import HERE
from scalar_orthogonal_walls import native, reduce_native
from f4_subgroup_equation import integral_rows
from single_trace_flag import native26
from source_symmetry import centralizer, character, partitions, singlet_trace
from trace_ideal import monomials, order
from trace_source_cycles import GLUE, WEIGHT, quartic_pivots, reduce_tree


def terminal_dimension(q):
    return sum(Fraction(character((2,) * q, g) * character((2, 2, 2), b)
                        * singlet_trace(tuple(sorted(g + b))), centralizer(g) * centralizer(b))
               for g in partitions(2 * q) for b in partitions(6))


def full_reduction(row, basis, index, extra):
    out = Counter()
    for w, x in row.items():
        s, v = order((GLUE,) * extra + w)
        if s:
            out[index[v]] += s * x
    for j in sorted(basis):
        x = out.get(j, 0)
        if x:
            for k, v in basis[j].items():
                out[k] -= x * v
    return {k: v for k, v in out.items() if v}


def obstruction_basis(basis, index):
    src = json.loads(gzip.decompress((HERE/"trace_source_cycles.json.gz").read_bytes()))
    score = json.loads((HERE/"trace_score_descent.json").read_bytes())
    mon = [tuple(map(tuple, m)) for m in src["monomials"]]
    rows = []
    for j in score["singlet_basis_columns"]:
        vec = Counter()
        for pos, x in src["singlet_columns"][j]:
            for k, y in reduce_tree(mon[pos], basis, index).items():
                vec[k] += x * y
        rows.append({k: v for k, v in vec.items() if v})
    return rows, list(map(int, score["eta_Q1_values_at_selected_columns"]))


def calculate(q, workers, maximum=None, sample_terms=4, seed=850775):
    from sympy import Matrix

    weight = (q, q, 2, 2, 2)
    if maximum is None:
        shapes = native.words_by_shape(sum(weight), weight)
    else:
        shapes = defaultdict(set)
        for content in native.contents(weight, sum(weight)):
            for word in native.multitraces_with_content(content, maxlen=maximum):
                sign, canonical = native.so_canon(word)
                if sign:
                    shapes[native.shape(canonical, "spin")].add(canonical)
    stable = [row for words in shapes.values() for row in native.singlets(sorted(words))]
    dim = terminal_dimension(q)
    assert dim.denominator == 1
    ev, frame = integral_rows(stable, "F4", 2*q, min(int(dim), len(stable)), seed,
                              workers, basis=native26(), sample_terms=sample_terms)
    native_rank = ev.rank()
    result = dict(dotted_degree=q, removed_glueball_factors=4-q, stable_singlets=len(stable),
                  maximum_trace_length=maximum, terms_per_sampled_letter=sample_terms,
                  native_character=int(dim), trace_evaluation_rank=native_rank, frame=frame)
    print(q, len(stable), int(dim), native_rank, flush=True)
    basis, index = quartic_pivots()
    differentials, scores = obstruction_basis(basis, index)
    projected = [full_reduction(reduce_native(r), basis, index, 4-q) for r in stable]
    keys = sorted(set().union(*(r.keys() for r in differentials+projected)))
    dest = fmpq_mat([[str(r.get(k, 0)) for k in keys] for r in projected])
    # A genuine identity must lie in this evaluated kernel, even on an incomplete frame.
    ker_columns, nullity = ev.transpose().nullspace()
    if nullity:
        coeff = fmpz_mat([[ker_columns[i, j] for i in range(len(stable))] for j in range(nullity)])
    else:
        coeff = fmpz_mat(0, len(stable))
    assert (coeff * ev).is_zero() and coeff.rank() == nullity
    rel = fmpq_mat(coeff) * dest
    rr, rank = rel.rref()
    compact = Matrix([[Fraction(str(x)) for x in r] for r in rr.tolist()[:rank]]) if rank else Matrix.zeros(0, len(keys))
    five = Matrix([[r.get(k,0) for k in keys] for r in differentials])
    assert five.rank() == 5
    # Coefficients of native identities precede coefficients of the five source images.
    intersection = compact.col_join(-five).T.nullspace()
    coord = [list(z[rank:, 0]) for z in intersection]
    pairings = [sum(x*y for x,y in zip(z,scores)) for z in coord]
    words = monomials(WEIGHT)
    result.update(evaluated_nullity=nullity,native_frame_complete=(native_rank==int(dim)),
                  scaled_projected_relation_rank=rank, terminal_word_coordinates=[words[k] for k in keys] if rank else [],
                  projected_relation_echelon=[[str(x) for x in r] for r in compact.tolist()],
                  intersection_dimension=len(intersection),
                  source_coordinates=[[str(x) for x in z] for z in coord],
                  eta_Q1_pairings=list(map(str,pairings)))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dotted", type=int, choices=(0,1,2,3), required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-trace", type=int)
    parser.add_argument("--sample-terms", type=int, default=4)
    parser.add_argument("--seed", type=int, default=850775)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1 or not 1 <= args.sample_terms <= 52:
        parser.error("assertions and a positive worker count are required")
    data = calculate(args.dotted, args.workers, args.max_trace, args.sample_terms, args.seed)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(data["scaled_projected_relation_rank"],data["intersection_dimension"],data["eta_Q1_pairings"])
