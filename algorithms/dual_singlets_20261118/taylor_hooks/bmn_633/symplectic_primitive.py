"""Solve Q0 K = P in the six-dimensional defining representation of USp(6)."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from flint import nmod_mat
import mixed_incoming as m
from dual633_pairing import integral
from dual_charge_bands import coordinates
from product_replay import reconstruct
from so8_level24_chain import solve_image
from terminal_interval_ranks import sample
from trace_descent import highest_weights, multiply_rows


HERE = Path(__file__).resolve().parent
PAGE = HERE.parents[1] / "sp3_omega/lowered_source.json.gz"
WEIGHT = (0, 0, 6, 3, 3)


def source_and_target():
    page = json.loads(gzip.decompress(PAGE.read_bytes()))
    assert (page["weight"], page["length"], page["group"], page["colours"], page["prime"]) == (
        list(WEIGHT), 9, "Sp", 6, m.P)
    source = highest_weights(WEIGHT, 9)
    basis = [source[int(degree)][i] for degree, frame in page["centre_frames"].items()
             for i in frame["minor_rows"]]
    assert len(basis) == 1861
    pivots = page["outgoing"]["rows"]
    assert len(pivots) == len(set(pivots)) == 1149
    polynomial, _ = reconstruct()
    return page, [basis[i] for i in pivots], polynomial


def modular_solution(page, rows, polynomial, library, workers):
    images = [m.q(row) for row in rows] + [polynomial]
    alphabet = sorted({l for content in m.contents(WEIGHT, 10) for l in content})
    blocks = [coordinates(images, int(degree), alphabet, frame["samples"], 6, "Sp", library, workers)
              for degree, frame in page["following_frames"].items()]
    table = np.column_stack(blocks)
    indices = page["outgoing"]["columns"]
    square = table[:-1, indices]
    assert hashlib.sha256(square.astype("<u4").tobytes()).hexdigest() == page["outgoing"]["minor_sha256"]
    matrix = nmod_mat(square.tolist(), m.P)
    assert int(matrix.det()) == page["outgoing"]["determinant"] == 995218
    solution = matrix.transpose().solve(nmod_mat([[int(v)] for v in table[-1, indices]], m.P))
    values = np.asarray(list(map(int, solution.entries())), dtype=np.int64)
    assert len(values) * (m.P-1)**2 < 2**63
    assert np.array_equal(values @ table[:-1] % m.P, table[-1])
    return values


def signed_evaluations(rows, batch, seed):
    """One set of integral evaluation functionals in each of the two target degrees."""
    pieces = []
    alphabet = sorted({l for content in m.contents(WEIGHT, 10) for l in content})
    for degree in (0, 2):
        restricted = [{w: c for w, c in row.items() if m.nferm(w) == degree} for row in rows]
        words = sorted({w for row in restricted for w in row})
        if words:
            values = sample(words, degree, 6, batch, seed + 317*degree, "Sp", alphabet)
            pieces.append(multiply_rows(restricted, words, values))
    return np.column_stack(pieces)


def reconstruct_primitive(page, basis, polynomial, residues, out):
    indices = list(map(int, np.flatnonzero(residues)))
    images = [m.q(basis[i]) for i in indices]
    scaled, denominator = integral(polynomial)
    blocks = []
    width = 0
    rank = 0
    for i in range(512):
        block = signed_evaluations(images + [scaled], 24, 450107 + 4099*i)
        blocks.append(block)
        width += block.shape[1]
        if width < len(images) or i % 4:
            continue
        table = np.column_stack(blocks)
        rank = nmod_mat((table[:-1] % m.P).tolist(), m.P).rank()
        print("integer probe coordinates", width, "rank", rank, "/", len(images), flush=True)
        if rank == len(images):
            break
    if rank != len(images):
        raise ValueError("signed probes did not separate the selected boundary columns")
    solution, info = solve_image(table[:-1], table[-1].tolist())
    coeff = {indices[i]: c/denominator for i,c in solution.items()}
    assert all(Fraction(c).denominator % m.P and
               c.numerator*pow(c.denominator, -1, m.P) % m.P == int(residues[i])
               for i,c in coeff.items())
    row = Counter()
    for i,c in coeff.items():
        for word, n in basis[i].items():
            row[word] += c*n
    row = {w:c for w,c in row.items() if c}
    formal = Counter(m.q(row))
    formal.subtract(polynomial)
    formal = {w:c for w,c in formal.items() if c}
    receipt = dict(prime=m.P, page_sha256=hashlib.sha256(PAGE.read_bytes()).hexdigest(),
                   modular_minor_determinant=995218, modular_equations=1832,
                   modular_coefficients=len(indices), signed_probe_rank=rank, signed_probe_coordinates=width,
                   coefficients_in_outgoing_pivots=[[i,str(c)] for i,c in sorted(coeff.items())],
                   polynomial_terms=len(polynomial), primitive_terms=len(row),
                   primitive_denominator=integral(row)[1], formal_residual_terms=len(formal), **info)
    m.write(out, receipt)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not __debug__: parser.error("enable assertions")
    args.work.mkdir(parents=True, exist_ok=True)
    page, basis, polynomial = source_and_target()
    library = m.library(args.work, 6, m.P)
    residues = modular_solution(page, basis, polynomial, library, args.workers)
    print("solved across 1832 coordinates with", np.count_nonzero(residues), "nonzero coefficients", flush=True)
    reconstruct_primitive(page, basis, polynomial, residues, args.out)
