"""Regenerate the B row over Z/(1009*1013) with bounded signed-integer contractions."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from hashlib import sha256
from itertools import combinations, permutations
import json
import os
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.sparse import csr_matrix

from verify_crossed import ARCHIVE_SHA256, DIRECTORY, GRAPHS, PIVOTS, parity

MODULUS = 1009 * 1013


def apply(matrix, tensor, axes):
    first = tuple(range(len(axes)))
    aligned = np.moveaxis(tensor, axes, first).reshape(26 ** len(axes), -1)
    value = (matrix @ aligned) % MODULUS
    return np.moveaxis(value.reshape((26,) * 4), first, axes)


@cache
def inputs(archive):
    def get(z, name):
        return json.loads(z.read(f"{DIRECTORY}/{name}"))

    with ZipFile(archive) as z:
        native = get(z, "fixtures/color-basis-native.json")
        cubic = get(z, f"{PIVOTS}/cubic.json")
        points = get(z, f"{PIVOTS}/selected-1009.json")["points"]
    matrices = np.zeros((52, 26, 26), dtype=np.int64)
    for i, entries in enumerate(native["representation26"]):
        for j, k, value in entries:
            matrices[i, j, k] = int(value)
    form = np.zeros((26, 26), dtype=np.int64)
    for i, j, value in native["representationInvariantForm"]:
        form[i, j] = int(value)
    inverse = form.copy()
    inverse[:2, :2] = np.array([[2, 1], [1, 2]]) * pow(3, -1, MODULUS) % MODULUS
    assert np.array_equal(form @ inverse % MODULUS, np.eye(26, dtype=np.int64))
    c = np.zeros((26,) * 3, dtype=np.int64)
    for i, j, k, v in cubic["cubic"]:
        c[i, j, k] = v
    assert cubic["k"] == "28/3"
    factor = 3 * pow(28, -1, MODULUS) % MODULUS

    # inverse and C are multiplied on the second and third indices by batched matmul.
    def corner(x):
        x = x.reshape(26, -1) % MODULUS
        return ((((x.T @ inverse) % MODULUS) @ x % MODULUS) * factor % MODULUS).reshape(
            (26,) * 4
        )

    low = corner(c)
    high = corner((inverse @ c % MODULUS) @ inverse % MODULUS).transpose(0, 2, 1, 3)
    exchange = csr_matrix(
        corner(inverse @ c % MODULUS).transpose(0, 2, 1, 3).reshape(676, 676)
    )
    return matrices, points, low, high, exchange


def half(matrices, brackets, tensor, endpoint, mid):
    total = np.zeros_like(tensor)
    for i, j in combinations(range(4), 2):
        first, second = [s for s in range(4) if s not in (i, j)]
        marked = apply(brackets[i, j], tensor, (endpoint[0],))
        for a, b in [(first, second), (second, first)]:
            t = apply(matrices[a], marked, (endpoint[1],))
            total += parity((i, j, a, b)) * apply(matrices[b], t, (endpoint[2],))
    total %= MODULUS
    result = np.zeros_like(total)
    for order in permutations(range(3)):
        term = total
        for i, axis in zip(order, mid):
            term = apply(matrices[i], term, (axis,))
        result += parity(order) * term
    return result % MODULUS


def evaluate(arg):
    archive, position = arg
    matrices, points, low, high, exchange = inputs(archive)
    coeff = np.random.Generator(np.random.PCG64(144336 + points[position])).integers(
        -3, 4, size=(4, 52)
    )
    x = (coeff @ matrices.reshape(52, -1)).reshape(4, 26, 26) % MODULUS
    brackets = {
        (i, j): (x[i] @ x[j] - x[j] @ x[i]) % MODULUS
        for i, j in combinations(range(4), 2)
    }
    lhs = half(
        [csr_matrix(t.T) for t in x],
        {a: csr_matrix(t.T) for a, t in brackets.items()},
        low,
        (1, 0, 3),
        (0, 1, 3),
    )
    rhs = half(
        [csr_matrix(t) for t in x],
        {a: csr_matrix(t) for a, t in brackets.items()},
        high,
        (0, 2, 3),
        (1, 2, 3),
    )
    for pair in ((1, 2), (0, 1)):
        rhs = apply(exchange, rhs, pair)
    # This also bounds all 26- and 676-term matmuls of reduced residues.
    assert 26**4 * (MODULUS - 1) ** 2 < 2**63
    return int(lhs.ravel() @ rhs.ravel()) % MODULUS


if __name__ == "__main__":
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("archive")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--positions", type=int, nargs="+")
    a = ap.parse_args()
    if (
        not __debug__
        or sha256(Path(a.archive).read_bytes()).hexdigest() != ARCHIVE_SHA256
    ):
        raise SystemExit("use the pinned archive with assertions enabled")
    archive = os.path.abspath(a.archive)
    positions = a.positions if a.positions is not None else range(703)
    saved = {}
    with ZipFile(archive) as z:
        for p in (1009, 1013):
            saved[p] = json.loads(z.read(f"{DIRECTORY}/{GRAPHS}/outgoing-{p}.json"))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        found = pool.map(evaluate, ((archive, i) for i in positions), chunksize=2)
        for i, value in zip(positions, found, strict=True):
            assert all(value % p == saved[p][i][1] for p in saved), (i, value)
    print(json.dumps(dict(points=len(positions), primes=list(saved), modulus=MODULUS)))
