"""Recover the entire cubic response in the eight traceless SU(3) colour directions."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from hashlib import sha256
from math import prod
import gzip
import json

from source_fluctuations import SOURCE, response

HERE = Path(__file__).resolve().parent


def linear(*terms):
    result = Counter()
    for factor, row in terms:
        result.update({key: factor * value for key, value in row.items()})
    return {key: value for key, value in result.items() if value}


def direction(coords):
    matrix = [
        list(coords[:3]),
        list(coords[3:6]),
        [*coords[6:], -coords[0] - coords[4]],
    ]
    assert len(coords) == 8 and not sum(matrix[i][i] for i in range(3))
    return matrix


def source():
    raw = gzip.decompress(SOURCE.read_bytes())
    assert (
        sha256(raw).hexdigest()
        == "ca59e95ecf2ca092722cd612a4cbf51d562f26f06e58a6a7c9d362cdf126ed8a"
    )
    return json.loads(raw)["terms"]


def rebuild():
    terms = source()
    table = {}
    for indices in combinations_with_replacement(range(8), 3):
        coordinates = [indices.count(i) for i in range(8)]
        values = response(terms, direction(coordinates))
        table[indices] = {m: v for m, v in values.items() if m.bit_count() == 9}
    coefficients = {
        (i, i, i): linear((Fraction(1, 27), table[i, i, i])) for i in range(8)
    }
    for i, j in combinations(range(8), 2):
        a, b = coefficients[i, i, i], coefficients[j, j, j]
        u = linear((1, table[i, i, j]), (-8, a), (-1, b))
        v = linear((1, table[i, j, j]), (-1, a), (-8, b))
        coefficients[i, i, j] = linear((Fraction(1, 3), u), (Fraction(-1, 6), v))
        coefficients[i, j, j] = linear((Fraction(-1, 6), u), (Fraction(1, 3), v))
    for key in combinations(range(8), 3):
        lower = [(-1, v) for idx, v in coefficients.items() if set(idx) <= set(key)]
        coefficients[key] = linear((1, table[key]), *lower)
    return {idx: v for idx, v in coefficients.items() if v}


def load():
    raw = json.loads((HERE / "cubic_colours.json").read_text())
    data = {
        tuple(k): {int(m, 16): Fraction(v) for m, v in row.items()} for k, row in raw
    }
    assert len(data) == 22 and sum(map(len, data.values())) == 85
    assert len({m for row in data.values() for m in row}) == 46
    return data


def independent_checks(data):
    # Literal placement of all three insertions, evaluated outside the interpolation grid.
    terms = source()
    for coord in ((1, 2, 3, -1, -3, 1, -2, 4), (-2, 1, 4, 2, 3, -3, 1, -1)):
        expected = linear(
            *((prod(coord[i] for i in key), v) for key, v in data.items())
        )
        assert expected == response(terms, direction(coord), split=True)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    data = load()
    assert rebuild() == data
    independent_checks(data)
    print(
        "120 colour evaluations and two independent slot-by-slot polarizations agree."
    )
