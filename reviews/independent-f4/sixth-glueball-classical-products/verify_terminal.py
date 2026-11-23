"""Rebuild the ordered terminal minor from the two colour frames."""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import json

import numpy as np
from scipy.sparse import csr_matrix

from colour import bareiss_residue, invariant_ends, matrices
from exceptional_colour import projectors
from four_arm import point as strand_point
from open_pairs import Frame
from terminal import PolynomialFrame


HERE = Path(__file__).resolve().parent


def evaluate(prime):
    recipes = json.loads((HERE / "recipes.json").read_bytes())
    shape = (5, 4, 2)
    trace = PolynomialFrame(shape, recipes["trace_seeds"])
    paired = Frame(shape, recipes["paired_word_seeds"])
    frames = matrices(HERE)
    defining = frames[0]
    coordinates, inverse = projectors(HERE, defining, prime)[12]
    low, high, _ = invariant_ends(defining, coordinates, inverse, prime)
    low = low.astype(np.float64)
    high = high.astype(np.float64)
    coordinates = csr_matrix(coordinates.astype(np.float64))
    columns = []
    for seed in recipes["points"]:
        ordinary = trace.evaluate(frames, seed, prime)
        projected = paired.evaluate(defining, coordinates, inverse, seed, prime)
        strands = strand_point(seed, defining, low, high, prime, recipes["strand_layouts"])
        columns.append([*ordinary, *projected, *strands])
    matrix = np.array(columns, dtype=np.int64).T
    assert matrix.shape == (55, 55)
    reference = HERE / f"selected-{prime}.json"
    if reference.exists():
        assert np.array_equal(matrix, json.loads(reference.read_bytes()))
    determinant = bareiss_residue(matrix, prime)
    assert determinant
    return dict(prime=prime, evaluated_rows=[51, 2, 2], dimension=55,
                determinant_residue=determinant,
                ordered_le_u16_sha256=sha256(matrix.astype("<u2").tobytes()).hexdigest())


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, required=True, choices=(1009, 1013))
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable assertions.")
    print(json.dumps(evaluate(args.prime), indent=2))
