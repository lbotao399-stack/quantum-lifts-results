"""Evaluate classical boundaries with a second jet at each end of an Albert network."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache, partial
from itertools import combinations, permutations
from multiprocessing import get_context
from pathlib import Path
import argparse
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "rational_section"))
from integer_section import CROSSED, SIGNS, action, native26, tensors


PAIRS = tuple(combinations(range(4), 2))


def apply(matrix, tensor, axes):
    front = tuple(range(len(axes)))
    array = np.moveaxis(tensor, axes, front).reshape(26 ** len(axes), -1)
    acted = matrix @ array
    return np.moveaxis(acted.reshape((26,) * 4), front, axes)


def balanced(tensor, prime):
    return (tensor + prime // 2) % prime - prime // 2


def half(scalars, marks, corner, end, middle, prime):
    scalar_actions = [action(x)[0] for x in scalars]
    a, b, c = end
    total = np.zeros_like(corner)
    # The three inputs to each term have integral absolute value at most prime/2.
    # Reduce once after the first action to bound all three matrix products.
    for pos, (i, j) in enumerate(PAIRS):
        k, l = (h for h in range(4) if h not in (i, j))
        marked = balanced(apply(action(marks[pos])[0], corner, (a,)), prime)
        forward = apply(scalar_actions[l], apply(scalar_actions[k], marked, (b,)), (c,))
        reverse = apply(scalar_actions[k], apply(scalar_actions[l], marked, (b,)), (c,))
        total += SIGNS[i, j, k, l] * (forward - reverse)
    total = balanced(total, prime)
    combined = np.zeros_like(total)
    for order in permutations(range(3)):
        first, second, third = (scalar_actions[i] for i in order)
        value = balanced(apply(first, total, (middle[0],)), prime)
        combined += SIGNS[order] * apply(third, apply(second, value, (middle[1],)), (middle[2],))
    return balanced(combined, prime)


@cache
def colour_tensors():
    low, high, exchange, scales = tensors()
    assert scales == (28, 252, 84)
    assert max(np.max(np.abs(low)), np.max(np.abs(high))) * 26 * 1013 < 2**53
    assert max(np.abs(exchange).sum(axis=0).max(), np.abs(exchange).sum(axis=1).max()) * 1013 < 2**53
    return low, high, action(exchange)[0], action(exchange.T)[0]


def samples(number, prime, closed=False):
    native = native26().reshape(52, -1)
    rng = np.random.default_rng(14337 + number)
    x = (rng.integers(-3, 4, size=(10, 52)) @ native).reshape(10, 26, 26)
    assert np.max(np.abs(x)) <= prime // 2
    brackets = np.stack([x[i] @ x[j] - x[j] @ x[i] for i, j in PAIRS])
    return x[:4], balanced(brackets, prime), balanced(brackets if closed else x[4:], prime)


def make_halves(scalars, brackets, marks, side, prime, pattern):
    corners = colour_tensors()
    if side == 0:
        scalars, brackets, marks = (a.transpose(0, 2, 1) for a in (scalars, brackets, marks))
    end, middle = pattern[2 * side : 2 * side + 2]
    return [half(scalars, m, corners[side], end, middle, prime) for m in (brackets, marks)]


def features(halves, exchange, prime):
    bank = []
    for arr in halves:
        bank.append(np.stack([arr.ravel(), *(balanced(apply(exchange, arr, pair), prime).ravel()
                                             for pair in PAIRS)]))
    return bank


def pure_column(tensor, matrices, axes, prime):
    total = np.zeros_like(tensor)
    for order in permutations(range(len(axes))):
        term = tensor
        for j, (i, axis) in enumerate(zip(order, axes)):
            term = apply(matrices[i], term, (axis,))
            if j % 2 == 0:
                term = balanced(term, prime)
        total += SIGNS[order] * term
    return balanced(total, prime)


def jacobi_columns(scalars, brackets, marks, prime):
    low, high, exchange, transposed = colour_tensors()
    xs = [action(x.T)[0] for x in scalars]
    ys = {pair: marks[j] for j, pair in enumerate(PAIRS)}
    total = np.zeros_like(low)
    for l in range(4):
        i, j, k = (n for n in range(4) if n != l)
        jac = np.zeros((26, 26), dtype=np.int64)
        for a, b, c in ((i, j, k), (j, k, i), (k, i, j)):
            y = ys[min(b, c), max(b, c)]
            jac += (1 if b < c else -1) * (scalars[a] @ y - y @ scalars[a])
        term = apply(action(balanced(jac.T, prime))[0], low, (0,))
        total += SIGNS[i, j, k, l] * apply(xs[l], term, (2,))
    left = pure_column(balanced(total, prime), xs, (1, 3, 2), prime)
    two = np.zeros_like(low)
    for pos, (i, j) in enumerate(PAIRS):
        k, l = (n for n in range(4) if n not in (i, j))
        tail = PAIRS.index((k, l))
        first = apply(action(brackets[pos].T)[0], low, (0,))
        second = apply(action(marks[pos].T)[0], low, (0,))
        two += SIGNS[i, j, k, l] * (apply(action(marks[tail].T)[0], first, (2,))
                                    - apply(action(brackets[tail].T)[0], second, (2,)))
    left_two = pure_column(balanced(two, prime), xs, (1, 3, 2), prime)
    xs = [action(x)[0] for x in scalars]
    right = pure_column(high, xs, (0, 2, 1, 3), prime)
    right = pure_column(right, xs, (3, 1, 0), prime)
    a, c = features([left, left_two], transposed, prime)
    b = features([right], exchange, prime)[0]
    with np.errstate(all="ignore"):
        values = np.stack([a @ b.T, c @ b.T])
    assert np.array_equal(values, values.astype(np.int64))
    return np.asarray(values.astype(np.int64) % prime, dtype=np.uint16).ravel()


def evaluate(number, prime, closed=False, network="reference"):
    assert prime in (1009, 1013)
    assert 12 * 26**2 * (prime // 2 + 1)**3 < 2**53
    assert 26**4 * prime**2 < 2**53
    scalars, brackets, marks = samples(number, prime, closed)
    if network == "jacobi":
        return jacobi_columns(scalars, brackets, marks, prime)
    _, _, exchange, transposed = colour_tensors()
    pattern = ((0, 2, 1), (0, 0, 2), (2, 3, 0), (1, 3, 2)) if network == "alternate" else CROSSED
    left = features(make_halves(scalars, brackets, marks, 0, prime, pattern), transposed, prime)
    right = features(make_halves(scalars, brackets, marks, 1, prime, pattern), exchange, prime)
    with np.errstate(all="ignore"):
        values = left[0] @ right[1].T - left[1] @ right[0].T
    assert np.all(np.isfinite(values)) and np.array_equal(values, values.astype(np.int64))
    if number in (12, 219):
        for i, j in ((0, 0), (2, 5), (5, 2), (6, 6)):
            literal = sum(sign * int(a[i].astype(np.int64) @ b[j].astype(np.int64))
                          for sign, a, b in ((1, left[0], right[1]), (-1, left[1], right[0])))
            assert literal == values[i, j]
    result = [values]
    if network != "reference":
        for pair in PAIRS:
            middle = [np.stack([balanced(apply(exchange, a.reshape((26,) * 4), pair), prime).ravel()
                                for a in leg]) for leg in right]
            with np.errstate(all="ignore"):
                result.append(left[0] @ middle[1].T - left[1] @ middle[0].T)
    all_values = np.stack(result)
    assert np.all(np.isfinite(all_values)) and np.array_equal(all_values, all_values.astype(np.int64))
    return np.asarray(all_values.astype(np.int64) % prime, dtype=np.uint16).ravel()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--points", type=Path, help="JSON list of integer point IDs")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--network", choices=("reference", "third", "alternate", "jacobi"), default="reference")
    args = parser.parse_args()
    if args.points is None:
        for i in (12, 219):
            row = evaluate(i, args.prime, network=args.network)
            assert row.any() and not evaluate(i, args.prime, True, args.network).any()
            print(i, len(row), int(np.count_nonzero(row)), row[:7].tolist())
        return
    numbers = json.loads(args.points.read_bytes())
    assert args.output and len(numbers) == len(set(numbers))
    width = {"reference": 49, "jacobi": 98, "third": 343, "alternate": 343}[args.network]
    matrix = np.empty((len(numbers), width), dtype=np.uint16)
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        for i, row in enumerate(pool.map(partial(evaluate, prime=args.prime, network=args.network), numbers, chunksize=2)):
            matrix[i] = row
            if (i + 1) % 200 == 0:
                print("evaluated", i + 1, flush=True)
    np.savez_compressed(args.output, points=numbers, values=matrix)


if __name__ == "__main__":
    main()
