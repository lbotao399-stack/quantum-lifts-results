"""Choose physical Taylor polarizations that detect several Levi source diagrams."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import cache
from itertools import product
import json
from pathlib import Path
import subprocess

from levi_source_minors import branches, taylor_weight


@cache
def remove_vertical(shape):
    strips = []
    for removed in product((0, 1), repeat=len(shape)):
        shorter = tuple(r - k for r, k in zip(shape, removed))
        if any(a < b for a, b in zip(shorter, shorter[1:])):
            continue
        strips.append((tuple(i for i in shorter if i), sum(removed)))
    return strips


@cache
def components(shape):
    """Skew Schur restriction to three even and two odd directions."""
    weights = Counter()
    for middle, negative in remove_vertical(shape):
        for scalar, positive in remove_vertical(middle):
            if len(scalar) <= 3:
                scalar = scalar + (0,) * (3 - len(scalar))
                weights[scalar + (positive, negative)] += 1
    result = {}
    for sector, dimension in weights.items():
        *scalar, positive, negative = sector
        if positive < negative:
            continue
        multiplicity = dimension - weights[(*scalar, positive + 1, negative - 1)]
        assert multiplicity >= 0
        if multiplicity:
            result[sector] = multiplicity
    assert result[taylor_weight(shape)] == 1
    return result


def dimensions(engine, group, sector):
    degree = sum(sector)
    command = [str(engine), "dims", group, *map(str, sector), str(degree - 2), str(degree)]
    data = json.loads(subprocess.check_output(command))
    assert data["gauge"] == group and tuple(data["sector"]) == sector
    return [data["d_Y"].get(str(length), 0) for length in range(degree - 2, degree + 1)]


def calculate(engine, groups, workers):
    out = {}
    for group in groups:
        shapes = branches(1 if group in ("SO7", "Sp3") else 2)
        weights = sorted({w for shape in shapes for w in components(shape)})
        with ThreadPoolExecutor(workers) as pool:
            sizes = dict(zip(weights, pool.map(lambda w: dimensions(engine, group, w), weights)))
        rows = []
        for shape in shapes:
            def cost(w):
                parent, source, target = sizes[w]
                return (source * (parent + target + source), source, w)
            best = min(components(shape), key=cost)
            direct = taylor_weight(shape)
            rows.append(dict(shape=shape, usual=direct, usual_sizes=sizes[direct],
                             selected=best, selected_sizes=sizes[best],
                             branching_multiplicity=components(shape)[best]))
        out[group] = dict(diagrams=rows, dimensions=[dict(weight=w, sizes=sizes[w]) for w in weights])
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--groups", nargs="+", choices=("SU3", "Sp2", "SO7", "Sp3"), default=("SO7", "Sp3"))
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("enable assertions and supply a positive worker count")
    result = calculate(args.engine.resolve(), args.groups, args.workers)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
