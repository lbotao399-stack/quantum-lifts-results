"""List lower BMN products and single raising actions reaching weight (10,4,4), Y=13."""
from collections import Counter
from itertools import product
import gzip
import json
from pathlib import Path


def remove_horizontal(shape, boxes):
    """Pieri with three rows. Consecutive rows of the two partitions interlace."""
    ranges = (range(lo, hi+1) for lo, hi in zip((*shape[1:], 0), shape))
    return Counter({mu:1 for mu in product(*ranges) if sum(shape)-sum(mu)==boxes})


def candidates():
    root = Path(__file__).resolve().parent
    public = json.loads((root/"public_so5_j28.json").read_bytes())
    light = [(r["j"], r["lam"], r["h"]) for r in public["columns"] if r["j"] <= 6 and r["h"]]
    assert sorted(light)==[(4, [2,0,0], 1), (6, [2,1,0], 1)]
    odd = json.loads(gzip.decompress((root/"odd_through_j28.json.gz").read_bytes()))
    assert odd["rational_odd_cohomology_through_level_28"]==0
    target = (10,4,4)
    adjacency = remove_horizontal(target, 1)
    quadratic = remove_horizontal(target, 2)
    # s_(2,1) = s_(2) s_(1) - s_(3) for the degree-six [1,1] graviton.
    adjoint = Counter()
    for weight, multiplicity in adjacency.items():
        for parent, n in remove_horizontal(weight, 2).items():
            adjoint[parent] += multiplicity*n
    adjoint.subtract(remove_horizontal(target, 3))
    assert min(adjoint.values()) >= 0
    families = [("degree-four product",quadratic,11), ("degree-six product",+adjoint,11),
                ("one BMN raising supercharge",adjacency,13)]
    return [dict(operation=operation, weight=list(w), level=2*sum(w), length=y, multiplicity=n)
            for operation, weights, y in families for w,n in sorted(weights.items())]


if __name__ == "__main__":
    print(json.dumps(candidates(),indent=2))
