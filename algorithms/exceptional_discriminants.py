"""Check trace coordinates on the regular Cartans of G2 and F4 in characteristic zero."""
import argparse
from itertools import product
import json
from pathlib import Path

from sympy import Matrix, diff, eye, factor, prod, symbols


def g2():
    x, y = symbols("x y")
    z = -x-y
    weights = (0, x, -x, y, -y, z, -z)
    traces = [sum(weight**n for weight in weights) for n in (2, 6)]
    second = x*y + y*z + z*x
    third = x*y*z
    positive = x*y*z*(x-y)*(x-z)*(y-z)
    jacobian = Matrix([[diff(t, var) for var in (x, y)] for t in traces]).det()
    discriminant = prod((a-b)**2 for i, a in enumerate(weights) for b in weights[i+1:])
    assert factor(traces[0]+4*second) == factor(traces[1]+4*second**3-6*third**2) == 0
    assert factor(jacobian+48*positive) == 0
    assert factor(discriminant-64*third**10*(-4*second**3-27*third**2)**2) == 0
    return dict(representation_dimension=7, degrees=[2, 6], root_linear_factors=6,
                jacobian_over_ordered_root_product=-48,
                jacobian_at_x1_y2=int(jacobian.subs({x: 1, y: 2})))


def f4():
    basis = [eye(4)[:, i] for i in range(4)]
    short = [s*e for e in basis for s in (-1, 1)]
    short.extend(Matrix(signs)/2 for signs in product((-1, 1), repeat=4))
    long = [s*e+t*f for i, e in enumerate(basis) for f in basis[i+1:]
            for s, t in product((-1, 1), repeat=2)]
    roots = short + long
    positive = [w for w in roots if next(v for v in w if v) > 0]
    weights = {tuple(w) for w in short}
    assert len(weights) == len(positive) == 24 and len(set(map(tuple, roots))) == 48
    # All four degrees are positive, so the two zero weights of the 26 drop out.
    degrees = (2, 6, 8, 12)
    assert sum(n-1 for n in degrees) == len(positive)
    for alpha in positive:
        reflection = eye(4)-2*alpha*alpha.T/alpha.dot(alpha)
        assert reflection.det() == -1 and {tuple(reflection*w) for w in short} == weights

    records = []
    for point in ((1, 2, 3, 5), (1, 2, 4, 8)):
        h = Matrix(point)
        rows = [[n*sum(w[i]*w.dot(h)**(n-1) for w in short) for i in range(4)] for n in degrees]
        determinant = Matrix(rows).det()
        root_product = prod(alpha.dot(h) for alpha in positive)
        assert root_product and determinant
        records.append(dict(point=list(point), jacobian=str(determinant), root_product=str(root_product),
                            ratio=str(determinant/root_product)))
    assert records[0]["ratio"] == records[1]["ratio"]
    return dict(representation_dimension=26, degrees=list(degrees), short_weights=24, zero_weights=2,
                positive_roots=24, checks=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    data = {"G2": g2(), "F4": f4()}
    args.out.write_text(json.dumps(data, indent=2)+"\n")
