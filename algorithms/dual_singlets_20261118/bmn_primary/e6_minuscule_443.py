"""Exact Cartan evaluations and all possible lower odd factors of the E6 scalar."""
import argparse
from fractions import Fraction
from itertools import permutations, product
import json
from math import lcm, prod
from pathlib import Path

from exceptional_scalar_pairs import colour_data, inverse_killing_form
from exceptional_weyl_pairs import DEGREES


HERE = Path(__file__).resolve().parent


def schur_trace(shape, complete):
    total = 0
    for perm in permutations(range(3)):
        sign = (-1)**sum(i>k for j,i in enumerate(perm) for k in perm[j+1:])
        degrees = [n-i+perm[i] for i,n in enumerate(shape)]
        if min(degrees) >= 0:
            total += sign*prod(complete[n] for n in degrees)
    return total


def factor_multiplicities():
    classes = json.loads((HERE/"exceptional_weyl_classes.json").read_bytes())["groups"]["E6"]
    order = prod(DEGREES["E6"])
    assert sum(size for _,size in classes) == order
    targets = [(i,j,k) for i in range(5) for j in range(i+1)
               for k in range(min(j,3)+1) if (i+j+k)%2]
    numerators = dict.fromkeys(targets, 0)
    for denominator, size in classes:
        complete = [1]
        for k in range(1,7):
            complete.append(-sum(denominator[i]*complete[k-i]
                                 for i in range(1,min(k+1,len(denominator)))))
        for shape in targets:
            numerators[shape] += size*schur_trace(shape,complete)
    assert all(v%order==0 for v in numerators.values())
    assert [w for w,v in numerators.items() if v] == [(4,4,3)]
    assert numerators[4,4,3] == order and len(targets)==16
    return [[*w,v//order] for w,v in numerators.items()]


def evaluated_weights(weights, point):
    return [tuple(sum(a*b for a,b in zip(w,axis)) for axis in point) for w in weights]


def determinant(first, second, last):
    return sum(first[i]*(second[(i+1)%3]*last[(i+2)%3]-second[(i+2)%3]*last[(i+1)%3])
               for i in range(3))


def gram(weights):
    return [[sum(w[i]*w[j] for w in weights) for j in range(6)] for i in range(6)]


def compute():
    roots, _, representations = colour_data("E6")
    weights = representations[5]
    form = inverse_killing_form(roots)
    assert gram(roots)==[[4*n for n in row] for row in gram(weights)]
    dual_weights = [[sum(form[i][j]*w[j] for j in range(6)) for i in range(6)] for w in weights]
    inner = [[sum(x*y for x,y in zip(u,v)) for v in dual_weights] for u in weights]
    denominator = lcm(*(n.denominator for row in inner for n in row))
    coefficients = [[int(n*denominator) for n in row] for row in inner]
    assert denominator==72
    assert all(set(row)=={-2,1,4} for row in coefficients)
    points = (((1,0,0,0,0,1),(0,0,1,0,0,0),(0,0,0,1,0,0)),
              ((-3,5,-1,3,-2,2),(-5,-2,4,-2,0,4),(-1,-1,0,-3,-5,-5)))
    records = []
    for point in points:
        evaluated = evaluated_weights(weights,point)
        root_values = evaluated_weights(roots,point)
        total = 0
        for i,j,k in product(range(27),repeat=3):
            first, middle, last = evaluated[i], evaluated[j], evaluated[k]
            pair = first[0]*last[1]-first[1]*last[0]
            total += coefficients[i][j]*coefficients[j][k]*determinant(first,middle,last)**3*pair
        product_value = sum(u[0]**3*determinant(u,v,w)**2
                            for u in evaluated for v in root_values for w in root_values)
        records.append(dict(simple_coroot_coordinates=point, split_casimir_numerator=total,
                            split_casimir_divisor=denominator**2,
                            scalar_443=str(Fraction(total,denominator**2)), product_522=product_value))
    assert [r["scalar_443"] for r in records] == ["1/12","69349822/3"]
    assert records[1]["product_522"] == -410705095680
    return dict(group="E6",weight=[4,4,3],level=22,outer_parity=-1,
                number_of_minuscule_weights=len(weights),weyl_group_order=prod(DEGREES["E6"]),
                all_odd_subshapes=factor_multiplicities(), evaluations=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(compute(),indent=2)+"\n")
