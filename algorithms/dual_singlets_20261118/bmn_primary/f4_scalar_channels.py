"""Check all lower scalar factors and all full-letter singlet partners of the F4 scalar."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mixed_incoming as m
from bmn_primary_index import cochains
from f4_formal_cycles import probe
from f4_singlet_characters import singlets as dimensions
sys.path.insert(0, str(HERE.parents[1]))
from verify_g2_s3_currents import components, diagonal_probe


def current_product_minor():
    # Octonion derivations act on the seven imaginary units in each of the
    # three off-diagonal entries. The remaining six Albert coordinates are fixed.
    rows = []
    for matter, mask in ((True, 31), (False, 13)):
        fields = {}
        for letter, diagonal in diagonal_probe(matter).items():
            fields[letter] = {(7*copy+i, 7*copy+i): v for copy in range(3)
                              for i, v in enumerate(diagonal) if v}
        values = [v.get(mask, 0) for v in components(fields, -1)]
        if any(x != int(x) for x in values):
            raise ArithmeticError("Unexpected denominator on the Albert Cartan jets")
        rows.append(list(map(int, values)))
    determinant = rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]
    assert rows == [[-1728,0],[432,15552]] and determinant == -26873856
    return dict(rows=rows, determinant=determinant)


def image(weight, length, colour_dimensions, lib, prime, seed):
    blocks = m.words_by_shape(length, weight)
    images = [m.q(x) for _, block in sorted(blocks.items()) for x in m.singlets(block)]
    if any(m.q(x) for x in images):
        raise ArithmeticError("The classical supercharge did not square to zero")
    words, rows = m.flatten(images)
    if not words:
        return {"rank": 0, "frames": []}
    by_degree = {m.nferm(w): colour_dimensions.get(m.nferm(w), 0) for w in words}
    evaluated, frames = probe(rows, words, by_degree, prime, lib, 2, seed)
    result = m.minor(evaluated, prime)
    result["frames"] = frames
    return result


def calculate(work, prime):
    lib = m.library(work, 27, prime)
    shapes, levels = cochains(9, "F4")
    factors = []
    for shape in shapes:
        degree = sum(shape)
        if shape[0] > 3 or degree >= 9 or degree % 2 == 0:
            continue
        size = levels[degree].get(shape, [0])[0]
        weight = (0, 0, *shape)
        certificate = image(weight, degree-1, {0: size}, lib, prime,
                            710000+100*shape[0]+10*shape[1]+shape[2]) if size else {"rank": 0}
        factors.append(dict(weight=shape, degree=degree, dimension=size, incoming=certificate))
        if certificate["rank"] != size:
            raise ArithmeticError(f"Unresolved lower scalar factor {shape}")
    partners = []
    for t in range(4):
        weight = (t, t) + (3-t,)*3
        length = 8-t
        size = dimensions("F4", weight, length)
        assert m.sector(weight, length) == (1, 1, 6, 6, 6)
        assert sum(weight) == length+1
        before = image(weight, length-1, size, lib, prime, 1046001+7*t)
        after = image(weight, length, dimensions("F4", weight, length+1), lib, prime, 2875900+7*t)
        upper = sum(size.values())-before["rank"]-after["rank"]
        if upper < 0:
            raise ArithmeticError("Incompatible consecutive image ranks")
        partners.append(dict(weight=weight, length=length, dimensions_by_fermions=size,
                             incoming=before, outgoing=after, cohomology_upper=upper))
    return dict(group="F4", prime=prime, lower_scalar_factors=factors,
                physical_incoming_doubled_charges=[1,1,6,6,6], source_blocks=partners,
                exact_current_product_minor=current_product_minor())


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=1000003)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args=parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    data=calculate(args.work, args.prime)
    args.out.write_text(json.dumps(data, indent=2)+"\n")
    print("Physical source upper bounds", [r["cohomology_upper"] for r in data["source_blocks"]])
