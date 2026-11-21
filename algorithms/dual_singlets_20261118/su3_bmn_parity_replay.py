"""Rebuild the selected BMN products on rational traceless diagonal matrices."""
import argparse
import json
from pathlib import Path

from rank_two_bmn_product_replay import highest, witness
from sp3_omega_products import weight as factor_weight
from su3_bmn_parity_products import labelled
from trace_descent import dimensions


def replay(path):
    record = json.loads(path.read_bytes())
    weight = tuple(record["weight"])
    assert record["group"] == "SU(3)" and record["level"] == 2*sum(weight)
    rows = record["rows"]
    assert [r["length"] for r in rows] == list(range(max(weight),sum(weight)+1))
    euler, protected, minors = [0,0], [0,0], {}
    for row in rows:
        length = row["length"]
        dimensions_by_colour = [sum(dimensions((0,0,*weight),length,3,group).values())
                                for group in ("SU+","SU-")]
        assert dimensions_by_colour == row["even_odd_cochains"]
        for i, dimension in enumerate(dimensions_by_colour):
            euler[i] += (-1)**length * dimension
        products = [tuple(map(tuple,recipe)) for recipe in row["inputs"]]
        assert len(products) == row["rational_product_rank"]
        for recipe in products:
            assert sum(sum(g[:5]) for g in recipe) == length
            assert tuple(map(sum,zip(*(factor_weight(g) for g in recipe)))) == (0,0,*weight)
        if products:
            # First reproduce the saved columns and determinant in the cyclic eigenbasis.
            cyclic = witness([highest(labelled(recipe),weight) for recipe in products])
            assert json.loads(json.dumps(cyclic)) == row["minor"]
            # X=diag(x,y,-x-y) now gives coefficients over Q without a root of unity.
            ordinary = witness([highest(labelled(recipe,True),weight) for recipe in products])
            assert ordinary["rank"] == len(products)
            minors[length] = ordinary
        else:
            assert row["minor"] is None
        protected[length%2] += len(products)
    lower = [max(0,protected[0]-euler[0]),max(0,protected[1]+euler[1])]
    assert euler == record["euler_with_length_sign_by_C"]
    assert protected == record["independent_products_by_C"]
    assert lower == record["opposite_graviton_parity_lower_by_C"]
    assert sum(lower) == record["opposite_graviton_parity_total_lower"]
    return dict(weight=weight,ordinary_diagonal_minors=minors,
                euler_by_C=euler,protected_by_C=protected,opposite_graviton_parity_lower_by_C=lower)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("path",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        raise RuntimeError("verification requires enabled assertions")
    record=replay(a.path)
    a.out.write_text(json.dumps(record,indent=2)+"\n")
