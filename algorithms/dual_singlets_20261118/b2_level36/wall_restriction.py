"""Reconstruct the long-root restriction of the published level-36 class."""
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ALG = HERE.parents[1]
sys.path.insert(0, str(ALG / "su3_j34_20261118"))
from match_threshold import paper_core
from verify import load_words, pairing, polynomial


def public_table(group, digest):
    raw = gzip.decompress((HERE / f"public_su{group}.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == digest
    table = json.loads(raw)
    assert table["N"] == group
    return table["columns"]


def spin_five(columns, level):
    selected = [(r, s) for r in columns if r["j"] == level for s in r["sectors"]
                if s["twoJ"] == 5]
    assert selected and all(r["all_checks"] and s["check"] for r, s in selected)
    return dict(columns=len(selected), cohomology=sum(s["h"] for _, s in selected),
                gravitons=sum(s["grav"] for _, s in selected),
                excess=[[r["lam"], s["BH"]] for r, s in selected if s["BH"]])


def exchange_matrix():
    # Basis h^k t^(6-k), with h=hat phi and t=tilde phi from (2.15) and (4.14).
    # Swapping colours 3 and 4 sends (h,t) to ((h+4t)/3,(2h-t)/3).
    columns = []
    for k in range(7):
        column = Counter()
        for i in range(k+1):
            for j in range(7-k):
                column[i+j] += (comb(k,i) * 4**(k-i) * comb(6-k,j)
                                * 2**j * (-1)**(6-k-j))
        columns.append([column[i] for i in range(7)])
    return list(map(list, zip(*columns)))


def null_line(matrix):
    rows = [[Fraction(x) for x in r] for r in matrix]
    pivots = []
    for j in range(7):
        i = len(pivots)
        chosen = next((k for k in range(i,len(rows)) if rows[k][j]), None)
        if chosen is None:
            continue
        rows[i], rows[chosen] = rows[chosen], rows[i]
        rows[i] = [x/rows[i][j] for x in rows[i]]
        for k in range(len(rows)):
            if k != i:
                rows[k] = [x-rows[k][j]*y for x,y in zip(rows[k],rows[i])]
        pivots.append(j)
    free, = set(range(7)).difference(pivots)
    line = [Fraction(0)]*7
    line[free] = 1
    for i,j in enumerate(pivots):
        line[j] = -rows[i][free]
    assert all(sum(x*y for x,y in zip(row,line)) == 0 for row in matrix)
    return pivots, line


def calculate():
    two = public_table(2, "f486b664d48dcff80cc84de9a32b13cd702500069d71d120248ef04edcd69528")
    three = public_table(3, "cc4ec033e2dfecc3f4c7673a91859e6c28da342cc20bfa1eda032b50aef4b2c5")
    odd = []
    for row in two:
        if row["j"] > 24:
            continue
        assert row["all_checks"]
        for sector in row["sectors"]:
            assert sector["check"]
            length = row["j"]//2-sector["twoJ"]
            if length%2 and sector["h"]:
                odd.append([row["j"],length,row["lam"],sector["h"]])
    assert odd == [[24,7,[4,4,4],1]]
    stages = {level:spin_five(three,level) for level in (24,26,28,30)}
    assert all(not stages[j]["excess"] for j in (24,26,28))
    assert stages[30]["excess"] == [[[7,4,4],1]]

    words = load_words(ALG / "su3_independent_prime_20261118/runtime/work"
                       / "cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/bmn/words.py")
    cert = json.loads((ALG/"su3_j34_20261118/certificate.json").read_bytes())["channels"][2]
    ell = {(i,j):Fraction(v) for i,j,v in cert["witness"]["monomial_coefficients"]}
    core = paper_core(words)
    parents = words.words_of_sector((6,4,4,4),3)
    assert len(parents) == 197 and not polynomial(words.Q_expr(core))
    assert all(not pairing(ell,polynomial(words.Q_multi(w))) for w in parents)
    detected = pairing(ell,polynomial(core))
    assert detected == -20

    # h^0,h^1,h^2 correspond to outer powers 6,5,4. Their SU(3)
    # coefficients have C opposite to gravitons and vanish in cohomology.
    exchange = exchange_matrix()
    equations = [[int(i==j) for j in range(7)] for i in range(3)]
    equations += [[x-3**6*int(i==j) for j,x in enumerate(row)]
                  for i,row in enumerate(exchange)]
    pivots,line = null_line(equations)
    assert line == [0,0,0,64,48,12,1]
    # On diag(0,0,z,-z), h=-z/3 and t=z/3. Use a ratio so (4.15)
    # supplies only nonvanishing and the parent normalization cancels.
    wall = sum(Fraction((-1)**k,3**6)*c for k,c in enumerate(line))
    ratio = wall / line[3]
    assert ratio == Fraction(-1,1728)
    return dict(published_input="arXiv:2610.20340v2, Table 5 and Eqs. (2.15), (4.14), (4.15)",
                published_mixed_coefficient="nonzero; absolute value not supplied",
                su3_spin_five=stages, su2_odd_through_24=odd, su2_parents_checked=len(parents),
                su2_core_pairing=int(detected), exchange_numerator=exchange, exchange_denominator=3**6,
                constraint_rank=len(pivots), homogeneous_line=[int(c) for c in line],
                long_wall_relative_to_mixed_coefficient=str(ratio),
                restricted_class={"group":"USp(4)","j":36,"length":13,"weight":[10,4,4],
                                  "flavour_dynkin":[6,0],"BMN_spin":5},
                possible_lower_scalar_dressing={"group":"USp(4)","j":32,"length":11,
                                                "weight":[8,4,4],"flavour_dynkin":[4,0]},
                quantum_obstructions_evaluated=[])


if __name__ == "__main__":
    assert __debug__
    print(json.dumps(calculate(), indent=2))
