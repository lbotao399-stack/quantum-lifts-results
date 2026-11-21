"""Verify the level-36 G2 index and the adjacent classical ranks through level 28."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from bmn_primary_index import cochains
from g2_primary_targets import calculate as target_cochains


PUBLIC = {
    "g2_hw_cohomology_j12_extra.jsonl": "116e0316f5c34625d42a878ae72cd861e91c8428",
    "g2_hw_cohomology_j18.jsonl": "91e1eda7c189b1bb6bc18678bfe68c4e1e9f3065",
    "g2_hw_cohomology_j20.jsonl": "1b88cd63bc4e28290c268b3b5dea115f76d58d99",
    "g2_hw_cohomology_j22_j24.jsonl": "dcbf1a0d5f68c3dfc7b5b0963608b92a70790d4d",
    "g2_hw_cohomology_j26.jsonl": "856861637b81c867178781690f96798d1fa27ece",
    "g2_hw_cohomology_j28.jsonl": "02638ffe6526bf975756c4d7db0ace948dbdc1aa",
}


def read(path):
    return json.loads(path.read_bytes())


def public_ranks(path):
    incoming = {}
    for name, digest in PUBLIC.items():
        raw = (path/name).read_bytes()
        header = f"blob {len(raw)}\0".encode()
        if hashlib.sha1(header+raw).hexdigest() != digest:
            raise ValueError(f"upstream Git blob mismatch for {name}")
        for line in raw.splitlines():
            row = json.loads(line)
            if row["group"] != "G2" or not row["check"]:
                continue
            weight = tuple(row["lam"])
            # r_up in the Fortuity data denotes the image entering this length.
            key = weight, sum(weight)-row["n"]
            incoming[key] = max(incoming.get(key,0),row["r_up"])
    return incoming


def extract_odd(rows, incoming):
    records = []
    for length, weights in sorted(rows.items()):
        if length%2:
            for weight, (dim,) in sorted(weights.items()):
                if 2*sum(weight) <= 28:
                    records.append([list(weight),length,dim,
                                    incoming.get((weight,length),0),incoming.get((weight,length+1),0)])
    return records


def receipt(job):
    folder = HERE/job
    record = read(folder/"receipt.json")
    assert record["job"] == job
    for name, digest in record["files"].items():
        assert Path(name).name == name and hashlib.sha256((folder/name).read_bytes()).hexdigest() == digest


def verify(upstream=None):
    for job in ("5690901", "5690811", "5691072"):
        receipt(job)
    levels = cochains(18,"G2")[1]
    saved = read(HERE/"odd_through_28.json")
    if upstream is not None:
        assert extract_odd(levels,public_ranks(upstream)) == saved["rows"]
    expected = {(w,y): v[0] for y,weights in levels.items() if y%2 for w,v in weights.items() if sum(w)<=14}
    seen = set()
    for weight,y,dim,incoming,outgoing in saved["rows"]:
        key = tuple(weight), y
        assert key not in seen and expected[key] == dim == incoming+outgoing
        assert min(incoming,outgoing) >= 0
        seen.add(key)
    assert seen == set(expected) and len(seen) == 293
    ordinary = [sum((-1)**y * chars.get((6,6,k),[0])[0] for y,chars in levels.items()) for k in range(7)]
    assert ordinary == [1,3,7,11,11,9,3]
    table = {tuple(row["weight"]):row for row in read(HERE/"5690901/exceptional.json")["rows"]}
    shape = (6,6,6)
    counts = {y:sum((-1)**(6-k)*chars.get((6,6,k),[0])[0] for k in range(7)) for y,chars in levels.items()}
    assert {str(y):[n] for y,n in counts.items() if n} == table[shape]["cochains_by_length"]
    assert sum((-1)**y * n for y,n in counts.items()) == table[shape]["primary_euler"][0] == -1
    for order in range(6):
        result = target_cochains(shape, order)
        assert json.loads(json.dumps(result)) == read(HERE/"5691072"/f"order{order}.json")
        if order == 0:
            assert result["cochains_by_length"] == {
                y:levels[y].get(shape,[0])[0] for y in range(1,19)}
    even_ancestors = [table[w]["primary_euler"][0] for w in ((5,5,5),(6,5,5),(6,6,5))]
    assert even_ancestors == [1,1,2] and sum(even_ancestors)-ordinary[-1] == 1
    records = [read(HERE/"5690811"/f"level{n}.json") for n in (32,34)]
    # Each original column and returned column is expressed in the same pair of rational functionals.
    def matrix(field):
        return list(zip(*(map(Fraction,record[field]) for record in records)))
    (a,b),(c,d) = matrix("original_pairings")
    (e,f),(g,h) = matrix("returned_pairings")
    determinant = a*d-b*c
    result = [[(d*e-b*g)/determinant,(d*f-b*h)/determinant],
              [(a*g-c*e)/determinant,(a*h-c*f)/determinant]]
    assert result == [[Fraction(35,3),Fraction(221,108)],[0,Fraction(25,4)]]
    return dict(ordinary_indices=ordinary,primary_even_cochains=sum(n for y,n in counts.items() if y%2==0),
                primary_odd_cochains=sum(n for y,n in counts.items() if y%2),odd_blocks_saturated=len(seen),
                G2_odd_primary_at_level_36_at_least=1,SU3_return_matrix=[[str(x) for x in r] for r in result])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream",type=Path,help="data/g2 directory at Fortuity commit 4256f207")
    args=parser.parse_args()
    if not __debug__:
        parser.error("Enable the certificate checks")
    print(json.dumps(verify(args.upstream),indent=2))
