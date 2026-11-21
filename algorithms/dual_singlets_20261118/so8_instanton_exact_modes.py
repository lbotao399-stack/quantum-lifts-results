"""One-instanton exact-mode constraints for the bare level-24 SO(8) two-point function."""
import argparse
from collections import Counter
import gzip
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path

from so8_instanton_modes import standard


def determinant(rows):
    size=len(rows)
    total=0
    for perm in permutations(range(size)):
        value=(-1)**sum(x>y for i,x in enumerate(perm) for y in perm[i+1:])
        for i,j in enumerate(perm):
            value*=rows[i][j]
        total+=value
    return total


def two_site_ranks():
    # At the first insertion the two zeta components and the single allowed
    # projected bare S-bar mode span three directions of a four-dimensional space.
    left=((1,0,0,0),(0,1,0,0),(0,0,1,0))
    # One invertible separation suffices to exhibit allowed four-forms.
    right=((1,0,1,0),(0,1,0,1))
    return {k:any(determinant(a+b) for a in combinations(left,k)
                  for b in combinations(right,4-k)) for k in range(5)}


def run(path):
    words=json.loads(gzip.decompress(path.read_bytes()))
    content=Counter()
    chiral=[(0,1),(0,2),(0,3),(0,2,3),(0,1,3),(0,1,2)]
    for coefficient,traces in words:
        assert coefficient
        labels=[i for tr in traces for i in tr]
        assert len(labels)==10 and sum(i<3 for i in labels)==8 and all(i in range(6) for i in labels)
        count=tuple(sum(j in chiral[i] for i in labels) for j in range(4))
        assert count==(10,4,4,4) and tuple(10-v for v in count)==(0,6,6,6)
        content[tuple(labels.count(i) for i in range(6))]+=1
    assert len(words)==13449
    allowed=two_site_ranks()
    assert {i for i,v in allowed.items() if v} == {2,3}
    rows=[]
    for left_extra in range(3):
        for right_extra in range(3-left_extra):
            quartic=2-left_extra-right_extra
            totals_left=tuple(i+left_extra for i in (10,4,4,4))
            totals_right=tuple(i+right_extra for i in (0,6,6,6))
            kept=[]
            for left_exact in product((2,3),repeat=4):
                left=tuple(n-k for n,k in zip(totals_left,left_exact))
                right=tuple(n-(4-k) for n,k in zip(totals_right,left_exact))
                if min(*left,*right)<0 or max(*left,*right)>8:
                    continue
                assert all(a+b+quartic==8 for a,b in zip(left,right))
                a,b=standard(left),standard(right)
                if a and b:
                    kept.append(dict(left_exact=left_exact,left_lifted=left,right_lifted=right,
                                     left_standard_multiplicity=a,right_standard_multiplicity=b))
            rows.append(dict(left_extra_quartets=left_extra,right_extra_quartets=right_extra,
                             action_quartic_insertions=quartic,assignments=kept))
    assert all(not r["assignments"] for r in rows if r["right_extra_quartets"]==0)
    return dict(operator_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),terms=len(words),
                contents=[[list(k),n] for k,n in sorted(content.items())],
                left_minimal=[10,4,4,4],right_minimal=[0,6,6,6],
                exact_modes_per_flavour=4,lifted_modes_per_flavour=8,
                left_exact_degrees_with_four_form=[k for k in range(5) if allowed[k]],
                standard_channels=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operator",type=Path,default=Path(__file__).with_name("so8_level24")/"even_bmn.json.gz")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    record=run(args.operator)
    args.out.write_text(json.dumps(record,indent=2)+"\n")
    for row in record["standard_channels"]:
        print(row["left_extra_quartets"],row["right_extra_quartets"],row["action_quartic_insertions"],len(row["assignments"]))
