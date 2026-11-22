"""Return the G2 outgoing product after lowering the same dotted direction twice."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys


def run(repo, endpoint):
    sys.path.insert(0,str(repo/'algorithms/dual_singlets_20261118'))
    import mixed_incoming as m
    from sp3_omega_product_replay import trace_factor
    from so8_level24_descendant import substitute

    recipe_file=endpoint/'g2-d2-333-rows.json.gz'
    raw=recipe_file.read_bytes()
    packet=json.loads(gzip.decompress(raw))
    assert packet['weight']==[3,3,3,0,0] and packet['length']==8
    steps=[]
    for index,recipe in enumerate(packet['recipes']):
        filled=Counter()
        for factors,scale in recipe:
            term={():Fraction(scale)}
            for g in factors:
                pol=trace_factor(tuple(g))
                assert not m.q(pol)
                product=Counter()
                for w,c in term.items():
                    for tr,n in pol.items():
                        sg,new=m.so_canon(w+tr)
                        if sg:
                            product[new]+=sg*c*n
                term={w:c for w,c in product.items() if c}
            filled.update(term)
        x={w:c for w,c in filled.items() if c}
        assert x or index==3
        assert not m.q(x)
        expr=x
        supports=[len(expr)]
        for flavor,increase in ((2,False),(1,False),(1,True),(2,True)):
            expr=substitute(expr,1,flavor,increase)
            assert not m.q(expr)
            supports.append(len(expr))
        assert expr=={w:6*c for w,c in x.items()}
        steps.append(supports)
    return dict(recipes_sha256=sha256(raw).hexdigest(),
                selected_nonzero_recipe=0,return_coefficient=6,supports=steps)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--endpoint',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error('enable assertions')
    print(json.dumps(run(a.repo.resolve(),a.endpoint.resolve()),indent=2))
