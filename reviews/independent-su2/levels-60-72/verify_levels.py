"""Verify both supercharge functionals on every finite-colour singlet predecessor."""
from argparse import ArgumentParser
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

import gauge_frame as gauge
from first_frame import d
from first_character import calculate as character
from su2_bmn_character50 import characters

HERE=Path(__file__).resolve().parent


def check(n,all_columns=False):
    w=3+2*n
    y=7+2*n
    level=24+12*n
    leading_weight=(w+1,)*3
    src=d.frame.run(leading_weight,0,y-1,y+1)
    c=characters([leading_weight])[','.join(map(str,leading_weight))]
    assert all(c[str(i)]==src[i]['dimension'] for i in range(y-1,y+2))
    source=src[y]
    assert source['cohomology_upper']==2
    data=json.loads(gzip.decompress((HERE/f'duals{level}-native.json.gz').read_bytes()))
    duals=[{(tuple(a),tuple(b),mask):c for (a,b,mask),c in h['terms']} for h in data]
    count=0
    for profile in gauge.profiles(w,y-1):
        basis=gauge.invariant(*profile,0)
        for i,row in enumerate(basis):
            image=gauge.differential_twice(row,profile[-1])
            if i==0 or all_columns:assert not gauge.full_tree(image,w,y)
            for dual in duals:
                assert sum(c*image.get(k,0) for k,c in dual.items())==0
        count+=len(basis)
    assert count==character(n)[y-1]['dimension']
    matrix=[[],[]]
    for poly in (d.core(n),d.dressed(n-1)):
        assert not d.act('Q',0,poly)
        quantum=gauge.frame.quantum(poly)
        assert not gauge.frame.tree(quantum)
        image=gauge.scalar_field_slice(quantum)
        assert not gauge.full_tree(image,w,y)
        for j,dual in enumerate(duals):
            matrix[j].append(sum(c*image.get(k,0) for k,c in dual.items()))
    assert [[str(v) for v in row] for row in matrix]==[h['scores'] for h in data]
    det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    assert det
    return dict(level=level,flavour=list(leading_weight),length=y,
                classical=[source['dimension'],source['incoming'],source['out'],2],
                singlet_predecessors=count,dual_terms=[len(v) for v in duals],
                first_matrix=[[str(v) for v in row] for row in matrix],
                determinant=str(det),second_page_dimension=0)


if __name__=='__main__':
    a=ArgumentParser(description=__doc__)
    a.add_argument('n',nargs='*',type=int,default=[2,3,4],choices=(2,3,4))
    a.add_argument('--all-columns',action='store_true')
    opts=a.parse_args()
    if not __debug__:a.error('enable assertions')
    print(json.dumps([check(n,opts.all_columns) for n in opts.n],indent=2))
