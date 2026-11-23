"""Native Weyl characters, rank-one Ward rows and Levi bounds for the focused powers."""
from functools import cache
from itertools import combinations
import json
from pathlib import Path
import sys
from wall_test import run as ward
from source_sizes import rows
from b2_schur_complex import transpose
from fractions import Fraction

HERE=Path(__file__).resolve().parent
page=HERE.parent/'first_page'
sys.path.insert(0,str(page))
from fifth_cartan_mixing import schur,weyl_spectra
from commuting_source import simple_coordinates
from colour import Colour

def pieces(shape):
    out=set()
    for i,j in combinations(range(len(shape)),2):
        m=tuple(v-int(k in (i,j)) for k,v in enumerate(shape))
        if sorted(m,reverse=True)==list(m):
            out.add(tuple(x for x in m if x))
    return sorted(out)

def count(s,mark):
    numerator=sum(n*schur(s,t)*(t[0] if mark else 1) for t,n in weyl_spectra().items())
    assert numerator%1152==0
    return numerator//1152

def result():
    system=simple_coordinates(Colour())[4:];pos={r for r in system if min(r)>=0}
    assert len(pos)==24
    maxheight=[]
    for n in range(5):
        choices=list(combinations(range(4),n))
        maxheight.append(max((sum(r) for r in pos for chosen in choices
                              if all(j in chosen or not a for j,a in enumerate(r))),default=0))
    assert maxheight==[0,1,3,5,11]
    walls=ward()
    matrix=[[Fraction(v) for v in walls[i]['equations']] for i in (0,2)]
    det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    assert det==126
    records=[]
    for p in (6,7):
        for info in rows(p):
            if info['order']==1:continue
            s=transpose(info['shape']);co=pieces(s)
            mult=[count(q,True) for q in co]
            minima=[min(sum(q[4-l:]) for q in co) for l in range(5)]
            assert all(minima[l]>maxheight[l] for l in (2,3,4))
            known_zero=sum(mult)==0 or (p,info['order'])==(6,3)
            records.append(dict(**info,ordinary_shape=s,scalar_shapes=co,marked_cartan_multiplicities=mult,
                                terminal_cartan_dimension=count(s,False),min_positive_slots=minima,
                                classical_cycles_vanish_on_reduced_commuting_locus=known_zero,
                                regular_chart_denominator_exponents=[sum(s)-j for j in s] if known_zero else None))
    assert [r['marked_cartan_multiplicities'] for r in records]==[[0,0],[2],[0,0,0],[1,0],[0]]
    return dict(weyl_order=sum(weyl_spectra().values()),levi_maximum_heights=maxheight,root_Ward=walls,
                long_short_determinant=str(det),source_tests=records)
if __name__=='__main__':
    print(json.dumps(result(),indent=2))
