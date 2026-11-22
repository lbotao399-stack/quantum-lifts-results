"""Full source-symbol gate for any SUN supported by the preserved engine."""
import sys,os,pickle,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'rational_c4'))
import exact_product as E
sys.path.insert(0,str(ROOT/'su3_rational'))
import fast_wick_edges,fast_selected_direct
fast_wick_edges.install(E.T);fast_selected_direct.install(E.T)
E.S.UP['q']=True
p=('p',0,0);q=('q',0,0);b=('b',0,0)
b1=('b',1,0);b2=('b',0,1);c1=('c',1,0);c2=('c',0,1)
def commtrace(a,c,d,db):
    return E.add(E.cyc_poly([a,c,d,db]),E.sc(E.cyc_poly([a,c,db,d]),-1),
                 E.sc(E.cyc_poly([c,a,d,db]),-1),E.cyc_poly([c,a,db,d]))
left=E.T.qn(1,E.mul(E.cyc_poly([p,b]),E.cyc_poly([q,b])))
right=E.sc(E.sub(E.add(commtrace(p,c1,q,b2),commtrace(q,c1,p,b2)),
                    E.add(commtrace(p,c2,q,b1),commtrace(q,c2,p,b1))),-1)
assert left==right
result={'status':'PASS','N':E.S.N,'component_terms':len(left),'independent_even_probes':2,
        'source_normalization':'unchanged','full_component_residual':0}
Path(__file__).with_name(f'symbol_SU{E.S.N}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
