"""Two dual target coordinates for the N=4 G2 instanton candidate.

These maps evaluate a *supplied* perturbative target vector. They compute
neither a physical instanton insertion nor the quantum mixing of its basis.
"""
import argparse,gzip,hashlib,json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from sympy import Matrix,zeros
from verify_g2_current_product import FIRST,SECOND
HERE=Path(__file__).resolve().parent
def exterior(x,y):
    r=Counter()
    for p,c in x.items():
        for q,v in y.items():
            if p&q:continue
            sg=(-1)**sum((q&((1<<k)-1)).bit_count() for k in range(p.bit_length()) if (p>>k)&1)
            r[p|q]+=sg*c*v
    return {p:c for p,c in r.items() if c}
def run(out):
    path=HERE.parent/'data/N4_G2_np_partner_replay_20261118.json.gz';raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest();assert digest=='a3d95d83b6e8ea512a972730f1747dcfa172a47068ff0dbc4b11e086bb7d7e5a';z=json.loads(gzip.decompress(raw));record={}
    for key,n in [('S2_evaluations',2),('Gram_evaluations',4)]:
        duals=[];common=None
        for saved in z[key]:
            a=Matrix([[Fraction(i) for i in row] for row in saved]);frame=a[:n,:];relations=a[n:,:];piv=frame.rref()[1];assert len(piv)==n and relations.rows==n
            co=relations[:,piv]*frame[:,piv].inv();assert co*frame==relations
            assert common is None or common==co;common=co
            unit=zeros(n,1);unit[n-1,0]=1;f=co.inv()*unit;ev=relations[:,piv].inv()*unit
            assert relations[:,piv]*ev==unit and frame[:,piv]*ev==f
            duals.append(dict(coordinate_indices=list(piv),coordinate_count=frame.cols,determinant=str(relations[:,piv].det()),evaluation_weights=[str(k) for k in ev],basis_coordinate_weights=[str(i) for i in f]))
        assert len(duals)==2 and duals[0]['basis_coordinate_weights']==duals[1]['basis_coordinate_weights'];record[key]=dict(boundary_and_normalized_target_rows=[[str(i) for i in row] for row in common.tolist()],two_independent_evaluations=duals)
    # Odd variables occur in order first-axis gaugino, second-axis gaugino.
    m=Matrix([[sum(x*y for x,y in zip(a,b)) for b in (FIRST,SECOND)] for a in (FIRST,SECOND)]);assert m==Matrix([[4,0],[0,12]])
    S=Counter()
    for f,g in zip(FIRST,SECOND):S.update(exterior({1:f,2:g},{4:f,8:g}))
    S={k:c for k,c in S.items() if c};S2=exterior(S,S);assert S=={5:4,10:12} and S2=={15:-96}
    # Every 3 by 3 Gram is A M A^T for an arbitrary 3 by 2 A; use
    # independent symbolic variables to check the determinant coefficientwise.
    from sympy import symbols,expand
    A=Matrix(3,2,symbols('x0:6'));assert expand((A*m*A.T).det())==0
    result=dict(input_sha256=digest,targets=record,cartan_trace_metric=m.tolist(),exterior_masks_of_trace_bilinear={str(k):c for k,c in S.items()},top_coefficient_of_squared_trace_bilinear=-96,scalar_Gram_Cartan_determinant=0,target_dimension=2,target_Cartan_rank=1,target_Cartan_kernel_dimension=1,scope='Independent duals for the two archived classical lengths. The Cartan map uses trace_7(lambda1 lambda2) squared and the scalar Gram determinant. Neither target functional has been applied to an instanton insertion.');result['cartan_trace_metric']=[[int(v) for v in row] for row in result['cartan_trace_metric']];out.write_text(json.dumps(result,indent=2)+'\n');print('Cartan:',-96,0,'target coordinates:',[record[k]['two_independent_evaluations'][0]['basis_coordinate_weights'] for k in record])
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--out',type=Path,required=True);assert __debug__;run(a.parse_args().out)
