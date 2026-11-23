"""Record the exact coefficient probes of an already computed minor."""
from fast_paths import *
import argparse,json,pickle
import abelian_obstruction as ab
from flint import fmpz_mat

ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,default=6)
a=ap.parse_args();rank=a.rank
folder=HERE/f'small_sp{rank}'
certificate=json.loads((folder/'certificate.json').read_text())
coordinates=json.loads((folder/'rational_coordinates.json').read_text())
selected=coordinates['combined_basis'];rows=[]
for i in selected:
    _,row=pickle.loads((HERE/f'frames_sp{rank}/L{2*rank+4}/row_{i}.pkl').read_bytes())
    _,q=pickle.loads((HERE/f'frames_sp{rank}/Q{2*rank+4}/row_{i}.pkl').read_bytes())
    restriction=ab.restrict(row)
    rows.append({**{('Q0',w):c for w,c in q.items()},**{('Cartan',w):c for w,c in restriction.items()}})
cert=ab.rank_certificate(rows,1000003);assert cert['rank']==len(selected)==61
det=int(fmpz_mat([[row.get(w,0) for w in cert['columns']] for row in rows]).det())
assert str(det)==certificate['combined_exact_minor_determinant']
certificate['combined_minor_rows']=selected
certificate['combined_minor_columns']=[(kind,pf.m.encode(w)) for kind,w in cert['columns']]
(folder/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print('captured',rank,len(selected),'coefficient probes; determinant',det,flush=True)
