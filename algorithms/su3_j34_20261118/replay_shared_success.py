"""Reconstruct the shared rational correction from its five terms and replay both integral witnesses and the seven polynomials."""
import argparse, hashlib
from fractions import Fraction
from pathlib import Path
import shared_success_coupling as c
q=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--exe',type=Path,required=True);ap.add_argument('--work',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--jobs',type=int,default=2);a=ap.parse_args();assert __debug__
s=c.s; b=q/'joint_70';F=b/'shared_witness';cert=s.read(F/'preserve_303.certificate.json.gz');rec=s.read(F/'a8_5665405.json')['receipt'];audit=s.read(F/'preserve_303.audit.json');sw=s.read(b/'centralizer_sweep.json.gz');source,digest=c.scan.old.vector(b/'vector_shared.json.gz');x=dict(source);changes={};assert digest==cert['source']; assert s.read(F/'cost_at_most_three.certificate.json.gz')['frame_solution']==cert['frame_solution'];
for g,w,v in cert['frame_solution']:
 base=s.unpack_hex(w); assert bytes(sorted(b''.join(base))).hex()==g
 from su_trace_adapter_v1 import parity_project
 proj=parity_project({s.expanded(base):1},1);c.j.l.p.axpy(changes,{s.packed(u):k for u,k in proj.items()},Fraction(v))
c.j.l.p.axpy(x,changes,1);vv=c.j.l.load_vector(F/'preserve_303.vector.json.gz');assert vv==x and len(x)==cert['vector_words']==2728 and rec['files']['cost_at_most_three.vector.json.gz']==rec['files']['preserve_303.vector.json.gz'] and rec['files']['cost_at_most_three.audit.json']==rec['files']['preserve_303.audit.json']
side,fw,back,cost,n,sc,ds=c.graph(b/'census_2724.json');ok=set(s.read(b/'census_2724.json')['exact_zero']);changed={bytes(sorted(b''.join(w))).hex() for w in changes};affected=set.union(*(fw[g] for g in changed));assert affected&ok==set(side)
hv,rows=c.scan.old.residuals(q/'negative_Q1.jsonl.gz',b/'vector_shared.json.gz',c.FAIL)
for w,v in changes.items():
 for u,k in c.j.side.replay.Q_multi(s.expanded(w)).items():
  z=s.packed(u);g=bytes(sorted(b''.join(z))).hex()
  if g in rows:c.j.l.p.axpy(rows[g],{z:k},v)
vals=[]
for t in c.FAIL:
 w=c.rc.test(sw['targets'][t],2)['witness'];bits=[i for i in range(7) if w['mask']>>i&1];keys=set(bytes.fromhex(t));M={(l,i):w['matrices'][f'{l:02x}:{bits[i] if i>=0 else -1}'] for l in keys for i in (range(len(bits)) if c.f.ODD[l] else [-1])};v=sum(a*c.scan.ft.exterior_value(z,M,len(bits)) for z,a in rows[t].items());assert v==0;vals.append(str(v))
new=c.j.side.replay.replay(q/'negative_Q1.jsonl.gz', F/'preserve_303.vector.json.gz','joint',a.exe,a.work,a.jobs,map(bytes.fromhex,c.FAIL),set(map(bytes.fromhex,side)));assert new==audit
result=dict(vector_words=len(x),frame_terms=len(cert['frame_solution']),changed_parent_contents=sorted(changed),affected_previous_successes=sorted(affected&ok),unaffected_previous_successes=len(ok-affected),witnesses_over_Q=vals,polynomial_check=new);assert result==s.read(F/'independent.json');s.write(a.out,result);print(len(x), vals,[i['nonzero_monomials'] for i in new['holdouts']])
