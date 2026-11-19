"""Pointwise annihilators of every target-6 incoming block, also tested separately on target 10.

All projected one-letter replacements vanish as integer Grassmann matrices. A zero RHS on these functionals does not imply exactness.
"""
import argparse,gzip,hashlib,json,sys
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import numpy as np
import projected_q0 as maps
import edge_annihilator as edge
from coupled_first_replay import PRIME,RHS_SHA

def populate(alpha,n,tag,seed):
    fr,v=edge.int_frame(alpha,3,'diagonal_scalars',n,seed);rng=np.random.default_rng(seed+24);A=np.diag([1,1,-2]);L=np.zeros((3,3),dtype=np.int64);L[0,1]=1
    Es=[(1,2),(0,2),(0,1)];scalar=[(0,0,(i,)) for i in range(3)];order=[(0,0), (0,1), (1,0)]
    # Both the zeroth-order and first scalar jets are included. Q(f) contains the flavour signs (-,+,-).
    for k in range(n):
        if tag in ('regular_active','nilpotent_Lambda'):
            for i in range(3):
                for d in order:v[(*d,(i,))][k,0]=A*int(rng.integers(1,5))
            # Matrices in the simultaneous centralizer of the scalars, all jets, and lambda.
            if tag=='regular_active':
                for l in alpha:
                    if len(l[2]) in (2,3):v[l][k,:,2,:2]=0;v[l][k,:,:2,2]=0
            else:
                for l in alpha:
                    if len(l[2]) in (2,3):
                        for j in range(len(v[l][k])):
                            u,vv,w,y=rng.integers(-4,5,4);v[l][k,j]=[[u,vv,w],[0,u,0],[0,y,-2*u]]
                # All psi's in the centralizer of E12, with the weighted sum phi_i psi^i also in Z(A).
                # Choose independent psi_1,psi_2 and solve the off-block constraints by polynomial scaling.
                coeff=[int(v[(*order[0],(i,))][k,0,0,0]) for i in range(3)]
                # First jets are proportional to their base scalars, hence give the same moment equation.
                for d in order[1:]:
                    c=int(rng.integers(1,5))
                    for i,u in enumerate(coeff):v[(*d,(i,))][k,0]=c*u*A
                for xy in ((0,2),(2,1)):
                    x,y=xy;vals=[v[(0,0,s)][k,:,x,y].copy() for s in Es[:2]]
                    for i in range(2):v[(0,0,Es[i])][k,:,x,y]=coeff[2]*vals[i]
                    v[(0,0,Es[2])][k,:,x,y]=-coeff[0]*vals[0]+coeff[1]*vals[1]
                # derivative fermions enter target 10 via the derivative of Qf; impose the separate undifferentiated-scalar component there.
                for d in order[1:]:
                    for x,y in ((0,2),(2,1)):
                        vals=[v[(*d,s)][k,:,x,y].copy() for s in Es[:2]]
                        for i in range(2):v[(*d,Es[i])][k,:,x,y]=coeff[2]*vals[i]
                        v[(*d,Es[2])][k,:,x,y]=-coeff[0]*vals[0]+coeff[1]*vals[1]
            for l in alpha:
                if len(l[2])==0:
                    for j in range(len(v[l][k])):v[l][k,j]=(A if tag=='regular_active' else L)*int(rng.integers(-4,5))
        elif tag=='regular_moment':
            # phi, Dphi all diagonal; Dphi_i=c D^0phi_i. lambda=0. Solve (-,+,-) sum [phi_i,psi^i]=0 for generic off-diagonals.
            for d in order[1:]:
                c=int(rng.integers(1,5))
                for i in range(3):v[(*d,(i,))][k,0]=c*v[(0,0,(i,))][k,0]
            for l in alpha:
                if not l[2]:v[l][k]*=0
            for d in order:
                for x in range(3):
                    for y in range(3):
                        if x==y:continue
                        k0,k1,k2=[int(v[l][k,0,x,x]-v[l][k,0,y,y]) for l in scalar];u,w=[v[(*d,s)][k,:,x,y].copy() for s in Es[:2]]
                        v[(*d,Es[0])][k,:,x,y]=k2*u;v[(*d,Es[1])][k,:,x,y]=k2*w;v[(*d,Es[2])][k,:,x,y]=-k0*u+k1*w
        else:raise ValueError(tag)
    for l,x in v.items():fr.L[l]=(int(edge.is_ferm(l)),x%PRIME)
    return fr,v

def validate(v,s,t,alpha):
    from itertools import combinations
    k=v[(0,0,(0,1))].shape[1];rr=edge.matching(alpha,maps.pattern(t),maps.pattern(s));assert rr
    for L,terms in rr:
        total={};
        for c,a,b in terms:
            for i in range(v[a].shape[1]):
                x=1<<i if edge.is_ferm(a) else 0
                for j in range(v[b].shape[1]):
                    y=1<<j if edge.is_ferm(b) else 0
                    if x&y:continue
                    coef=-c if x and y and i>j else c
                    total[x|y]=total.get(x|y,0)+coef*(v[a][:,i]@v[b][:,j])
        assert all(not np.any(m) for m in total.values()),(s['index'],t['index'],L)
    return len(rr)

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));import sp_environment
    from fortuity.full.letters import Q_letter,is_ferm
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents
    from native_su_values_v1 import evaluate
    edge.Q_letter,edge.is_ferm,edge.Frame=Q_letter,is_ferm,Frame
    plan=json.loads((a.root/'classical_order1.json').read_bytes())['incoming'];mu=tuple(plan['mu']);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});targets=(6,10);by={maps.pattern(plan['targets'][t]):t for t in targets};records=[];frames={};values={};terms={t:0 for t in targets};buf={t:{} for t in targets}
    for tag in ('regular_active','nilpotent_Lambda','regular_moment'):
        seed=731321+109*len(frames);f,v=populate(alpha,a.points,tag,seed);maxentry=max(np.max(np.abs(m)) for m in v.values());assert 3*int(maxentry)**2*8*len(alpha)**2<2**63
        for t in targets:
            if tag=='regular_moment' and t==6:continue
            src=[s for s in plan['sources'] if t in s['targets']];counts={str(s['index']):validate(v,s,plan['targets'][t],alpha) for s in src};k=(tag,t);frames[k]=f;values[k]=np.zeros(f.ncols,dtype=np.int64);records.append(dict(frame=tag,target=t,seed=seed,points=a.points,integer_vanishing_replacement_counts=counts,matrix_sha256=hashlib.sha256(b''.join(v[l].astype('<i8').tobytes() for l in alpha)).hexdigest()))
    @lru_cache(None)
    def coef(c):x=Fraction(c);return x.numerator*pow(x.denominator,-1,PRIME)%PRIME
    def flush(t):
        for k,f in frames.items():
            if k[1]==t:values[k]=(values[k]+evaluate([buf[t]],f)[0])%PRIME
        buf[t]={}
    h=hashlib.sha256()
    with gzip.open(a.root/'negative_Q1.jsonl.gz','rb') as src:
        for line in src:
            h.update(line);c,ww=json.loads(line);w=maps.decode(ww);t=by.get(edge.pattern(l for tr in w for l in tr))
            if t is None:continue
            terms[t]+=1;buf[t][w]=coef(c)
            if len(buf[t])==512:flush(t)
    assert h.hexdigest()==RHS_SHA and terms=={6:31720,10:27780}
    for t in targets:
        if buf[t]:flush(t)
    for r in records:
        k=r['frame'],r['target'];v=values[k];r.update(width=len(v),rhs_nonzero=int(np.count_nonzero(v)),rhs_values=v.tolist());print(k,r['rhs_nonzero'],v[:3].tolist(),flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(dict(prime=PRIME,rhs_stream_sha256=RHS_SHA,records=records),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=maps.HERE);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--out',type=Path,required=True);p.add_argument('--points',type=int,default=6);run(p.parse_args())
