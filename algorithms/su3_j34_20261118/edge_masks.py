"""Systematically remove subsets of the five remaining target-6 blocks by exact, constrained matrix substitutions."""
import argparse,gzip,hashlib,json,sys,itertools
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
import projected_q0 as maps
import edge_annihilator as edge
from coupled_first_replay import PRIME,RHS_SHA

def make(alpha,mask,n,seed):
    f,v=edge.int_frame(alpha,3,'diagonal_scalars',n,seed);rng=np.random.default_rng(seed+1); sids=[3,4,7,8,9];ids=[s for s in sids if mask&(1<<sids.index(s))]; Es=[(1,2),(0,2),(0,1)];base=[(0,0,(i,)) for i in range(3)];A=np.diag([1,1,-2]);
    # Potentially nonzero rules by index: 3=[phi,psi], 4=[Dphi,psi], 7=[lambda,phi], 8={lambda,psi}, 9=[lambda,f].
    for k in range(n):
        if 7 not in ids or 8 not in ids:
            # Use a regular Cartan element if [lambda,phi]=0 is requested; in the other case choose a conjugate so phi need not commute with lambda.
            H=A.copy()
            if 7 in ids:
                U=np.eye(3,dtype=np.int64);U[0,2]=1;Ui=np.eye(3,dtype=np.int64);Ui[0,2]=-1;H=U@H@Ui
            for l in alpha:
                if not l[2]:
                    for j in range(len(v[l][k])):v[l][k,j]=H*int(rng.integers(1,5))
        else:H=None # lambda arbitrary
        if 8 not in ids:
            # For conjugate H, impose an integer centralizer; moment conditions below use independent diagonal phi, so instead choose diagonal psi when intersection is required.
            for se in Es:
                l=(0,0,se)
                for j in range(5):
                    M=v[l][k,j];M[:2,2]=0;M[2,:2]=0
                    if 7 in ids:M[:]=U@M@Ui
        if 9 not in ids:
            if H is None:
                H=A.copy()
                for l in alpha:
                    if not l[2]:
                        for j in range(len(v[l][k])):v[l][k,j]=H*int(rng.integers(1,5))
            l=(0,0,(0,1,2)); M=v[l][k,0];M[:2,2]=0;M[2,:2]=0
            if H[0,2]:M[:]=U@M@Ui
        # Impose the relevant moment equations through proportional scalar jets, choosing the base and first derivatives independently where useful.
        for d in ((1,0),(0,1)):
            if 3 not in ids and 4 not in ids:
                c=int(rng.integers(1,5));
                for i,l in enumerate(base):v[(*d,(i,))][k,0]=c*v[l][k,0]
            elif 4 not in ids:
                # All jets zero would destroy target 6; set common scalar colour A with one flavour vector.
                aa=(2,3,1)
                for i in range(3):v[(*d,(i,))][k,0]=aa[i]*A
        if 3 not in ids or 4 not in ids:
            K=(base if 3 not in ids else [(1,0,(i,)) for i in range(3)])
            # Safe polynomial elimination. If H is not diagonal, multiplication by different root coefficients may change [H,psi]. The exact rule test decides which subsets were really killed.
            for x in range(3):
                for y in range(3):
                    if x==y:continue
                    c0,c1,c2=[int(v[l][k,0,x,x]-v[l][k,0,y,y]) for l in K];
                    if not (c0 or c1 or c2):continue
                    w0,w1=[v[(0,0,s)][k,:,x,y].copy() for s in Es[:2]];v[(0,0,Es[0])][k,:,x,y]=c2*w0;v[(0,0,Es[1])][k,:,x,y]=c2*w1;v[(0,0,Es[2])][k,:,x,y]=-c0*w0+c1*w1
    for l,x in v.items():f.L[l]=(int(edge.is_ferm(l)),x%PRIME)
    return f,v

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));import sp_environment
    from fortuity.full.letters import Q_letter,is_ferm
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents
    from native_su_values_v1 import evaluate
    from edge_closed import validate
    edge.Q_letter,edge.is_ferm,edge.Frame=Q_letter,is_ferm,Frame
    root=a.root;plan=json.loads((root/'classical_order1.json').read_bytes())['incoming'];target=plan['targets'][6];mu=tuple(plan['mu']);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});pat=maps.pattern(target);ff=[];vv=[];records=[]
    for mask in range(32):
        seed=879122+131*mask;f,v=make(alpha,mask,a.points,seed); ff.append(f);vv.append(np.zeros(f.ncols,dtype=np.int64));possible=[]
        assert 3*int(max(np.max(np.abs(m)) for m in v.values()))**2*8*len(alpha)**2<2**63
        for s in plan['sources']:
            if 6 not in s['targets']:continue
            try:validate(v,s,target,alpha)
            except AssertionError as e:
                assert e.args and e.args[0][:2]==(s['index'],6);possible.append(s['index'])
        assert set(possible)<=set((3,4,7,8,9));records.append(dict(mask=mask,seed=seed,points=a.points,possible=possible,matrix_sha256=hashlib.sha256(b''.join(v[l].astype('<i8').tobytes() for l in alpha)).hexdigest()))
    # Concatenate point batches: trace prefixes shared across all substitutions.
    frame=ff[0];frame.P=a.points*len(ff);frame.ncols=frame.P*frame.cpp;frame.L={l:(frame.L[l][0],np.concatenate([f.L[l][1] for f in ff])) for l in alpha};val=np.zeros(frame.ncols,dtype=np.int64);buf={};terms=0;h=hashlib.sha256()
    @lru_cache(None)
    def coef(s):q=Fraction(s);return q.numerator*pow(q.denominator,-1,PRIME)%PRIME
    with gzip.open(root/'negative_Q1.jsonl.gz','rb') as src:
        for line in src:
            h.update(line);c,ww=json.loads(line);w=maps.decode(ww)
            if edge.pattern(l for tr in w for l in tr)!=pat:continue
            terms+=1;buf[w]=coef(c)
            if len(buf)==512:val=(val+evaluate([buf],frame)[0])%PRIME;buf={}
    if buf:val=(val+evaluate([buf],frame)[0])%PRIME
    assert terms==31720 and h.hexdigest()==RHS_SHA
    for r,q in zip(records,val.reshape(32,-1)):
        r['values']=q.tolist();r['nonzero']=int(np.count_nonzero(q))
    pos={tuple(x['possible']) for x in records if x['nonzero']};minima=sorted(k for k in pos if not any(set(j)<set(k) for j in pos));print('minimal necessary alternatives:',minima,flush=True);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(dict(prime=PRIME,rhs_stream_sha256=RHS_SHA,minimal_necessary_alternatives=minima,records=records),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=maps.HERE);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--points',type=int,default=2);p.add_argument('--out',type=Path,required=True);run(p.parse_args())
