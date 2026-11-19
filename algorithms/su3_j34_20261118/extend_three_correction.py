"""Extend the three-content primitive using parents with exactly five copies of 00.
Their classical differential is allowed only in the other 67 equations. Use
independent connected components, and arbitrary-size rational right-hand sides.
"""
import argparse, hashlib
from fractions import Fraction
from math import lcm
from pathlib import Path
import numpy as np
from flint import fmpz_mat,nmod_mat
import scalar_filtration as l
from physical_image_solve import canonical_digest,pivots
s,f,p,t=l.s,l.f,l.p,l.tri

def load(work,z):
    # Do not put a common-denominator RHS into scipy's int64 sparse storage.
    A,keys=t.load_columns(work,dict(z,rhs=[]));A=A[:-1,:].tocsc();out=[0]*len(keys);rhs=dict(z['rhs']);ind={int(k):j for j,k in enumerate(keys)}
    with (work/'polys.txt').open() as fp:
        for w,line in enumerate(fp):
            if w not in rhs:continue
            n,*row=map(int,line.split());assert n*2==len(row);a=rhs[w]
            for k,v in zip(row[::2],row[1::2]):out[ind[k]]+=a*v
    return A,keys,out

def solve(work,z,first_width=896):
    A,keys,b=load(work,z);n,m=A.shape;X=A.tocsr();order=np.random.default_rng(714293).permutation(m);must=[];test=[];width=min(m,first_width)
    while True:
        use=list(dict.fromkeys(must+order[:width].tolist()));F=A[:,use].toarray();bi=[b[i] for i in use];E=np.vstack([F,np.asarray([v%t.PRIME for v in bi],dtype=np.int64)]);cs=pivots(F);r=len(cs);aug=len(pivots(E));rows=pivots(F[:,cs].T);assert len(rows)==r;test.append(dict(width=len(use),rank=r,augmented_rank=aug));print('sample',len(use),r,aug,flush=True)
        if r!=aug:
            s.write(work/'restricted_gap.json',dict(samples=test,source_words=n,coordinates=m,scope='Only the selected parents at the indicated prime.'));return None
        a=F[rows,:][:,cs].T;B=fmpz_mat([[bi[j]] for j in cs]);M=fmpz_mat(a.tolist());u=M.solve(B);co=[Fraction(str(u[j,0])) for j in range(r)];den=lcm(*(c.denominator for c in co));res=[-den*y for y in b];
        for i,c in zip(rows,co):
            k=int(c*den)
            if not k:continue
            lo,hi=X.indptr[i:i+2]
            for j,v in zip(X.indices[lo:hi],X.data[lo:hi]):res[j]+=k*int(v)
        rem=[i for i,y in enumerate(res) if y];print('all integer coefficients',m,'nonzero',len(rem),flush=True)
        if not rem:
            det=int(nmod_mat((a%t.PRIME).tolist(),t.PRIME).det());assert det;ans={s.unpack_hex(z['source'][i]):c/z['den'] for i,c in zip(rows,co) if c};l.save_vector(work/'vector.json.gz',ans);s.write(work/'exact.json',dict(recipe=canonical_digest(z),contents=z['contents'],source_words=n,trace_words=len(z['word']),verified_coordinates=m,nonzero=0,samples=test,determinant_mod_1000003=det,vector=canonical_digest(s.read(work/'vector.json.gz')),used=len(ans)));return ans
        must=list(dict.fromkeys(must+rem[:24]));assert width<m,'inconsistent exact system despite modular equality';width=min(m,width+384)

def prepare(work,rhsfile,initial):
    hi={k for k in l.contents(11,0,l.CHARGES) if k.count(0)>=6};keep=hi|set(t.CONTENTS);src=s.read(initial);ans={s.unpack_hex(w):Fraction(c) for w,c in src['solution']};assert len(hi)==67 and len(keep)==70;b=p.rhs_rows(rhsfile,hi)
    for w,c in ans.items():p.axpy(b,{u:k for u,k in f.q(w).items() if p.content(u) in hi},-c)
    parents=sorted({w for u in b for w in f.parents(u) if p.content(w).count(0)==5});im={}
    for w in parents:
        v={u:k for u,k in f.q(w).items() if p.content(u) in keep}
        if v and all(p.content(u) in hi for u in v):im[w]=v
    rows=[]
    for i,(cc,ys,rhs) in enumerate(l.connected(im,b)):
        d=work/('component-%02d'%i);d.mkdir(parents=True,exist_ok=True);cert,co=l.formal(ys,rhs);entry=dict(component=i,contents=sorted(k.hex() for k in cc),formal=cert,source=len(ys),rhs=len(rhs));
        if cert['residual_terms']:l.recipe(d,ys,rhs,cc,'Physical equations on these complete high contents, leaving the earlier three contents fixed word for word.')
        else:l.save_vector(d/'vector.json.gz',co)
        rows.append(entry);print(i,len(cc),len(ys),cert,flush=True)
    s.write(work/'preparation.json',dict(components=rows,initial=canonical_digest(src['solution']),target=sorted(c.hex() for c in keep),nonzero_rhs=len(b),parent_rule='five 00; nonzero projected images must avoid all of the three original contents.'))

def combine(work,initial):
    prep=s.read(work/'preparation.json');src=s.read(initial)['solution'];assert canonical_digest(src)==prep['initial'];ans={s.unpack_hex(w):Fraction(c) for w,c in src};missing=[]
    for row in prep['components']:
        d=work/('component-%02d'%row['component']);file=d/'vector.json.gz'
        if not file.exists():missing.extend(row['contents']);continue
        co=l.load_vector(file);p.axpy(ans,co,1)
    l.save_vector(work/'combined_vector.json.gz',ans);passed=set(prep['target'])-set(missing);s.write(work/'combined.json',dict(source_terms=len(ans),enforced=sorted(passed),missing_selected=sorted(missing)));print('combined:',len(ans),len(passed),'contents',len(missing),'excluded',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('action',choices=['prepare','solve','combine']);a.add_argument('--work',type=Path,required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--initial',type=Path,default=s.HERE/'scalar_layers/three_solution.json.gz');a.add_argument('--component',type=int);a.add_argument('--jobs',type=int,default=2);a.add_argument('--exe',type=Path);a.add_argument('--width',type=int,default=896);a.add_argument('--reuse-columns',action='store_true');args=a.parse_args();assert __debug__
    if args.action=='prepare':prepare(args.work,args.rhs,args.initial)
    elif args.action=='combine':combine(args.work,args.initial)
    else:
        d=args.work/('component-%02d'%args.component);z=s.read(d/'recipe.json.gz');
        if not args.reuse_columns:t.columns(d,z,args.exe,args.jobs)
        solve(d,z,args.width)
