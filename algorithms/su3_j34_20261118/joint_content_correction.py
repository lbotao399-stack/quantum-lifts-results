"""Couple the middle equation to side blocks 0, 2 and 5. Run in three steps.

prepare: enumerate every unordered predecessor by reverse merge, and independently
by forward Q_letter on all 581 allowed length-ten contents. Save the fixed C-even
residual. block: export integer entries on one target's detecting monomial frame.
solve: compare (i) parents confined to the 12 equations among the original 70 and
(ii) ALL parents of the 12, allowing spill into the other 58. A modular gap alone
is never reported as a rational obstruction. Reconstructed vectors and separating
functionals must pass exact integer multiplication over every row or column.
"""
import argparse,sys
from collections import defaultdict
from fractions import Fraction
from math import isqrt,lcm
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,vstack,save_npz,load_npz,hstack
from flint import nmod_mat
import side_content_completion as side
E,l,s,P=side.E,side.l,side.s,side.P

def prepare(rhs,top,prep,out):
    universe=side.replay.domain('joint'); raw=s.read(prep)['components']; selected={l.tri.MID}
    for i in (0,2,5):selected.update(map(bytes.fromhex,raw[i]['contents']))
    backward=defaultdict(set)
    for t in universe:
        for g in l.tri.out.incident(t):backward[g].add(t)
    forward={}
    gs=l.contents(10,0,l.CHARGES);assert len(gs)==581
    for g in gs:
        targets=set()
        for j,k in enumerate(g):
            for _,pair in l.f.SPLITS[k]:
                t=bytes(sorted(g[:j]+g[j+1:]+pair))
                if t in universe:targets.add(t)
        if targets:forward[g]=targets
    assert forward==dict(backward)
    groups=sorted(g for g,ts in forward.items() if not ts.isdisjoint(selected));internal=sorted(g for g in groups if forward[g]<=selected);ts=sorted(selected)
    assert len(ts)==12 and len(internal)==19 and len(groups)==36
    x=side.even(l.load_vector(top));b=l.p.rhs_rows(rhs,selected)
    for w,c in x.items():l.p.axpy(b,{v:n for v,n in l.f.q(w).items() if l.p.content(v) in selected},-c)
    assert b and side.even(b)==b;den=lcm(*(v.denominator for v in b.values()));assert den%P
    for t in ts:s.write(out/(t.hex()+'.rhs.json.gz'),[[s.hexes(w),str(int(den*v))] for w,v in sorted(b.items()) if l.p.content(w)==t]);
    z=dict(prime=P,rhs=s.INPUT,base=side.canonical_digest([[s.hexes(w),str(c)] for w,c in sorted(x.items())]),denominator=den,targets=[t.hex() for t in ts],internal=[g.hex() for g in internal],groups={g.hex():[t.hex() for t in sorted(forward[g])] for g in groups},dimensions={g.hex():E.parity_frame_dims(g)[0] for g in ts+groups});assert sum(z['dimensions'][g.hex()] for g in groups)==11311 and sum(z['dimensions'][g.hex()] for g in internal)==9249 and sum(z['dimensions'][t.hex()] for t in ts)==5475;s.write(out/'graph.json',z);return z

def sparse(rows,n):
    ix=[];data=[];bounds=[0]
    for rr in rows:
        for j,v in rr:ix.append(j);data.append(v)
        bounds.append(len(ix))
    return csr_matrix((np.asarray(data,dtype=np.int64),ix,bounds),shape=(len(bounds)-1,n))

def matrix_for_target(target,frames,exe,out):
    z=s.read(out/'graph.json');assert target in z['targets'];t=bytes.fromhex(target);gs=sorted(g for g,ts in z['groups'].items() if target in ts);fr={k:s.read(frames/(k+'.json.gz')) for k in gs+[target]}
    for k,a in fr.items():assert a['content']==k and a['determinant_mod_prime'] and z['dimensions'][k]==len(a['words'])==len(set(a['polynomial_keys']))==a['C_even_dimension']
    keys=fr[target]['polynomial_keys'];jx={v:j for j,v in enumerate(keys)};b={s.unpack_hex(w):int(v) for w,v in s.read(out/(target+'.rhs.json.gz'))};images=[];pool=set(b)
    for g in gs:
        for i,w in enumerate(fr[g]['words']):
            w=s.unpack_hex(w);im={}
            for h,n in E.projected(w).items():
                full=l.f.q(h);l.p.axpy(im,{u:m for u,m in full.items() if l.p.content(u)==t},n)
                if i%191==0:assert full=={s.packed(v):a for v,a in side.replay.Q_multi(s.expanded(h)).items()}
            images.append(im);pool.update(im)
    words=sorted(pool);ix={w:i for i,w in enumerate(words)}
    def row_iter():
        for a in range(0,len(words),320):
            for rr in E.poly_iter(words[a:a+320],exe,out):yield ((jx[k],v) for k,v in rr if k in jx)
    V=sparse(row_iter(),len(keys));I=sparse(([(ix[w],v) for w,v in im.items()] for im in images),len(words));bound=max((sum(abs(v) for v in r.values()) for r in images),default=0)*max((abs(int(x)) for x in V.data),default=0);assert bound<2**63
    A=(I@V).tocsr();A.eliminate_zeros();A.sort_indices();y=[0]*len(keys)
    for w,c in b.items():
        a,h=V.indptr[ix[w]:ix[w]+2]
        for j,v in zip(V.indices[a:h],V.data[a:h]):y[j]+=c*int(v)
    save_npz(out/(target+'.npz'),A);s.write(out/(target+'.map.json'),dict(groups=gs,y=[str(n) for n in y],cols=len(keys),trace_words=len(words),entries=A.nnz,frames={k:side.canonical_digest(v) for k,v in fr.items()}));print('matrix',target,A.shape,A.nnz,flush=True)

def assemble(out,scope,frames):
    z=s.read(out/'graph.json');ts=z['targets'];gs=sorted(z['groups']) if scope=='all' else z['internal'];ds=z['dimensions'];blocks=[];rhs=[];pool=[];keys=[]
    for g in gs:pool.extend((g,w) for w in s.read(frames/(g+'.json.gz'))['words'])
    for t in ts:
        m=s.read(out/(t+'.map.json'));a=load_npz(out/(t+'.npz'));rhs.extend(map(int,m['y']));keys.extend((t,k) for k in s.read(frames/(t+'.json.gz'))['polynomial_keys']);offs={};i=0
        for g in m['groups']:offs[g]=i;i+=ds[g]
        assert a.shape==(i,ds[t]);blocks.append(vstack([a[offs[g]:offs[g]+ds[g],:] if g in offs else csr_matrix((ds[g],ds[t]),dtype=np.int64) for g in gs],format='csr'))
    A=hstack(blocks,format='csr');assert A.shape==(len(pool),len(rhs)) and len(keys)==len(set(keys))==len(rhs);return z,A,rhs,pool,keys

def prime_matrix(A):
    rows=A.toarray() if hasattr(A,'toarray') else A
    return nmod_mat((rows%P).tolist(),P)

def piv(A):
    R=prime_matrix(A);R,k=R.rref(inplace=True);cs=[];j=0
    for i in range(k):
        while not R[i,j]:j+=1
        cs.append(j);j+=1
    return cs

def reconstruct(a,mod):
    lim=isqrt(mod//2);n,r=mod,a%mod;v,w=0,1
    while r>lim:
        q=n//r;n,r=r,n-q*r;v,w=w,v-q*w
    if not w or abs(w)>lim:return None
    t=Fraction(r,w)
    return t if (a*t.denominator-t.numerator)%mod==0 else None

def multiply(A,x):
    # Never coerce arbitrary-sized integer weights to a fixed-width dtype.
    return [sum(int(A.data[k])*x[A.indices[k]] for k in range(i,j)) for i,j in zip(A.indptr,A.indptr[1:])]

def exact_solve(M,b,steps,accept):
    # Reuse the same p-adic inverse at each digit. Reconstruction is only a
    # proposal; accept checks the original, rectangular integer system.
    if not len(b):return [] if accept([]) else None
    S=prime_matrix(M).inv();current=list(b);mod=1;coeff=[0]*len(b)
    for j in range(steps):
        digits=[int(a) for a in (S*nmod_mat(len(b),1,[x%P for x in current],P)).entries()];
        md=multiply(M,digits);assert all((x-v)%P==0 for x,v in zip(current,md));current=[(x-v)//P for x,v in zip(current,md)];coeff=[x+mod*v for x,v in zip(coeff,digits)];mod*=P
        if j in (0,1,3,7,15,31,63,127,255) or j+1==steps:
            cand=[reconstruct(c,mod) for c in coeff]
            if None not in cand and accept(cand):return cand
    return None

def solve(out,frames,scope,steps):
    z,A,rhs,pool,keys=assemble(out,scope,frames);n,m=A.shape
    # Sparse-first ordering selects a better-conditioned integer minor. We keep
    # original word/monomial labels throughout, including for exact replay.
    nr=np.argsort(np.diff(A.indptr),kind='stable');nc=np.argsort(np.asarray(A.getnnz(axis=0)),kind='stable');At=A[nr,:][:,nc];print('rank',scope,n,m,A.nnz,flush=True)
    c=[int(nc[j]) for j in piv(At)];rr=[int(nr[i]) for i in piv(A[nr,:][:,c].T.tocsr())];assert len(rr)==len(c);r=len(rr);minor=A[rr,:][:,c].tocsr();bmod=[v%P for v in rhs];field=prime_matrix(minor.T);sol=[int(v) for v in field.solve(nmod_mat(r,1,[bmod[j] for j in c],P)).entries()];del field;pick=A[rr,:].T.tocsr();yp=multiply(pick,sol);nz=[j for j,(v,y) in enumerate(zip(yp,bmod)) if (v-y)%P];a=r+bool(nz);rec=dict(scope=scope,prime=P,parents=n,target_coordinates=m,projected_rank=r,projected_augmented_rank=a,rhs_sha=z['rhs'],base=z['base'],parent_contents=len(set(g for g,w in pool)),residual_coordinates_mod_prime=len(nz),exact_status='not attempted');dest=out/(scope+'.result.json');s.write(dest,rec);print('ranks',scope,r,a,flush=True)
    if not steps:return rec
    if nz:
        j=min(nz,key=lambda j:A[:,j].nnz);b=A[rr,j].toarray().ravel().tolist();C=A[:,c].tocsr();whole=A[:,j].toarray().ravel().tolist()
        def accept(vals):
            scale=lcm(*(x.denominator for x in vals));y=[int(x*scale) for x in vals];return multiply(C,y)==[v*scale for v in whole] and sum(rhs[t]*u for t,u in zip(c,y))!=rhs[j]*scale
        co=exact_solve(minor,b,steps,accept)
        if co is not None:
            pairs=[(keys[j],Fraction(1))]+[(keys[t],-x) for t,x in zip(c,co) if x];bvalue=rhs[j]-sum(rhs[t]*u for t,u in zip(c,co));assert bvalue
            s.write(out/(scope+'.separator.json.gz'),[[t,str(k),str(v)] for (t,k),v in pairs]);rec.update(exact_status='verified rational separator',nonzero_coefficients=len(pairs),rhs_pairing=str(Fraction(bvalue,z['denominator'])))
    else:
        def accept(vals):
            scale=lcm(*(x.denominator for x in vals));return multiply(pick,[int(x*scale) for x in vals])==[y*scale for y in rhs]
        co=exact_solve(minor.T.tocsr(),[rhs[j] for j in c],steps,accept)
        if co is not None:
            add={}
            for i,x in zip(rr,co):
                if x:l.p.axpy(add,E.projected(s.unpack_hex(pool[i][1])),x/z['denominator'])
            assert side.even(add)==add;l.save_vector(out/(scope+'.delta.json.gz'),add);rec.update(exact_status='verified rational frame solution',nonzero_trace_words=len(add))
    if co is None:rec.update(exact_status='no verified rational reconstruction at digit cap',digit_cap=steps)
    s.write(dest,rec);print(rec,flush=True);return rec

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','block','solve']);p.add_argument('--out',type=Path,required=True);p.add_argument('--frames',type=Path);p.add_argument('--exe',type=Path);p.add_argument('--target');p.add_argument('--scope',choices=('internal','all'));p.add_argument('--digits',type=int,default=32);p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');p.add_argument('--vector',type=Path,default=s.HERE/'scalar_layers/top_vector.json.gz');p.add_argument('--prep',type=Path,default=s.HERE/'scalar_layers/extension_preparation.json');a=p.parse_args();assert __debug__;a.out.mkdir(exist_ok=True,parents=True)
    if a.mode=='prepare':prepare(a.rhs,a.vector,a.prep,a.out)
    if a.mode=='block':matrix_for_target(a.target,a.frames,a.exe,a.out)
    if a.mode=='solve':solve(a.out,a.frames,a.scope,a.digits)
