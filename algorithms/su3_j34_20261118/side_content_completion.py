"""Repair three disjoint sets of equations without disturbing the other 70-content rows.

Start at the C-even part of the 1201-word partial correction. Every fixed-content
source space whose letter splits touch the selected set but no other one of the
70 equations is included. Combinations from other source contents may also
cancel; we make no exclusion for them. Use independent character-saturated
frames for source and target, then check rational solutions on every coefficient.
"""
import argparse,hashlib
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from math import lcm
import numpy as np
from scipy.sparse import csr_matrix
from flint import fmpz_mat,nmod_mat
import even_content_frames as E
import replay_scalar_layers as replay
from physical_image_solve import canonical_digest,pivots
s,f,l,P=E.s,E.f,E.l,E.P

def sets(prep,which):
    universe=replay.domain('joint');cc={bytes.fromhex(c) for c in prep['components'][which]['contents']};near=defaultdict(set)
    for t in universe:
        for g in l.tri.out.incident(t):near[g].add(t)
    groups=sorted(g for g,q in near.items() if q<=cc);assert all(near[g] and g.count(0)==5 for g in groups) and len(groups)=={0:6,2:6,5:3}[which];return sorted(cc),groups

def even(v):
    z={}
    for w,c in v.items():l.p.axpy(z,E.projected(w),Fraction(c,2))
    return z

def run(rhs,vec,prep,which,exe,work,frames):
    cc,groups=sets(s.read(prep),which);work.mkdir(parents=True,exist_ok=True);zz={g:s.read(frames/(g.hex()+'.json.gz')) for g in cc+groups};vs={g:[s.unpack_hex(w) for w in zz[g]['words']] for g in zz};dims={g:E.parity_frame_dims(g)[0] for g in zz};assert all(zz[g]['content']==g.hex() and dims[g]==zz[g]['C_even_dimension']==len(vs[g])==len(set(zz[g]['polynomial_keys'])) and zz[g]['determinant_mod_prime'] for g in zz)
    co=even(l.load_vector(vec));b=l.p.rhs_rows(rhs,set(cc))
    for w,c in co.items():l.p.axpy(b,{t:n for t,n in f.q(w).items() if l.p.content(t) in cc},-c)
    assert b and even(b)==b;den=lcm(*(v.denominator for v in b.values()));ims=[];sources=[];pool=set(b)
    for g in groups:
        for w in vs[g]:
            v={}
            for u,n in E.projected(w).items():l.p.axpy(v,{t:m for t,m in f.q(u).items() if l.p.content(t) in cc},n)
            assert even(v)==v;ims.append(v);sources.append(w);pool.update(v)
    targets=sorted(pool);idx={w:i for i,w in enumerate(targets)};cols=[k for t in cc for k in zz[t]['polynomial_keys']];assert len(cols)==len(set(cols));colidx={k:i for i,k in enumerate(cols)}
    def sparse(rows,n):
        ind=[0];idx=[];data=[]
        for r in rows:
            for k,v in r:idx.append(k);data.append(v)
            ind.append(len(data))
        return csr_matrix((np.array(data,dtype=np.int64),idx,ind),shape=(len(ind)-1,n),dtype=np.int64)
    # Limit raw text on disk by streaming batches; RHS coefficients stay unbounded.
    def rows():
        for i in range(0,len(targets),400):
            for rr in E.poly_iter(targets[i:i+400],exe,work):yield [(colidx[k],v) for k,v in rr if k in colidx]
    V=sparse(rows(),len(cols));I=sparse(([(idx[w],v) for w,v in r.items()] for r in ims),len(targets));assert int(abs(V.data).max())*max(sum(abs(int(v)) for v in row.values()) for row in ims)<2**63;A=(I@V).toarray();y=[0]*len(cols)
    for w,c in b.items():
        k=int(c*den);i=idx[w];lo,hi=V.indptr[i:i+2]
        for j,v in zip(V.indices[lo:hi],V.data[lo:hi]):y[j]+=k*int(v)
    a=pivots(A);r=len(a);aug=len(pivots(np.vstack((A,np.array([v%P for v in y],dtype=np.int64)))));record=dict(component=which,prime=P,rhs_sha=s.INPUT,base_even_vector=canonical_digest([[s.hexes(w),str(c)] for w,c in sorted(co.items())]),target_dimensions={t.hex():dims[t] for t in cc},parent_dimensions={g.hex():dims[g] for g in groups},trace_words=len(targets),coordinates=len(cols),source_columns=len(sources),rank_mod_prime=r,augmented_rank_mod_prime=aug,scope='All C-even spaces of the listed contents. Other contents and combinations with nonzero individual side images are not excluded.');s.write(work/'rank.json',record);print('C-even rank',which,A.shape,r,aug,flush=True)
    if aug!=r:return record
    rows=pivots(A[:,a].T);assert len(rows)==r;minr=A[np.ix_(rows,a)].T;M=fmpz_mat(minr.tolist());res=M.solve(fmpz_mat([[y[j]] for j in a]));sol=[Fraction(str(res[i,0])) for i in range(r)];scale=lcm(*(x.denominator for x in sol));result=[-scale*v for v in y]
    for i,x in zip(rows,sol):
        m=int(scale*x)
        if m:
            for j in np.flatnonzero(A[i,:]):result[j]+=m*int(A[i,j])
    diff=sum(bool(x) for x in result);record['nonzero_exact_coordinates']=diff;assert diff==0,'rational verification failed despite modular equality';out={}
    for i,x in zip(rows,sol):l.p.axpy(out,E.projected(sources[i]),x/den)
    assert even(out)==out;l.save_vector(work/'delta.json.gz',out);total=dict(co);l.p.axpy(total,out,1);fil=work/'vector.json.gz';l.save_vector(fil,total);K=replay.domain('joint');audit=replay.replay(rhs,fil,'joint',exe,work/'audit',3,K-set(cc),cc);s.write(work/'replay.json',audit);allrows=audit['contents']+audit['holdouts'];record.update(minor_det_mod_prime=int(nmod_mat((minr%P).tolist(),P).det()),delta_words=len(out),combined_words=len(total),direct_replay_passed=sum(x['nonzero_monomials']==0 for x in allrows),replay=canonical_digest(audit),delta=canonical_digest(s.read(work/'delta.json.gz')),independent_rows=r);s.write(work/'rank.json',record);return record

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--vector',type=Path,default=s.HERE/'scalar_layers/extension_1201_vector.json.gz');a.add_argument('--prep',type=Path,default=s.HERE/'scalar_layers/extension_preparation.json');a.add_argument('--component',type=int,choices=[0,2,5],required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--frames',type=Path,required=True);a.add_argument('--exe',type=Path,required=True);b=a.parse_args();assert __debug__;print(run(b.rhs,b.vector,b.prep,b.component,b.exe,b.work,b.frames),flush=True)
