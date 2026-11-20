"""Exact finite-colour test of a correction on one length-eleven content.

Use only parents whose complete split support avoids all previous successes.
Coordinates are coefficients in the original three-colour trace polynomial.
The RHS is the archived -Q1 O0; the linear system solves Q0 delta = RHS-Q0 O1.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from math import lcm
from pathlib import Path
import numpy as np
from flint import nmod_mat
import inhomogeneous_frontier as front
import complete_joint_minor as dixon
j,s=front.front.j,front.s
TARGETS=('0000000000010105050716','0000000001010106060d0f')

def paths(base):
    return (base/'negative_Q1.jsonl.gz',base/'joint_70/support_cut/preserve_296.vector.json.gz',base/'joint_70/support_cut/preserve_296.census.json')

def polynomials(rows,exe,work):
    """Expand all collected words. Coefficients are Python integers, with one denominator per row."""
    den=[lcm(*(v.denominator for v in row.values())) for row in rows];out=[Counter() for row in rows];words=sorted(set().union(*rows));
    for start in range(0,len(words),48):
        chunk=words[start:start+48]
        for w,v in zip(chunk,j.E.poly_iter(chunk,exe,work),strict=True):
            coefficients=list(v)
            for row,z,d in zip(rows,out,den):
                if w in row:
                    a=int(row[w]*d)
                    for k,b in coefficients:z[k]+=a*b
    result=[]
    for row,z,d in zip(rows,out,den):
        nz={k:a for k,a in z.items() if a};k=min(nz,default=None)
        result.append(dict(formal_words=len(row),collected_monomials=len(z),nonzero_monomials=len(nz),first_key=k,first_coefficient=None if k is None else str(Fraction(nz[k],d))))
    return result

def read_equations(rhs,vector,census,targets):
    digest,nr,rr=front.residuals(rhs,vector,targets);b={t:{} for t in targets};h=hashlib.sha256();n=0
    for w,v in s.stream_rows(rhs,h):
        n+=1;t=front.content(w)
        if t in b:b[t][w]=Fraction(v)
    assert n==517468 and h.hexdigest()==s.INPUT and dict(nr)=={t:len(b[t]) for t in targets}
    cc=s.read(census);ok=set(cc['exact_zero']);assert cc['vector_sha']==digest and len(ok)==313 and not cc['unresolved'];F,B=front.front.incidence(); groups={}
    for t in targets:
        assert t in cc['mathematically_nonzero'];groups[t]=sorted(g.hex() for g in B[bytes.fromhex(t)] if ok.isdisjoint(k.hex() for k in F[g]));assert groups[t]
    return b,rr,groups,digest

def certify(A,y,progress):
    """Select an integer minor mod p. Lift on the square and verify on the full rectangle.

    The archived correction has large rational denominators: neither the scaled
    RHS nor its first p-adic carries can be cast to int64. Only digits and matrix
    products with an explicit overflow bound enter machine arithmetic.
    """
    P=j.P;oi=np.argsort(np.diff(A.indptr),kind='stable');oj=np.argsort(np.diff(A.T.tocsr().indptr),kind='stable');C=[int(oj[i]) for i in j.piv(A[oi][:,oj])];R=[int(oi[i]) for i in j.piv(A[oi][:,C].T.tocsr())];assert len(C)==len(R)>0;B=A[R].T.tocsr();S=j.prime_matrix(B[C]);det=int(S.det());assert det
    sol=[int(k) for k in S.solve(nmod_mat(len(C),1,[y[k]%P for k in C],P)).entries()];del S;bad=[i for i,v in enumerate(j.multiply(B,sol)) if (y[i]-v)%P];rec=dict(prime=P,image_rank_mod_prime=len(R),augmented_rank_mod_prime=len(R)+bool(bad),minor_determinant_mod_prime=det);progress(rec)
    if bad:
        ncols=np.diff(A.T.tocsr().indptr);k=min(bad,key=lambda i:(ncols[i],i));T=A[:,C].tocsr();rhs=[int(t) for t in A[:,k].toarray().ravel()];info,pair=dixon.solve(T,rhs,R,progress)
        if pair is None:return dict(**rec,kind='uncertified',reason=info),{}
        nums,den=pair;assert j.multiply(T,nums)==[v*den for v in rhs];f={i:-Fraction(v,den) for i,v in zip(C,nums) if v};f[k]=Fraction(1);p=sum(v*y[i] for i,v in f.items());assert p;return dict(**rec,kind='exact_separator',integer_rhs_pairing=str(p),digits=info['digits']),f
    M=B[C].tocsr();b=[y[i] for i in C];H=dixon.bound(M,b);r=len(C);mass=max(sum(abs(int(z)) for z in B.data[i:k]) for i,k in zip(B.indptr,B.indptr[1:]));assert mass*(P-1)<2**63 and r*(P-1)**2<2**63
    inv=j.prime_matrix(M).inv();F=np.array([[int(inv[i,k]) for k in range(r)] for i in range(r)],dtype=np.int64);del inv;digits=0;xs=[0]*r;mod=1;cur=list(y);next_try=1
    while True:
        v=F@np.array([cur[i]%P for i in C],dtype=np.int64)%P;u=(B@v).tolist();diff=[a-bb for a,bb in zip(cur,u)];digits+=1
        assert all(diff[i]%P==0 for i in C)
        if any(a%P for a in diff):return dict(**rec,kind='uncertified',reason='chosen minor fails a rectangular row',digits=digits),{}
        cur=[a//P for a in diff];xs=[a+mod*bb for a,bb in zip(xs,v.tolist())];mod*=P
        if digits==next_try or mod>2*H*H:
            ans=dixon.recover(xs,mod,M,b,H)
            if ans:
                nums,den,_=ans;assert j.multiply(B,nums)==[a*den for a in y];return dict(**rec,kind='exact_solution',digits=digits),{i:Fraction(n,den) for i,n in zip(R,nums) if n}
            assert mod<=2*H*H;next_try=min(2*digits,digits+128)

def setup(target,rhs,vector,census,out,exe=None):
    out.mkdir(exist_ok=True,parents=True);b,r,gs,digest=read_equations(rhs,vector,census,[target]);ts=[target];dims={g:j.E.parity_frame_dims(bytes.fromhex(g))[0] for g in gs[target]+ts};D=lcm(*(a.denominator for a in r[target].values()));s.write(out/(target+'.rhs.json.gz'),[[s.hexes(w),str(int(-a*D))] for w,a in sorted(r[target].items())]);doc=dict(targets=ts,groups={g:ts for g in gs[target]},internal=gs[target],dimensions=dims,rhs=s.INPUT,base=digest,denominator=str(D));s.write(out/'graph.json',doc)
    if exe is not None:s.write(out/'physical_before.json',dict(zip(('negative_Q1','Q0_O1_plus_Q1_O0'),polynomials([b[target],r[target]],exe,out))))
    return doc

def solve(target,rhs,vector,census,out,exe):
    from su_trace_adapter_v1 import parity_project
    z,A,y,pool,keys=j.assemble(out,'all',out/'frames');assert len(keys)==z['dimensions'][target] and len(pool)==sum(z['dimensions'][g] for g in z['internal']);D=int(z['denominator']);history=[]
    def progress(x):history.append(x);s.write(out/'progress.json',history);print(target,x,flush=True)
    rec,vals=certify(A,y,progress);rec.update(target=target,target_dimension=len(keys),parent_dimension=len(pool),allowed=z['internal'],preserved_by_support=313,base=z['base'],rhs=z['rhs'],physical_before=s.read(out/'physical_before.json'));s.write(out/'result.json',rec)
    if rec['kind']=='exact_separator':
        rec.update(functional=[[keys[k][1],str(v)] for k,v in sorted(vals.items())],pairing_with_negative_residual=str(Fraction(rec['integer_rhs_pairing'])/D));s.write(out/'result.json',rec);return rec
    if rec['kind']!='exact_solution':return rec
    coeff={};frame_terms=[]
    for i,a in sorted(vals.items()):
        g,w=pool[i];v=a/D;frame_terms.append([g,w,str(v)]);terms=parity_project({s.expanded(s.unpack_hex(w)):1},1);j.l.p.axpy(coeff,{s.packed(u):q for u,q in terms.items()},v)
    base=j.l.load_vector(vector);j.l.p.axpy(base,coeff,1);j.l.save_vector(out/'new.vector.json.gz',base);j.l.save_vector(out/'delta.json.gz',coeff);s.write(out/'frame_terms.json.gz',frame_terms)
    F,_=front.front.incidence();ok=set(s.read(census)['exact_zero']);changed={front.content(w) for w in coeff};possible={k.hex() for g in changed for k in F[bytes.fromhex(g)]};assert changed<=set(rec['allowed']) and not ok&possible and target in possible
    _,_,rr=front.residuals(rhs,out/'new.vector.json.gz',[target]);rec.update(changed_parents=sorted(changed),affected_contents=sorted(possible),delta_words=len(coeff),vector_words=len(base),frame_terms=len(frame_terms),physical_after=polynomials([rr[target]],exe,out)[0]);assert rec['physical_after']['nonzero_monomials']==0;s.write(out/'result.json',rec);return rec

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('step',choices=('setup','solve'));p.add_argument('--target',required=True,choices=TARGETS);p.add_argument('--rhs',type=Path,default=paths(s.HERE)[0]);p.add_argument('--vector',type=Path,default=paths(s.HERE)[1]);p.add_argument('--census',type=Path,default=paths(s.HERE)[2]);p.add_argument('--out',required=True,type=Path);p.add_argument('--exe',required=True,type=Path);a=p.parse_args();assert __debug__;fn={'setup':setup,'solve':solve}[a.step];result=fn(a.target,a.rhs,a.vector,a.census,a.out,a.exe);print(a.target,result.get('kind','setup'),flush=True)
