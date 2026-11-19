"""Rational solve with a Cramer bound, adaptive denominators and full rectangular checks.

No conclusion about all source columns is drawn from a failure on a selected minor.
All archived integers refer to the original unscaled row and content order.
"""
import argparse, hashlib, sys
from fractions import Fraction
from math import isqrt, lcm, prod
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
import lift_joint_minor as old
j=old.jc
P=j.P


def bound(M,b):
    def squared_norms(S):
        return [sum(int(x)**2 for x in S.data[a:e]) for a,e in zip(S.indptr,S.indptr[1:])]
    rows=squared_norms(M); cols=squared_norms(M.T.tocsr()); assert min(cols)>0
    # For a numerator minor, replace one column by b. Both expressions bound
    # its determinant squared and the original determinant squared (Hadamard).
    sq=min(prod(a+int(v)**2 for a,v in zip(rows,b)),
           prod(cols)//min(cols)*max(min(cols), sum(int(v)**2 for v in b)))
    return isqrt(sq)


def recover(xs,m,M,b,guaranteed_bound):
    lim=isqrt(m//2); den=1; trials=0
    while True:
        v=[(x*den)%m for x in xs]; v=[x-m if 2*x>m else x for x in v]
        i=max(range(len(v)),key=lambda k:abs(v[k])); trials+=1
        if abs(v[i])<=lim and old.numerator_check(M,v,b,den):return v,den,trials
        r=j.reconstruct(xs[i],m)
        if r is None or den%r.denominator==0:break
        den=lcm(den,r.denominator)
        if den>lim:break
    if m<=2*guaranteed_bound**2:return None
    # Here coordinatewise reduced reconstruction is unique by Cramer's rule.
    rats=[j.reconstruct(x,m) for x in xs]; assert all(r is not None for r in rats)
    den=lcm(*(x.denominator for x in rats));v=[int(x*den) for x in rats]
    assert max([den]+list(map(abs,v)))<=guaranteed_bound and old.numerator_check(M,v,b,den)
    return v,den,trials


def solve(B,y,c,on_progress):
    r=len(c);M=B[c,:].tocsr();b=[y[i] for i in c];H=bound(M,b);goal=2*H*H;limit=1;g=P
    while g<=goal:limit+=1;g*=P
    row_mass=max(sum(abs(int(z)) for z in B.data[i:k]) for i,k in zip(B.indptr,B.indptr[1:]));cap=max(row_mass,max(map(abs,y)));assert cap+(P-1)*row_mass<2**63 and r*(P-1)**2<2**63
    on_progress(dict(stage='Cramer bound',square_order=r,coefficient_bound_bits=H.bit_length(),guaranteed_digits=limit,absolute_row_sum_max=row_mass))
    ref=j.prime_matrix(M);det=int(ref.det());assert det;R=ref.inv();del ref;inv=np.empty((r,r),dtype=np.int64)
    for i in range(r):inv[i]=[int(R[i,k]) for k in range(r)]
    del R;cur=np.array(y,dtype=np.int64);idx=np.array(c);xs=[0]*r;m=1;d=0;next_try=1
    def answer(ans):
        nums,den,trials=ans;assert den%P and old.numerator_check(M,nums,b,den)
        if not old.numerator_check(B,nums,y,den):
            v=j.multiply(B,nums);i=next(i for i,t in enumerate(v) if t!=y[i]*den);rest=y[i]*den-v[i];k=1
            while rest%P==0:rest//=P;k+=1
            return dict(status='selected minor excludes RHS on a rectangular row',digits=k,row=i,residual_digit=rest*pow(den,-1,P)%P,determinant_mod_prime=det,guaranteed_digits=limit),None
        return dict(status='exact rectangular solution',digits=d,determinant_mod_prime=det,guaranteed_digits=limit,coefficient_bound_bits=H.bit_length(),last_denominator_trials=trials,denominator_bits=den.bit_length(),max_numerator_bits=max(abs(v).bit_length() for v in nums)),(nums,den)
    while d<limit:
        digit=inv@(cur[idx]%P)%P;v=cur-B@digit;assert not np.any(v[idx]%P);bad=np.flatnonzero(v%P);d+=1
        if len(bad):
            # The selected nonsingular square system has a unique solution mod p^d.
            # The failing row has residual p^(d-1) times the nonzero digit below.
            i=int(bad[0]);return dict(status='selected minor excludes RHS on a rectangular row',digits=d,row=i,residual_digit=int(v[i]%P),determinant_mod_prime=det,guaranteed_digits=limit),None
        cur=v//P;xs=[x+m*t for x,t in zip(xs,digit.tolist())];m*=P
        if d==next_try or d==limit:
            ans=recover(xs,m,M,b,H)
            if ans:return answer(ans)
            on_progress(dict(stage='digit',digits=d));next_try=max(d+1,min(2*d,d+128))
    raise ArithmeticError('Cramer limit reached without a certified result')


def run(previous,top,out):
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(100000)
    z,A,y,pool,keys=j.assemble(previous,'internal',previous/'frames');assert A.shape==(9249,5475)
    nr=np.argsort(np.diff(A.indptr),kind='stable');nc=np.argsort(np.diff(A.T.tocsr().indptr),kind='stable');c=[int(nc[i]) for i in j.piv(A[nr,:][:,nc])];rr=[int(nr[k]) for k in j.piv(A[nr,:][:,c].T.tocsr())];assert len(c)==len(rr)==4308
    B=A[rr,:].T.tocsr();out.mkdir(exist_ok=True,parents=True);history=[]
    def progress(x):history.append(x);j.s.write(out/'precision.json',history);print(x,flush=True)
    res,ans=solve(B,y,c,progress);res.update(prime=P,base=z['base'],rhs=z['rhs'],parents=len(pool),target_coordinates=len(keys),pivot_order=len(c),pivot_sha256=hashlib.sha256(np.asarray(rr+c,dtype='<i8').tobytes()).hexdigest());j.s.write(out/'minor_indices.json.gz',dict(parent_rows=rr,coordinates=c))
    if ans:
        nums,den=ans;delta={};base=j.side.even(j.l.load_vector(top));assert j.side.canonical_digest([[j.s.hexes(w),str(v)] for w,v in sorted(base.items())])==z['base']
        for k,x in zip(rr,nums):
            if x:j.l.p.axpy(delta,j.E.projected(j.s.unpack_hex(pool[k][1])),Fraction(x,den*z['denominator']))
        assert j.side.even(delta)==delta;j.l.p.axpy(base,delta,1);j.l.save_vector(out/'vector.json.gz',base);D=lcm(*(q.denominator for q in base.values()));shared=dict(denominator=str(D),words=[[j.s.hexes(w),str(int(q*D))] for w,q in sorted(base.items())]);j.s.write(out/'vector_shared.json.gz',shared);res.update(integer_residual_coordinates=0,vector_words=len(base),delta_words=len(delta),vector_sha=j.side.canonical_digest(j.s.read(out/'vector.json.gz')),vector_denominator_bits=D.bit_length());j.s.write(out/'rational_pivot_solution.json.gz',dict(denominator=str(den),numerators=[str(n) for n in nums]))
    else:res['failed_coordinate_key']=keys[res['row']]
    j.s.write(out/'exact.json',res);return res


def controls():
    # The previous five probes and its weighted sum all miss the factor 7.
    n=10
    B=np.eye(n,dtype=np.int64);B[1,1]=B[3,3]=7;y=[0]*n;y[1]=2;y[3]=-1
    bs=csr_matrix(np.vstack((B,B[1]+B[3])));c=list(range(n));record,sol=solve(bs,y+[1],c,lambda v:None);assert sol and sol[1]==7 and sol[0][1]==2 and sol[0][3]==-1
    mismatch,other=solve(bs,y+[1+P**2],c,lambda v:None);assert other is None and mismatch['digits']==3 and mismatch['row']==10 and mismatch['residual_digit']==1
    return {'denominator':sol[1],'controlled_nonzero_numerators':[sol[0][1],sol[0][3]],'rectangular_mismatch_at_digit':mismatch['digits'],'mismatch_residue':mismatch['residual_digit']}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--previous',type=Path);p.add_argument('--vector',type=Path);p.add_argument('--out',type=Path);p.add_argument('--control',action='store_true');a=p.parse_args();assert __debug__
    if a.control:print(controls())
    else:print(run(a.previous,a.vector,a.out),flush=True)
