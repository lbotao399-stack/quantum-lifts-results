"""Lift one modular minor to Q and verify the unsampled rectangular equation.

The inverse and every base-p digit product fit signed 64-bit integers by explicit
bounds. Accumulated coordinates, common-denominator recovery, and final products
use Python integers. A precision stop is never evidence against a rational solution.
"""
import argparse, hashlib, json
from fractions import Fraction
from math import isqrt,lcm
from pathlib import Path
import numpy as np
from flint import nmod_mat
import joint_content_correction as jc

P=jc.P

def numerator_check(A,x,y,d):
    # A is integral CSR with target coordinates as rows. Stop on first discrepancy.
    for i,yi in enumerate(y):
        lo,hi=A.indptr[i:i+2]
        if sum(int(a)*x[j] for a,j in zip(A.data[lo:hi],A.indices[lo:hi]))!=yi*d:return False
    return True

def rational_digit_solution(M,b,accept,max_digits,progress):
    r=M.shape[0];assert M.shape==(r,r) and r==len(b) and r*(P-1)**2<2**63
    bound=max((sum(abs(int(a)) for a in M.data[i:j]) for i,j in zip(M.indptr,M.indptr[1:])),default=0)
    # If |v| <= C then |(v-Mt)/p| <= C, since each digit t is in [0,p-1].
    cap=max(bound,max(map(abs,b),default=0));assert cap+(P-1)*bound<2**63
    R=jc.prime_matrix(M).inv(); inv=np.empty((r,r),dtype=np.int64)
    for i in range(r): inv[i]=[int(R[i,j]) for j in range(r)]
    del R
    cur=np.array(b,dtype=np.int64);powers=1;xs=[0]*r;probes=sorted(set([0,r//5,r//2,4*r//5,r-1]));expanded={4,8,16,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096,6144,8192,max_digits}
    for n in range(1,max_digits+1):
        digit=(inv@(cur%P))%P;v=cur-M@digit;assert not np.any(v%P);cur=v//P
        xs=[u+powers*j for u,j in zip(xs,digit.tolist())];powers*=P
        if n not in expanded:continue
        # A common denominator can be exposed by one generic linear combination.
        # It is only a guess. Equality on ALL original coordinates decides acceptance.
        vals=probes+[None];den=1;seen=[]
        for j in vals:
            live=(sum((1+i%251)*a for i,a in enumerate(xs)) if j is None else xs[j])%powers
            q=jc.reconstruct(live,powers)
            if q is None:continue
            if q.denominator not in seen:seen.append(q.denominator)
        for d in seen:
            if d!=1 and (den%d==0): continue
            den=lcm(den,d)
            if den*den*2>powers:break
            nums=[(x*den)%powers for x in xs];nums=[u if 2*u<=powers else u-powers for u in nums]
            if max(map(abs,nums),default=0)>isqrt(powers//2):continue
            if accept(nums,den):
                progress({'digits':n,'verified':True,'denominator_bits':den.bit_length(),'max_numerator_bits':max(abs(x).bit_length() for x in nums)})
                return nums,den
        progress({'digits':n,'verified':False,'carry_bound':int(np.max(np.abs(cur)))})
    return None

def run(previous,frames,out,top,max_digits,alt=False):
    z,A,b,pool,keys=jc.assemble(previous,'internal',frames);assert A.shape==(9249,5475)
    # Choose sparse, small signed coordinate rows before pivoting. The equations and
    # their labels are preserved; this changes only the chosen square minor.
    def order(B):
        if not alt:return np.argsort(np.diff(B.indptr),kind='stable')
        mass=np.asarray(abs(B).sum(axis=1)).ravel();return np.lexsort((np.arange(B.shape[0]),np.diff(B.indptr),mass))
    nr=order(A);nc=order(A.T.tocsr());c=[int(nc[j]) for j in jc.piv(A[nr,:][:,nc])];rr=[int(nr[i]) for i in jc.piv(A[nr,:][:,c].T.tocsr())];r=len(rr);assert r==len(c)==4308
    B=A[rr,:].T.tocsr();M=B[c,:].tocsr();small=[b[j] for j in c];out.mkdir(exist_ok=True,parents=True);history=[]
    def progress(rec):history.append(rec);jc.s.write(out/'precision.json',history);print(rec,flush=True)
    ans=rational_digit_solution(M,small,lambda x,d:numerator_check(B,x,b,d),max_digits,progress)
    if ans is None:return
    nums,den=ans;assert numerator_check(B,nums,b,den);div=den*z['denominator'];co={}
    for i,x in zip(rr,nums):
        if x:jc.l.p.axpy(co,jc.E.projected(jc.s.unpack_hex(pool[i][1])),Fraction(x,div))
    assert jc.side.even(co)==co;start=jc.side.even(jc.l.load_vector(top));assert jc.side.canonical_digest([[jc.s.hexes(w),str(v)] for w,v in sorted(start.items())])==z['base'];jc.l.p.axpy(start,co,1);jc.l.save_vector(out/'vector.json.gz',start)
    # Single denominator on disk avoids repeating a long integer for each word.
    whole=lcm(*(v.denominator for v in start.values()));lst=[[jc.s.hexes(w),str(int(v*whole))] for w,v in sorted(start.items())];jc.s.write(out/'vector_shared.json.gz',dict(denominator=str(whole),words=lst));idx=hashlib.sha256(np.asarray(rr+c,dtype='<i8').tobytes()).hexdigest();jc.s.write(out/'exact.json',dict(prime=P,scopes='internal',parents=len(pool),target_coordinates=len(keys),pivot_order=r,pivot_sha256=idx,equations_with_integer_residual=0,rhs_sha=z['rhs'],base=z['base'],vector_sha=jc.side.canonical_digest(jc.s.read(out/'vector.json.gz')),vector_words=len(start),delta_words=len(co),denominator_bits=whole.bit_length(),digits=history[-1]['digits'],alternate_order=alt))

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--previous',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--top',type=Path,required=True);a.add_argument('--digits',type=int,default=4096);a.add_argument('--alternate-order',action='store_true');v=a.parse_args();assert __debug__ and v.digits>=4;run(v.previous,v.previous/'frames',v.out,v.top,v.digits,v.alternate_order)
