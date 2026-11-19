"""Bounded rational reconstruction of a specific pair of coupled modular solutions.

Each accepted primitive fraction satisfies both congruences. Uniqueness follows
from 2*bound**2 < modulus. For a solution n/d with 0<d<=bound and |n|<=bound,
|x/modulus-k/d| = |n|/(modulus*d) < 1/(2*d*d) for k=(d*x-n)/modulus.
Since gcd(d,modulus)=gcd(d,n)=1, k/d is reduced and a continued-fraction
convergent of x/modulus (including 0/1). Testing all convergents with d<=bound
therefore certifies failures too. No bound on an unrestricted correction is inferred.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from io import BytesIO
from math import gcd, isqrt
from pathlib import Path
from zipfile import ZipFile
import argparse, gzip, json
import numpy as np
from verify_su3_j30_solution_receipts import verify


def convergents(x,modulus,bound):
    a,b=0,1;c,d=1,0
    s,t=x,modulus
    while t:
        h,r=divmod(s,t); a,c=c,h*c+a; b,d=d,h*d+b
        if d>bound:return
        yield c,d
        s,t=t,r


def rational(x,m,bound):
    assert 0<=x<m and 2*bound**2<m
    for k,d in convergents(x,m,bound):
        n=d*x-k*m
        if abs(n)<=bound and gcd(n,d)==gcd(d,m)==1:return n,d
    return None


def control():
    for m in (5,15,35,77,143,187):
        bound=isqrt((m-1)//2)
        table={n*pow(d,-1,m)%m:(n,d) for d in range(1,bound+1) if gcd(m,d)==1 for n in range(-bound,bound+1) if gcd(d,n)==1}
        assert all(rational(x,m,bound)==table.get(x) for x in range(m))


def run(data):
    control();p,q=1000003,1000033;m=p*q;B=isqrt((m-1)//2)
    path=data/'N4_SU3_j30_solution_receipts_20261118.zip';summary=verify(path)
    meta=json.loads((data/'N4_SU3_j30_first_prime_residuals_20261118.json').read_bytes());other=summary['coupled_equations']
    left=gzip.decompress((data/'N4_SU3_j30_first_prime_solution.npy.gz').read_bytes());assert sha256(left).hexdigest()==meta['solution_sha256'] and meta['prime']==p
    with ZipFile(path) as archive:
        right=archive.read('second_solve_compact/particular.npy');free=archive.read('second_solve_compact/free_variables.npy')
    assert sha256(right).hexdigest()==other['particular_sha256'] and sha256(free).hexdigest()==meta['free_sha256']==other['free_variables_sha256']
    left,right,free=(np.load(BytesIO(x),allow_pickle=False) for x in (left,right,free));assert left.shape==right.shape==(77639,) and len(free)==27492 and (left<p).all() and (right<q).all() and not (left[free].any() or right[free].any())
    inv=pow(p,-1,q);ok=Counter();bad=Counter();first=None;finger=sha256()
    for i,(a,b) in enumerate(zip(left,right)):
        a,b=int(a),int(b);x=a+p*((b-a)*inv%q);rec=rational(x,m,B);block='first_correction' if i<28311 else 'second_correction'
        if rec is None:
            bad[block]+=1;finger.update(b'-\n')
            if first is None:
                candidates=[(abs(d*x-k*m),d) for k,d in convergents(x,m,B)]
                # A second, exhaustive certificate, independent of continued fractions.
                best=min((min(d*x%m,(-d*x)%m),d) for d in range(1,B+1) if gcd(d,m)==1)
                assert min(candidates)==best and best[0]>B
                first=dict(index=i,residues=[a,b],CRT_residue=x,minimum_absolute_numerator=best[0],denominator_attaining_minimum=best[1])
        else:
            n,d=rec;assert (n-a*d)%p==(n-b*d)%q==0;ok[block]+=1;finger.update(f'{n}/{d}\n'.encode())
    assert (sum(ok.values()),sum(bad.values()))==(58230,19409)
    return dict(primes=[p,q],modulus=m,equal_numerator_denominator_bound=B,free_coordinates_set_to_zero=len(free),ordered_pivots_in_agreement=summary['matching_ordered_pivot_lists'],reconstructed=dict(ok),unreconstructed=dict(bad),first_failure=first,fractions_and_failure_markers_sha256=finger.hexdigest(),rational_correction_excluded=False)

if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__);args.add_argument('--output',type=Path);args=args.parse_args();assert __debug__
    text=json.dumps(run(Path(__file__).resolve().parents[1]/'data'),indent=2)+'\n';print(text,end='')
    if args.output:args.output.write_text(text)
