"""Reconstruct only after matching ordered free coordinates and replay receipts at each prime.
Usage: python postsolve.py REPOSITORY_DATA THIRD_SECOND_SOLVE_DIR OUTPUT_JSON
The result is a congruence test, not a rational residual certificate.
"""
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
from math import gcd,isqrt
from hashlib import sha256
from collections import Counter
import argparse,gzip,json,numpy as np
P=(1000003,1000033,1000037)
FREE='a7421e462000b93c18376fbe5eda38c69e1cfa2b8642a7e5034e1824309b0e2d'

def rational(x,m,b):
    a,v=0,1;c,d=1,0;s,t=x,m
    while t:
        h,r=divmod(s,t);a,c=c,h*c+a;v,d=d,h*d+v
        if d>b:break
        n=d*x-c*m
        if abs(n)<=b and gcd(d,n)==gcd(d,m)==1:return (n,d)
        s,t=t,r
    return None

def control():
    for m in (5,15,35,77,143):
        b=isqrt((m-1)//2);table={n*pow(d,-1,m)%m:(n,d) for d in range(1,b+1) if gcd(m,d)==1 for n in range(-b,b+1) if gcd(d,n)==1}
        assert all(rational(x,m,b)==table.get(x) for x in range(m))

def run(data,third):
    control();p,q,r=P;m=p*q*r;b=isqrt((m-1)//2);b0=isqrt((p*q-1)//2);H=lambda x:sha256(x).hexdigest()
    a=gzip.decompress((data/'N4_SU3_j30_first_prime_solution.npy.gz').read_bytes());one=json.loads((data/'N4_SU3_j30_first_prime_residuals_20261118.json').read_bytes());assert one['prime']==p and H(a)==one['solution_sha256'] and one['free_sha256']==FREE
    with ZipFile(data/'N4_SU3_j30_solution_receipts_20261118.zip') as z:
        get=lambda n:z.read('second_solve_compact/'+n);bb=get('particular.npy');fb=get('free_variables.npy');two=json.loads(get('result.json'))
    get=lambda n:(third/n).read_bytes();cc=get('particular.npy');fc=get('free_variables.npy');three=json.loads(get('result.json'))
    for prime,x,f,s in ((q,bb,fb,two),(r,cc,fc,three)):
        assert s['prime']==prime and s['status']=='complete_modular_second_solution_with_all_first_parameters' and s['equation_rank']==50147 and s['variables']==77639 and s['homogeneous_dimension']==27492
        assert len(s['residuals'])==50 and all(e['nonzero']==0 for e in s['residuals']) and [sum(e['coordinates'] for e in s['residuals'] if e['equation']==i) for i in (1,2)]==[34901,65406]
        assert H(f)==s['free_variables_sha256']==FREE and H(x)==s['particular_sha256']
    assert fb==fc;load=lambda v:np.load(BytesIO(v),allow_pickle=False);f=load(fb);v=[load(x) for x in (a,bb,cc)];assert f.shape==(27492,) and all(w.shape==(77639,) and (w<prime).all() and not w[f].any() for w,prime in zip(v,P))
    inv1=pow(p,-1,q);inv2=pow(p*q,-1,r);stats=[Counter(),Counter()];stable=Counter();finger=sha256()
    for i,values in enumerate(zip(*v)):
        u,h,k=map(int,values);z=u+p*((h-u)*inv1%q);prev=rational(z,p*q,b0);candidate=rational(z+p*q*((k-z)*inv2%r),m,b);block=int(i>=28311)
        if prev is not None:stable['agrees' if (prev[0]-k*prev[1])%r==0 else 'disagrees']+=1
        if candidate is None:stats[block]['unreconstructed']+=1;finger.update(b'-\n');continue
        n,d=candidate;assert all((n-int(x)*d)%prime==0 for x,prime in zip(values,P));stats[block]['reconstructed']+=1;finger.update(f'{n}/{d}\n'.encode())
    assert sum(stable.values())==58230
    return dict(primes=P,modulus=m,bound=b,verified_zero_free_coordinates=len(f),prior_two_prime_candidates=stable,first_correction=stats[0],second_correction=stats[1],fraction_marker_sha256=finger.hexdigest(),rational_residual_verified=False)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('data',type=Path);p.add_argument('third',type=Path);p.add_argument('output',type=Path);opts=p.parse_args();assert __debug__;opts.output.write_text(json.dumps(run(opts.data,opts.third),indent=2)+'\n')
