#!/usr/bin/env python3
"""Certify the entire reconstructed rational trace-identity space by bounded CRT.

Every prime uses the same small integer symplectic matrices. One complete
finite-colour frame is certified by a nonzero2689minor. The exact integer
relation evaluations are zero because all prime residues vanish and their
product exceeds twice the proven absolute bound.
"""
from probe_stable import HERE,m,np,nmod_mat
from fractions import Fraction
from collections import Counter
import math,pickle,json,gzip,time,hashlib
import dual_charge_bands as bands
from terminal_interval_ranks import IntegralTerminalFrame
from flint import fmpz

BASE_PRIME=1000003

def reconstruct(value):
    if not value:return Fraction(0)
    r,s=BASE_PRIME,int(value);t,u=0,1;bound=math.isqrt((BASE_PRIME-1)//2)
    while abs(s)>bound:
        quotient=r//s;r,s=s,r-quotient*s;t,u=u,t-quotient*u
    assert u and abs(u)<=bound and math.gcd(s,u)==1
    result=Fraction(s,u)
    assert (result.numerator-value*result.denominator)%BASE_PRIME==0
    return result

def integer_frame(colours,group,points,degree,alphabet,seed):
    assert colours==6 and group=='Sp'
    frame=m.Frame(3,points,degree+1,[degree],alphabet,seed=seed,p=BASE_PRIME,group='Sp')
    j=np.block([[np.zeros((3,3),dtype=np.int64),np.eye(3,dtype=np.int64)],[-np.eye(3,dtype=np.int64),np.zeros((3,3),dtype=np.int64)]])
    integer={}
    for letter,(parity,values) in frame.L.items():
        centred=np.where(values>BASE_PRIME//2,values-BASE_PRIME,values)
        small=np.sign(centred)*(np.abs(centred)%3)
        assert np.max(np.abs(small))<=2
        assert not np.any(small.swapaxes(-1,-2)@j+j@small)
        integer[letter]=(parity,small)
        frame.L[letter]=(parity,small%m.P)
    frame.p=m.P
    frame.integer_letters=integer
    return frame

def main():
    start=time.monotonic()
    _,centre,_,_=pickle.loads((HERE/'stable_rows.pkl').read_bytes())
    stored=np.load(HERE/'relation_image.npz');residues=stored['relations'];basis=stored['basis_rows']
    lookup={int(v):reconstruct(int(v)) for v in np.unique(residues)}
    relations=[];denominators=[];polynomials=[]
    for row in residues:
        co={int(i):lookup[int(row[i])] for i in np.flatnonzero(row)}
        denominator=math.lcm(*(v.denominator for v in co.values()))
        poly=Counter()
        for i,v in co.items():
            for w,n in centre[i].items():poly[w]+=int(v*denominator)*n
        relations.append(co);denominators.append(denominator)
        polynomials.append({w:n for w,n in poly.items() if n})
    maximum=max(sum(abs(c) for c in row.values()) for row in polynomials)
    for row in polynomials:
        assert all(m.nferm(w)==4 and sum(map(len,w))==9 for w in row)
    bound=maximum*math.factorial(4)*12**9
    primes=[1000003,1000033,1000037]
    assert all(fmpz(prime).is_prime() for prime in primes)
    assert math.prod(primes)>2*bound
    page=json.loads(gzip.decompress((HERE/'inputs/algorithms/dual_singlets_20261118/sp3_omega/target.json.gz').read_bytes()))
    letters=sorted({l for content in m.contents(tuple(page['weight']),9) for l in content})
    samples=page['centre_frames']['4']['samples']
    bands.dense_frame=integer_frame
    checks=[];full_determinant=None
    words=sorted({w for row in polynomials for w in row});controls=[words[i] for i in np.linspace(0,len(words)-1,12,dtype=int)]
    for prime in primes:
        m.P=prime
        lib=m.library(HERE,6,prime)
        rows=([centre[i] for i in basis]+polynomials) if prime==BASE_PRIME else polynomials
        values=bands.coordinates(rows,4,letters,samples,6,'Sp',lib,1)
        if prime==BASE_PRIME:
            full_determinant=int(nmod_mat(values[:len(basis)].tolist(),prime).det())
            assert full_determinant
            values=values[len(basis):]
        assert not np.any(values)
        points,seed,_=samples[0]
        frame=integer_frame(6,'Sp',points,4,letters,seed)
        exact=IntegralTerminalFrame(6,points,5,[4],letters,seed=1)
        exact.L=frame.integer_letters
        exact_values=exact.eval_words(controls,dtype=np.int64)
        np.testing.assert_array_equal(m.values(controls,frame,lib),exact_values%prime)
        checks.append({'prime':prime,'zero_coordinates':int(values.size),'independent_integer_controls':int(exact_values.size)})
        print('prime',prime,'exact-bound zeros',values.size,'seconds',time.monotonic()-start,flush=True)
    out={'formal_dimension':len(centre),'finite_dimension':len(basis),'complete_relation_dimension':len(relations),'basis_coefficients':[[[i,str(c)] for i,c in sorted(row.items())] for row in relations],'integer_row_denominators':denominators,'free_rows':stored['free_rows'].tolist()}
    raw=(json.dumps(out,separators=(',',':'))+'\n').encode()
    (HERE/'rational_relation_basis.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    result={'passed':True,'source_basis_dimension':len(centre),'finite_character_dimension':len(basis),'rational_relation_dimension':len(relations),'reconstructed_nonzero_coefficients':int(np.count_nonzero(residues)),'distinct_reconstructed_coefficients':len(lookup),'denominators':sorted(set(denominators)),'total_integer_polynomial_terms':sum(map(len,polynomials)),'maximum_integer_polynomial_coefficient_l1':maximum,'per_word_absolute_bound':math.factorial(4)*12**9,'integer_evaluation_absolute_bound':bound,'crt_modulus':math.prod(primes),'full_finite_frame_minor_determinant':full_determinant,'full_finite_frame_minor_prime':BASE_PRIME,'checks':checks,'rational_basis_sha256':hashlib.sha256(raw).hexdigest(),'seconds':time.monotonic()-start,'scope':'All181 rational finite-colour trace identities certified exactly on a character-saturated2689-coordinate frame. Bounded CRT uses identical integer Sp3 matrix points at allthree primes, including all cleared denominators and full highest-weight coefficients.'}
    (HERE/'rational_relations_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
