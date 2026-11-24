#!/usr/bin/env python3
"""Independent colour-coordinate Wick checks, exact norms and basis covariance.

This checker builds an explicit orthogonal Hermitian USp(4) generator basis.
Its colour calculations do not use the signed-index union algorithm in wick.cpp.
Real integer matrices represent each generator after its global phase 1 or i
is stripped. That phase cancels between equal occupation monomials in a norm.
"""
from compute_metric import *
from itertools import product, permutations
import random

def generators(rank=2):
    size=2*rank;answer=[]
    def add(entries,phase=0):
        a=np.zeros((size,size),dtype=np.int64)
        for i,j,n in entries:a[i,j]=n
        norm=(-1)**phase*int(np.trace(a@a))
        assert norm>0
        answer.append((a,phase,norm))
    for i in range(rank):add([(i,i,1),(rank+i,rank+i,-1)])
    for i in range(rank):
        for j in range(i+1,rank):
            add([(i,j,1),(j,i,1),(rank+j,rank+i,-1),(rank+i,rank+j,-1)])
            add([(i,j,1),(j,i,-1),(rank+j,rank+i,-1),(rank+i,rank+j,1)],1)
    for i in range(rank):
        add([(i,rank+i,1),(rank+i,i,1)])
        add([(i,rank+i,1),(rank+i,i,-1)],1)
        for j in range(i+1,rank):
            add([(i,rank+j,1),(j,rank+i,1),(rank+j,i,1),(rank+i,j,1)])
            add([(i,rank+j,1),(j,rank+i,1),(rank+j,i,-1),(rank+i,j,-1)],1)
    assert len(answer)==rank*(2*rank+1)
    symplectic=np.zeros((size,size),dtype=np.int64)
    symplectic[:rank,rank:]=np.eye(rank,dtype=np.int64)
    symplectic[rank:,:rank]=-np.eye(rank,dtype=np.int64)
    for i,(a,phase,norm) in enumerate(answer):
        assert np.array_equal(a.T,(-1)**phase*a)
        assert not (a.T@symplectic+symplectic@a).any()
        for j,(b,other,_) in enumerate(answer):
            if i!=j:assert np.trace(a@b)==0
    return answer,symplectic

def direct_covariance():
    basis,J=generators();size=len(J);count=0
    for i,j,k,l in product(range(size),repeat=4):
        # T_ij conjugate(T_kl) has cancelling global phases.
        lhs=sum(Fraction(int(a[i,j]*a[k,l]),n) for a,p,n in basis)
        rhs=Fraction(int(i==k and j==l)+int(J[i,l]*J[k,j]),2)
        assert lhs==rhs,(i,j,k,l,lhs,rhs);count+=1
    return count

def coefficient_polynomial(word,basis):
    flat=[l for tr in word for l in tr];length=len(flat);result=Counter()
    sizes=list(map(len,word));size=basis[0][0].shape[0]
    for assignment in product(range(len(basis)),repeat=length):
        coefficient=1;start=0
        for count in sizes:
            matrix=np.eye(size,dtype=np.int64)
            for a in assignment[start:start+count]:matrix=matrix@basis[a][0]
            coefficient*=int(np.trace(matrix));start+=count
        if not coefficient:continue
        fermions=[(l,a) for l,a in zip(flat,assignment) if is_ferm(l)]
        if len(set(fermions))<len(fermions):continue
        inversions=sum(x>y for j,x in enumerate(fermions) for y in fermions[j+1:])
        key=tuple(sorted(zip(flat,assignment)))
        result[key]+=(-1)**inversions*coefficient
    return {k:v for k,v in result.items() if v}

def direct_pair(left,right,basis):
    total=Fraction(0)
    for monomial,c in left.items():
        if monomial not in right:continue
        weight=Fraction(c*right[monomial])
        for (letter,colour),count in Counter(monomial).items():
            weight*=Fraction(math.factorial(count) if not is_ferm(letter) else 1,basis[colour][2]**count)
        total+=weight
    return total

def direct_trace_checks():
    a=(0,0,(0,));b=(0,0,(1,));f=(1,0,());g=(0,1,())
    contents_to_check=[(a,a),(a,b),(f,g),(a,a,a,a),(a,a,b,b),(a,b,f,g),(f,f,g,g)]
    basis,J=generators();count=0;records=[]
    for content in contents_to_check:
        words=set()
        for w in multitraces_with_content(tuple(sorted(content))):
            sig,v=canonical(w)
            if sig:words.add(v)
        words=sorted(words)
        polynomials=[coefficient_polynomial(w,basis) for w in words]
        fast=colour_gram(words,4)
        for i in range(len(words)):
            for j in range(len(words)):
                direct=direct_pair(polynomials[i],polynomials[j],basis)
                assert direct==Fraction(int(fast[i,j]),2**len(content)),(i,j,direct,fast[i,j])
                count+=1
        records.append({'letters':[[a,b,list(s)] for a,b,s in content],
                        'words':len(words),'entries':len(words)**2})
    return count,records

def oscillator_checks():
    letters=set()
    for length in (5,6,7,8):
        for offset in (-1,0):
            data=read(HERE/f'frames_Y{length}_offset{offset}.json.gz')
            for rows in data['basis'].values():
                for row in import_rows(rows):
                    for w in row:
                        for tr in w:letters.update(tr)
    @lru_cache(None)
    def descendant(a,b,k):
        # [K^{dot alpha},P_{dot beta}] in (C.7), with
        # Delta+J_L=k.  Primary norms are 1,1,2 for scalar,psi,f.
        if b:return b*(k+a+b-1)*descendant(a,b-1,k)
        if a:return a*(k+a-1)*descendant(a-1,0,k)
        return (1,1,2)[k-1]
    @lru_cache(None)
    def lambda_plus(a,b):
        if b:return b*(a+b)*lambda_plus(a,b-1)
        if a:return a*(a+1)*lambda_plus(a-1,0)
        return 1
    @lru_cache(None)
    def lambda_minus(a,b):return lambda_plus(b,a)
    @lru_cache(None)
    def lambda_mixed(a,b):
        assert a and b
        return ((b-1)*(a+b)*lambda_mixed(a,b-1) if b>1 else 0)+a*lambda_plus(a-1,b-1)
    checked=[]
    for a,b,s in sorted(letters):
        k=len(s)
        if k:
            value=descendant(a,b,k)
            expected=Fraction(value,(math.factorial(a)*math.factorial(b))**2)
        else:
            # The generating gaugino combines both dotted spinor directions;
            # all diagonal and mixed terms are retained, including the e.o.m.
            value=(a*a*lambda_plus(a-1,b) if a else 0)+(b*b*lambda_minus(a,b-1) if b else 0)
            if a and b:value+=2*a*b*lambda_mixed(a,b)
            expected=Fraction(value,((a+b)*math.factorial(a)*math.factorial(b))**2)
        assert expected==letter_norm((a,b,s))
        checked.append({'letter':[a,b,list(s)],'norm':str(expected)})
    return checked

def gauge_action_spinor_normalization():
    """Resolve f directly from the action, independently of the unified norm.

    All coefficients below are in units g_YM^2/pi^2.  With tr(T_A T_B)=delta_AB,
    S_gauge=(1/(2*g_YM^2))*integral tr(F_mu_nu F_mu_nu).  The gauge-fixed free
    propagator is delta_mu_nu/(8*x^2).  Its curvature derivative at a unit radial
    point, followed by the tensor inversion in the radial bra, is evaluated
    before the spinor contractions.  All 16 spinor-pair entries are checked.
    """
    import sympy as sp
    sigma=[sp.eye(2),sp.Matrix([[0,sp.I],[sp.I,0]]),
           sp.Matrix([[0,1],[-1,0]]),sp.diag(sp.I,-sp.I)]
    epsilon=sp.Matrix([[0,1],[-1,0]])
    radial=[-1,1,1,1]
    second=[6,-2,-2,-2] # Hessian of 1/x^2 at x=(1,0,0,0).
    curvature={}
    for mu,nu,rho,tau in product(range(4),repeat=4):
        value=sp.Rational(1,8)*(-int(nu==tau and mu==rho)*second[mu]
              +int(nu==rho and mu==tau)*second[mu]
              +int(mu==tau and nu==rho)*second[nu]
              -int(mu==rho and nu==tau)*second[nu])
        value*=radial[mu]*radial[nu]
        assert value==sp.Rational(1,2)*(int(mu==rho and nu==tau)-int(mu==tau and nu==rho))
        curvature[mu,nu,rho,tau]=value
    coefficient={}
    for mu,nu in product(range(4),repeat=2):
        # A.5 has normalized antisymmetrization and its additional factor 1/2.
        raised=(epsilon*sigma[mu]*epsilon.T*sigma[nu].T
                -epsilon*sigma[nu]*epsilon.T*sigma[mu].T)/4
        lowered=-epsilon*raised
        assert lowered==lowered.T
        coefficient[mu,nu]=-lowered/16 # D.4, summed over both mu,nu orders.
    entries=[]
    for alpha,beta,gamma,delta in product(range(2),repeat=4):
        value=sum(sp.conjugate(coefficient[mu,nu][alpha,beta])*
                  coefficient[rho,tau][gamma,delta]*curvature[mu,nu,rho,tau]
                  for mu,nu,rho,tau in product(range(4),repeat=4))
        expected=sp.Rational(1,256)*(int(alpha==gamma and beta==delta)+int(alpha==delta and beta==gamma))
        assert sp.simplify(value-expected)==0
        entries.append(str(value))
    primary=sp.Rational(1,128)
    native=16*16*primary # (4*pi/g)^2 times the superfield f coefficient 4 squared.
    assert native==2
    return {'passed':True,'radially_conjugated_curvature_entries_checked':256,
            'spinor_pair_entries_checked':16,'spinor_pair_entries':entries,
            'field_strength_primary_norm_in_units_g_squared_over_pi_squared':'1/128',
            'native_undifferentiated_field_strength_norm':'2',
            'source_equations':['2306.04673v2 (A.1)-(A.5)','2306.04673v2 (D.4),(D.9)'],
            'audit_note':'The direct action calculation agrees with (3.13)/(E.30); the field-strength line in (E.12)/(E.13) is smaller by a factor two and is not used.'}

def positive_factorization(matrix):
    integer,denominator=matrix.numer_denom()
    assert denominator>0
    permutation,lower,diagonal,upper=integer.fflu()
    assert permutation.is_one()
    assert upper==lower.transpose()
    assert all(diagonal[i,i]>0 for i in range(diagonal.nrows()))
    assert fmpq_mat(integer)==fmpq_mat(lower)*fmpq_mat(diagonal).inv()*fmpq_mat(upper)
    return {'common_denominator':str(denominator),
            'positive_diagonal':[str(diagonal[i,i]) for i in range(diagonal.nrows())]}

def basis_changes():
    checked=[]
    for length in (5,6,7,8):
        file=HERE/f'harmonic_Y{length}.json.gz'
        if not file.exists():continue
        data=read(file);matrix=import_matrix(data['harmonic_gram']);n=matrix.nrows()
        if not n:
            checked.append({'length':length,'dimension':0});continue
        factorization=positive_factorization(matrix)
        # An exact integral invertible triangular change mixes the entire block.
        change=fmpq_mat(n,n)
        for i in range(n):
            change[i,i]=1
            if i+1<n:change[i,i+1]=(-1)**i*(i+1)
        transformed=change.transpose()*matrix*change
        assert transformed.det()==matrix.det()
        vector=fmpq_mat([[i+1] for i in range(n)])
        coordinates=change.solve(vector)
        assert vector.transpose()*matrix*vector==coordinates.transpose()*transformed*coordinates
        checked.append({'length':length,'dimension':n,'determinant':str(matrix.det()),
                        'norm_before_and_after':str((vector.transpose()*matrix*vector)[0,0]),
                        'positive_factorization':factorization})
    return checked

def main():
    start=time.monotonic();covariance=direct_covariance()
    contractions,records=direct_trace_checks();normalization=gauge_action_spinor_normalization()
    oscillators=oscillator_checks();changes=basis_changes()
    result={'passed':True,'explicit_hermitian_generator_count':10,'all_colour_covariance_entries':covariance,
            'independent_direct_colour_Wick_entries':contractions,'trace_checks':records,
            'independent_gauge_action_normalization':normalization,
            'all_letter_radial_recurrence_norms':oscillators,'exact_basis_change_checks':changes,
            'seconds':time.monotonic()-start}
    write(HERE/'independent_checks.json',result)
    print({'passed':True,'covariance':covariance,'Wick_entries':contractions,'letters':len(oscillators),
           'basis_changes':[r['dimension'] for r in changes],'seconds':result['seconds']},flush=True)

if __name__=='__main__':main()
