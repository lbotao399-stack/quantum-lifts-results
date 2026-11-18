"""Integer evaluation of SO(N) traces, with explicit int64 overflow bounds.

Input matrices have entries in {-1,0,1}; each Grassmann coefficient of
a length-L multitrace with d odd letters has absolute value <= d! N^L.
Expression accumulation uses Python integers. No modular reduction occurs.
"""
from highest import *
from fortuity.hwcoh.evalg import table
from math import factorial,isqrt,gcd
from collections import OrderedDict


def rational_reconstruct(value,p):
    value=int(value)%p
    if not value:return Fraction(0)
    bound=isqrt(p//2);r0,r1=p,value;t0,t1=0,1
    while abs(r1)>bound:
        q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
    if not t1 or abs(t1)>bound or gcd(r1,t1)!=1:raise ValueError((value,p))
    f=Fraction(r1,t1)
    if (f.numerator-value*f.denominator)%p:raise ValueError((value,p))
    return f


def rational_hw(sec,p=1000003):
    from hw_modular import modular_hw_basis
    result=[]
    for e in modular_hw_basis(sec,0,p):
        fractions={w:rational_reconstruct(c,p) for w,c in e.items()}
        den=lcm(*(c.denominator for c in fractions.values()))
        r={w:int(c*den) for w,c in fractions.items()}
        for i in (1,2):
            check=defaultdict(int)
            for w,c in r.items():
                for v,d in raising(w,i).items():check[v]+=c*d
            assert not any(check.values()),'Rational highest-weight reconstruction failed'
        result.append(r)
    return result


class IntegerFrame:
    def __init__(self,n,npts,k,seed):
        self.N,self.P,self.k=n,npts,k
        rng=np.random.default_rng(seed);self.L={}
        for l in (0,1,2,6,3,4,5):
            d=int(3<=l<=5);M=rng.integers(-1,2,(npts,k if d else 1,n,n),dtype=np.int64)
            M=np.triu(M,1);self.L[l]=(d,M-M.swapaxes(-1,-2))
        self.matrices=OrderedDict();self.matrix_bytes=0;self.scalars={}

    def matrix(self,w):
        d=sum(3<=l<=5 for l in w)
        assert factorial(d)*self.N**len(w)<2**62
        if w in self.matrices:
            self.matrices.move_to_end(w);return self.matrices[w]
        if len(w)==1:d,v=self.L[w[0]]
        else:
            A,dA=self.matrix(w[:-1]);dl,B=self.L[w[-1]]
            if not dl:v=A@B;d=dA
            else:
                I,J,S,st=table(self.k,dA,1)
                R=(A[:,I]@B[:,J])*S[None,:,None,None]
                v=np.add.reduceat(R,st,axis=1);d=dA+1
        while self.matrices and self.matrix_bytes+v.nbytes>32*1024**2:
            _,(old,_)=self.matrices.popitem(last=False);self.matrix_bytes-=old.nbytes
        if v.nbytes<=32*1024**2:self.matrices[w]=(v,d);self.matrix_bytes+=v.nbytes
        return v,d

    def trace(self,w):
        if w not in self.scalars:
            v,d=self.matrix(w);self.scalars[w]=(np.trace(v,axis1=-2,axis2=-1),d)
        return self.scalars[w]

    def invariant(self,w):
        paths,traces=w;assert not paths
        total=sum(map(len,traces));dtotal=sum(3<=l<=5 for t in traces for l in t)
        assert factorial(dtotal)*self.N**total<2**62
        A=np.ones((self.P,1),dtype=np.int64);d=0
        for t in traces:
            B,e=self.trace(t);I,J,S,st=table(self.k,d,e)
            A=np.add.reduceat(A[:,I]*B[:,J]*S[None,:],st,axis=1);d+=e
        return A.reshape(-1)


def integer_values(exprs,degrees_points,n=8,seed=2711):
    """degrees_points: {fermion_degree:(k,npoints)}; disjoint column blocks."""
    blocks=[];controls=[]
    for degree,(k,npoints) in sorted(degrees_points.items()):
        words=sorted({w for e in exprs for w in e if sum(count(w)[3:6])==degree})
        idx={w:i for i,w in enumerate(words)}
        coefficients=np.zeros((len(exprs),len(words)),dtype=object)
        for i,e in enumerate(exprs):
            for w,c in e.items():
                if w in idx:coefficients[i,idx[w]]=int(c)
        # A priori bound permits sparse int64 multiplication only when safe.
        bounds=[factorial(degree)*n**sum(map(len,w[1])) for w in words]
        safe=max((sum(abs(int(c))*b for c,b in zip(row,bounds)) for row in coefficients),default=0)<2**62
        if safe:
            sparse=csr_matrix(np.asarray(coefficients,dtype=np.int64))
        pieces=[]
        for start in range(0,npoints,8):
            fr=IntegerFrame(n,min(8,npoints-start),k,seed+100003*degree+7919*start)
            values=np.array([fr.invariant(w) for w in words],dtype=np.int64)
            pieces.append((sparse@values) if safe else coefficients@values.astype(object))
        blocks.append(np.concatenate(pieces,axis=1));controls.append(dict(degree=degree,k=k,npoints=npoints,int64_expression_safe=safe))
    return np.concatenate(blocks,axis=1),controls
