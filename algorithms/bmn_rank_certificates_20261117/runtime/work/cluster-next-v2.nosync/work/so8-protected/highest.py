"""Exact rational SU(3) highest-weight vectors on SO(8) invariant graphs."""
from protected import *
from flint import fmpq_mat,nmod_mat
from fractions import Fraction
from math import lcm

EMAP={1:{1:(1,0),3:(-1,4)},2:{2:(1,1),4:(-1,5)}}


def raising(w,i):
    out=defaultdict(int);paths,traces=w
    for typ,seq in enumerate((paths,traces)):
        for k,word in enumerate(seq):
            for pos,l in enumerate(word):
                if l not in EMAP[i]:continue
                c,nl=EMAP[i][l];new=seq[:k]+(word[:pos]+(nl,)+word[pos+1:],)+seq[k+1:]
                s,v=canon(new,traces) if typ==0 else canon(paths,new)
                if s:out[v]+=s*c
    return {w:c for w,c in out.items() if c}


def rational_left_kernel(A):
    A=np.asarray(A,dtype=object);m,n=A.shape
    if n==0:return np.eye(m,dtype=object)
    R,rank=fmpq_mat([[int(A[i,k]) for i in range(m)] for k in range(n)]).rref()
    pivots=[next(i for i in range(m) if R[k,i]) for k in range(rank)]
    rows=[]
    for free in range(m):
        if free in pivots:continue
        row=[Fraction(0) for _ in range(m)];row[free]=Fraction(1)
        for k,pivot in enumerate(pivots):row[pivot]=-Fraction(str(R[k,free]))
        den=lcm(*(x.denominator for x in row));rows.append([int(x*den) for x in row])
    return np.array(rows,dtype=object).reshape(-1,m)


def shape(inv):
    out=[]
    for typ,seq in enumerate(inv):
        ws=[]
        for w in seq:
            t=tuple(0 if x<3 else 3 if x<6 else 6 for x in w)
            if typ==0:ws.append(min(t,t[::-1]))
            else:ws.append(min(v[i:]+v[:i] for v in (t,t[::-1]) for i in range(len(v))))
        out.append(tuple(sorted(ws)))
    return tuple(out)


def hw_basis(sec,parity):
    groups=defaultdict(list)
    for w in sector_basis(sec,parity):groups[shape(w)].append(w)
    result=[]
    for words in groups.values():
        images=[{(i,v):c for i in (1,2) for v,c in raising(w,i).items()} for w in words]
        cols={v:k for k,v in enumerate(sorted({v for im in images for v in im}))}
        M=np.zeros((len(words),len(cols)),dtype=np.int64)
        for k,im in enumerate(images):
            for v,c in im.items():M[k,cols[v]]=c
        H=rational_left_kernel(M)
        for row in H:
            expr={w:int(c) for w,c in zip(words,row) if c}
            for i in (1,2):
                test=defaultdict(int)
                for w,c in expr.items():
                    for v,d in raising(w,i).items():test[v]+=c*d
                assert not any(test.values())
            result.append(expr)
    return result


def serialize(expr):
    return [[int(c),[list(x) for x in w[0]],[list(x) for x in w[1]]] for w,c in sorted(expr.items())]


def deserialize(terms):
    return {(tuple(tuple(x) for x in ps),tuple(tuple(x) for x in ts)):c for c,ps,ts in terms}


def main():
    for n in (0,1,2):
        basis=hw_basis((n,3,3,3),0)
        path=ROOT/('hw_j18_n%d.json'%n)
        path.write_text(json.dumps([serialize(e) for e in basis])+'\n')
        print('highest weight',n,len(basis),flush=True)


if __name__=='__main__':main()
