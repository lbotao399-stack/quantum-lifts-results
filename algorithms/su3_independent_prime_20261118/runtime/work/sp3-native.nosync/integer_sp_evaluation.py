"""Integer symplectic trace and Q0 evaluation with explicit overflow bounds.

The matrices use the pinned X=-Omega*S convention, with integral symmetric S
whose entries lie in {-1,0,1}. Q0 coefficients are imported without change.
"""
import sp_environment
import sp_census as census
from fortuity.hwcoh.evalg import table
from fortuity.bmn.words import QLET
from collections import OrderedDict
from math import factorial, comb
import numpy as np
from scipy.sparse import csr_matrix


def fermions(word):
    return sum(3 <= letter <= 5 for letter in word)


class IntegerSymplecticFrame:
    def __init__(self, rank, npts, k, degree, seed):
        self.rank, self.N, self.P, self.k, self.degree = rank, 2*rank, npts, k, degree
        self.ncols = npts*comb(k, degree)
        rng = np.random.default_rng(seed)
        omega = np.zeros((self.N,self.N), dtype=np.int64)
        omega[:rank,rank:] = np.eye(rank,dtype=np.int64)
        omega[rank:,:rank] = -np.eye(rank,dtype=np.int64)
        self.L = {}
        for letter in (0,1,2,6,3,4,5):
            odd = int(3 <= letter <= 5)
            upper = np.triu(rng.integers(-1,2,(npts,k if odd else 1,self.N,self.N),dtype=np.int64))
            symmetric = upper + np.swapaxes(upper,-1,-2)
            diagonal = np.arange(self.N)
            symmetric[...,diagonal,diagonal] = upper[...,diagonal,diagonal]
            matrix = -omega@symmetric
            assert np.max(np.abs(matrix),initial=0) <= 1
            assert not np.any(np.swapaxes(matrix,-1,-2)@omega + omega@matrix)
            self.L[letter] = odd,matrix

    def mat_mul(self,left,left_degree,right,right_degree):
        if not right_degree:
            return left@right
        assert right_degree == 1
        ii,jj,signs,starts = table(self.k,left_degree,1)
        raw = (left[:,ii]@right[:,jj])*signs[None,:,None,None]
        return np.add.reduceat(raw,starts,axis=1)

    def scal_mul(self,left,left_degree,right,right_degree):
        ii,jj,signs,starts = table(self.k,left_degree,right_degree)
        return np.add.reduceat(left[:,ii]*right[:,jj]*signs[None,:],starts,axis=1)


def trace_values(traces,frame,derivative=False):
    degree=frame.degree
    qletters={}
    if derivative:
        for letter,terms in QLET.items():
            parts={}
            for coefficient,word in terms:
                d=fermions(word)
                if d>degree:
                    continue
                left_degree,left=frame.L[word[0]]
                right_degree,right=frame.L[word[1]]
                raw=frame.mat_mul(left,left_degree,right,right_degree)*coefficient
                parts[d]=parts.get(d,0)+raw
            qletters[letter]=parts
    identity=np.broadcast_to(np.eye(frame.N,dtype=np.int64),(frame.P,1,frame.N,frame.N))
    ordinary,images={},{}
    previous,stack=(),[]
    for word in sorted(traces):
        length=len(word)
        assert factorial(degree)*frame.N**length*(6*length*frame.N if derivative else 1)<2**62
        common=0
        while common<min(len(word),len(previous)) and word[common]==previous[common]:
            common+=1
        del stack[common:]
        value,count,parts=stack[-1] if stack else (identity,0,{})
        for letter in word[common:]:
            letter_degree,letter_value=frame.L[letter]
            following={}
            for d,image in parts.items():
                if d+letter_degree<=degree:
                    following[d+letter_degree]=frame.mat_mul(image,d,letter_value,letter_degree)
            if value is not None:
                for d,image in qletters.get(letter,{}).items():
                    if count+d<=degree:
                        raw=frame.mat_mul(value,count,image,d)*(-1 if count%2 else 1)
                        following[count+d]=following.get(count+d,0)+raw
            value=frame.mat_mul(value,count,letter_value,letter_degree) if value is not None and count+letter_degree<=degree else None
            count+=letter_degree
            parts=following
            stack.append((value,count,parts))
        if value is not None:
            ordinary[word]=np.trace(value,axis1=-2,axis2=-1)
        if derivative:
            images[word]={d:np.trace(raw,axis1=-2,axis2=-1) for d,raw in parts.items()}
        previous=word
    return ordinary,images


def word_values(words,frame,ordinary,images=None):
    result=np.zeros((len(words),frame.ncols),dtype=np.int64)
    previous,stack=(),[]
    derivative=images is not None
    for index in sorted(range(len(words)),key=lambda i:words[i][1]):
        paths,multi=words[index]
        assert not paths
        length=sum(map(len,multi))
        assert factorial(frame.degree)*frame.N**length*(6*length*frame.N if derivative else 1)<2**62
        common=0
        while common<min(len(multi),len(previous)) and multi[common]==previous[common]:
            common+=1
        del stack[common:]
        value,count,parts=stack[-1] if stack else (np.ones((frame.P,1),dtype=np.int64),0,{})
        for trace in multi[common:]:
            trace_degree=fermions(trace)
            following={}
            if trace in ordinary:
                for d,image in parts.items():
                    if d+trace_degree<=frame.degree:
                        following[d+trace_degree]=frame.scal_mul(image,d,ordinary[trace],trace_degree)
            if value is not None and derivative:
                for d,image in images[trace].items():
                    if count+d<=frame.degree:
                        raw=frame.scal_mul(value,count,image,d)*(-1 if count%2 else 1)
                        following[count+d]=following.get(count+d,0)+raw
            value=frame.scal_mul(value,count,ordinary[trace],trace_degree) if value is not None and count+trace_degree<=frame.degree else None
            count+=trace_degree
            parts=following
            stack.append((value,count,parts))
        if derivative and frame.degree in parts:
            result[index]=parts[frame.degree].reshape(-1)
        elif not derivative and count==frame.degree:
            result[index]=value.reshape(-1)
        previous=multi
    return result


def evaluate(expressions,frame,ordinary=None,images=None,derivative=False,progress=None):
    def eligible(word):
        count=sum(fermions(t) for t in word[1])
        return abs(count-frame.degree)==1 if derivative else count==frame.degree
    expressions=[{w:c for w,c in expr.items() if eligible(w)} for expr in expressions]
    words=sorted({w for expr in expressions for w in expr})
    if not words:
        return np.zeros((len(expressions),frame.ncols),dtype=np.int64)
    bounds={w:factorial(frame.degree)*frame.N**sum(map(len,w[1]))*
        (6*sum(map(len,w[1]))*frame.N if derivative else 1) for w in words}
    maximum=max(sum(abs(int(c))*bounds[w] for w,c in expr.items()) for expr in expressions)
    assert maximum<2**62,('expression overflow bound',maximum)
    index={w:i for i,w in enumerate(words)}
    rr,cc,vv=[],[],[]
    for i,expr in enumerate(expressions):
        for w,c in expr.items():
            assert int(c)==c
            rr.append(i);cc.append(index[w]);vv.append(int(c))
    coefficients=csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(expressions),len(words)))
    if ordinary is None:
        ordinary,images=trace_values({t for _,multi in words for t in multi},frame,derivative)
    result=np.zeros((len(expressions),frame.ncols),dtype=np.int64)
    for start in range(0,len(words),2048):
        end=min(start+2048,len(words))
        result+=coefficients[:,start:end]@word_values(words[start:end],frame,ordinary,images if derivative else None)
        if progress is not None and start%(16*2048)==0:
            progress(end,len(words),maximum)
    return result
