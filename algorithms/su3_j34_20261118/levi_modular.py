"""Derivative-odd sparse elimination and p-adic reconstruction on the scalar restrictions."""
import argparse
import gzip
import heapq
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd, isqrt, lcm
from pathlib import Path

import numpy as np
from flint import ctx, nmod_mat
from scipy.sparse import csr_matrix

import levi_first as levi

P = 1000003
TRANSLATE = {i: (i % 8) + 16*((i//8)%2) + 8*(i//16) for i in range(32)}


def exchanged(word):
    sign = 1
    pairs = []
    for left, right in word:
        a,b = TRANSLATE[left], TRANSLATE[right]
        if a > b:
            sign *= (-1)**(levi.block.odd(a)*levi.block.odd(b))
            a,b = b,a
        pairs.append((a,b))
    v, pairs = levi.block.canonical_multi(tuple(pairs))
    assert v
    return v*sign, pairs


def read(blocks):
    words = levi.domain(0,0)
    lookup = {w:i for i,w in enumerate(words)}
    fold = {}
    pairs = []
    for i,w in enumerate(words):
        sign, partner = exchanged(w)
        j = lookup[partner]
        assert exchanged(partner) == (sign,w)
        if i < j or (i == j and sign < 0):
            column = len(pairs)
            terms = [(i,1)] if i==j else [(i,1),(j,-sign)]
            pairs.append(terms)
            for at, c in terms:
                fold[at] = column,c
    rows, values, labels = [], [], []
    for path in blocks:
        doc = json.loads(gzip.decompress(path.read_bytes()))
        for t, diagonal, key, row, b in doc['equations']:
            if tuple(t)>tuple(sorted(TRANSLATE[k] for k in t)):
                continue
            v = Counter()
            for at,c in row:
                if at in fold:
                    i,sign = fold[at]
                    v[i] += sign*c
            v={i:c for i,c in v.items() if c}
            if v or b:
                rows.append(v)
                values.append(b)
                labels.append([doc['scalar'],doc['power'],t,diagonal,key])
    return words,pairs,rows,values,labels


def sparse(rows, n):
    data, indices, ptr = [], [], [0]
    for row in rows:
        for j,c in row.items():
            indices.append(j)
            data.append(c)
        ptr.append(len(data))
    return csr_matrix((np.array(data,dtype=np.int64),indices,ptr),shape=(len(rows),n))


def pivots(mat):
    echelon, rank = mat.rref()
    result=[]
    j=0
    for i in range(rank):
        while not echelon[i,j]:
            j+=1
        result.append(j); j+=1
    return result


class Factorization:
    def __init__(self, integer, columns):
        rows={i:{j:c%P for j,c in row.items() if c%P} for i,row in enumerate(integer)}
        support=defaultdict(set)
        for i,row in rows.items():
            for j in row:support[j].add(i)
        self.m=len(rows);self.n=columns
        self.ops=[];self.tri=[]
        queue=[(len(v),j) for j,v in support.items() if v];heapq.heapify(queue)
        while queue:
            degree,k=heapq.heappop(queue)
            if degree!=len(support[k]):continue
            pivot=min(support[k],key=lambda i:(len(rows[i]),i))
            if degree>55 or (degree-1)*(len(rows[pivot])-1)>600:break
            row=rows.pop(pivot)
            inverse=pow(row[k],-1,P)
            self.tri.append((k,pivot,inverse,[(j,v) for j,v in row.items() if j!=k]))
            touched=set(row)
            for j in row:support[j].remove(pivot)
            for i in list(support[k]):
                current=rows[i];u=current[k]*inverse%P
                self.ops.append((i,pivot,u))
                for j,v in row.items():
                    old=current.get(j,0);new=(old-u*v)%P
                    if new:
                        current[j]=new
                        if not old:support[j].add(i)
                    elif old:
                        del current[j];support[j].remove(i)
            for j in touched:
                if support[j]:heapq.heappush(queue,(len(support[j]),j))
        ids=[i for i,r in rows.items() if r]
        col=sorted(j for j,rr in support.items() if rr)
        self.cr=[];self.cc=[];self.inv=None;self.invT=None
        print('modular sparse pivots',len(self.tri),'core',len(ids),len(col),flush=True)
        if col:
            idx={j:i for i,j in enumerate(col)}
            mat=np.zeros((len(ids),len(col)),dtype=np.int64)
            for i,r in enumerate(ids):
                for j,v in rows[r].items():mat[i,idx[j]]=v
            cs=pivots(nmod_mat(mat.tolist(),P))
            rs=pivots(nmod_mat(mat[:,cs].T.tolist(),P))
            assert len(rs)==len(cs)
            self.cr=[ids[i] for i in rs];self.cc=[col[i] for i in cs]
            self.inv=nmod_mat(mat[np.ix_(rs,cs)].tolist(),P).inv()
            self.invT=self.inv.transpose()
        self.rank=len(self.tri)+len(self.cr)
        print('rank mod',P,self.rank,flush=True)

    def apply(self, b, transpose):
        if transpose:
            x=[int(v)%P for v in b];v=[0]*self.m
            for k,p,inv,row in self.tri:
                u=x[k]*inv%P;v[p]=u
                for j,a in row:x[j]=(x[j]-a*u)%P
            if self.cc:
                a=self.invT*nmod_mat(len(self.cc),1,[x[j] for j in self.cc],P)
                for i,u in zip(self.cr,a.entries()):v[i]=int(u)
            for i,p,c in reversed(self.ops):v[p]=(v[p]-c*v[i])%P
            return v
        v=[int(a)%P for a in b]
        for i,p,c in self.ops:v[i]=(v[i]-c*v[p])%P
        x=[0]*self.n
        if self.cr:
            a=self.inv*nmod_mat(len(self.cr),1,[v[i] for i in self.cr],P)
            for j,u in zip(self.cc,a.entries()):x[j]=int(u)
        for k,p,inv,row in reversed(self.tri):
            x[k]=(v[p]-sum(a*x[j] for j,a in row))*inv%P
        return x


def reconstruct(x, modulus):
    cap=isqrt(modulus//2)
    a,b=modulus,x%modulus
    c,d=0,1
    while b>cap:
        q=a//b
        a,b=b,a-q*b
        c,d=d,c-q*d
    if d and abs(d)<=cap and gcd(b,d)==1:
        return Fraction(b,d)
    return None


def verify_product(matrix,x,b,scale=1):
    for i,v in enumerate(b):
        a,h=matrix.indptr[i:i+2]
        if sum(int(c)*x[int(j)] for j,c in zip(matrix.indices[a:h],matrix.data[a:h]))!=v*scale:return False
    return True


def dixon(matrix,b,fac,transpose,digits):
    mass=max((sum(abs(int(z)) for z in matrix.data[a:h]) for a,h in zip(matrix.indptr,matrix.indptr[1:])),default=0)
    assert mass*(P-1)<2**63
    x=[0]*matrix.shape[1];cur=list(b);modulus=1
    for step in range(1,digits+1):
        u=fac.apply(cur,transpose)
        prod=(matrix@np.array(u,dtype=np.int64)).tolist()
        diff=[a-v for a,v in zip(cur,prod)]
        failed=next((i for i,v in enumerate(diff) if v%P),None)
        if failed is not None:
            return None,failed,step
        cur=[v//P for v in diff]
        x=[z+modulus*d for z,d in zip(x,u)];modulus*=P
        if step not in (1,2,3,4,6,8) and step%8 and step!=digits:continue
        den=1;cap=isqrt(modulus//2)
        while den<=cap:
            nums=[z*den%modulus for z in x]
            nums=[z if z<=modulus//2 else z-modulus for z in nums]
            j=max(range(len(nums)),key=lambda i:abs(nums[i]))
            if abs(nums[j])<=cap and verify_product(matrix,nums,b,den):
                return {i:Fraction(z,den) for i,z in enumerate(nums) if z},None,step
            v=reconstruct(x[j],modulus)
            if v is None or den%v.denominator==0:break
            den=lcm(den,v.denominator)
    return None,None,digits


def solve(rows,values,n,digits):
    matrix=sparse(rows,n);fac=Factorization(rows,n)
    primal,failed,step=dixon(matrix,values,fac,False,digits)
    if primal is not None:return 'solution',primal,fac.rank,step
    if failed is None:return 'undetermined',{},fac.rank,step
    target=matrix[failed].toarray().ravel().tolist()
    z,_,d=dixon(matrix.T.tocsr(),target,fac,True,digits)
    if z is None:return 'undetermined',{},fac.rank,max(d,step)
    dual=Counter({failed:Fraction(1)})
    for i,v in z.items():dual[i]-=v
    dual={i:v for i,v in dual.items() if v}
    assert sum(v*values[i] for i,v in dual.items())
    return 'separator',dual,fac.rank,d


def run(blocks,out,digits,threads):
    ctx.threads=threads
    words,pairs,rows,b,labels=read(blocks)
    print('odd equations',len(rows),'unknowns',len(pairs),'nnz',sum(map(len,rows)),flush=True)
    kind,co,rank,steps=solve(rows,b,len(pairs),digits)
    result=dict(input=levi.block.RHS_SHA,full_pair_words=len(words),odd_pair_coordinates=len(pairs),
                folded_equations=len(rows),prime=P,rank_mod_prime=rank,digits=steps,kind=kind,
                selected_equations=[[levi.json.loads(gzip.decompress(f.read_bytes()))[k] for k in ('scalar','power')] for f in blocks])
    if kind=='solution':
        r={}
        for k,v in co.items():
            for i,sg in pairs[k]:r[i]=v*sg/120
        result['primitive']=[[words[i],str(v)] for i,v in sorted(r.items())]
    elif kind=='separator':
        result.update(pairing=str(sum(v*b[i] for i,v in co.items())/120),
                      functional=[[labels[i],str(v)] for i,v in sorted(co.items())],
                      negative_exchange_projection=True)
    out.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(kind,len(co),steps,flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--blocks',nargs='+',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--digits',type=int,default=64);p.add_argument('--threads',type=int,default=1)
    assert __debug__;run(**vars(p.parse_args()))
