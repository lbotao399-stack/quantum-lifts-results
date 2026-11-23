"""Rational ranks of consecutive terminal supercharge maps for orthogonal and symplectic trace cochains."""
import argparse
import hashlib
import math
import sys
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat
import mixed_incoming as m
from charge_maps import _multiplicity, content_patterns, multiplicity_parts
from rank_four_s5 import Q1d_multi, decode_rows
from rank_four_sixth import IntegralFrame


def dimension(colours, weight, length, group="SO"):
    value=0
    for spin,ss in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            w=(weight[0]+spin,weight[1]-spin)+tuple(weight[j+2]+perm[j]-j for j in range(3))
            sign=ss*(-1)**sum(x>y for k,x in enumerate(perm) for y in perm[k+1:])
            for bosons,fermions,count in content_patterns(w,length):
                mult=(_multiplicity(bosons,fermions,0,colours) if group=="SO"
                      else multiplicity_parts(bosons,fermions,colours//2))
                value+=sign*count*mult
    assert value>=0
    return value


@lru_cache(100000)
def first_word(word,colours,group):
    double=[l for tr in word for l in tr if sum(l[:2])+len(l[2])>1]
    assert len(double)==1 and sum(double[0][:2])+len(double[0][2])==2
    if len(double[0][2])<2:
        return ()
    out=Counter()
    if group=="Sp":
        from q1d_sp_exact import Q1d_multi as quantum
    else:
        quantum=Q1d_multi
    for term,n in quantum(word,colours,"SO_odd" if group=="SO" else "Sp").items():
        sg,w=m.so_canon(term)
        out[w]+=sg*n
    return tuple((w,n) for w,n in out.items() if n)


def first(row,colours,group):
    out=Counter()
    for w,c in row.items():
        for v,n in first_word(w,colours,group):
            out[v]+=c*n
    return {v:n for v,n in out.items() if n}


class IntegralTerminalFrame(IntegralFrame):
    def mat_mul_letter(self, arr, degree, letter):
        odd,b=self.L[letter]
        if odd:
            return super().mat_mul_letter(arr,degree,letter)
        return arr@b,degree


def sample(words,degree,colours,points,seed,group,alphabet=None):
    used={l for w in words for tr in w for l in tr}
    letters=sorted(used if alphabet is None else alphabet)
    assert used.issubset(letters)
    k=degree+1
    f=IntegralTerminalFrame(colours,points,k,[degree],letters,seed=1)
    rng=np.random.default_rng(seed)
    if group=="Sp":
        r=colours//2
        eye=np.eye(r,dtype=np.int64)
        zero=np.zeros_like(eye)
        omega=np.block([[zero,eye],[-eye,zero]])
    for l in letters:
        odd,template=f.L[l]
        vals=np.zeros_like(template)
        for a in vals.reshape(-1,colours,colours):
            for _ in range(2):
                perm=rng.permutation(colours)
                if group.startswith("SU"):
                    if rng.integers(2):
                        signs=rng.integers(-1,2,size=colours)
                        a[perm,np.roll(perm,1)]+=signs
                    else:
                        n=int(rng.integers(-1,2))
                        a[perm[0],perm[0]]+=n
                        a[perm[1],perm[1]]-=n
                    continue
                signs=rng.integers(-1,2,size=colours//2)
                i,j=perm[0:2*(colours//2):2],perm[1:2*(colours//2):2]
                if group=="SO":
                    a[i,j]+=signs
                    a[j,i]-=signs
                else:
                    # Symmetric signed partial permutations include fixed points.
                    swap=rng.integers(2,size=r).astype(bool)
                    a[i[swap],j[swap]]+=signs[swap]
                    a[j[swap],i[swap]]+=signs[swap]
                    a[i[~swap],i[~swap]]+=signs[~swap]
                    a[j[~swap],j[~swap]]+=rng.integers(-1,2,size=np.count_nonzero(~swap))
        if group=="Sp":
            vals=omega@vals
            assert not (vals.swapaxes(-1,-2)@omega+omega@vals).any()
        elif group.startswith("SU"):
            assert not vals.trace(axis1=-2,axis2=-1).any()
        else:
            assert not (vals+vals.swapaxes(-1,-2)).any()
        assert abs(vals).sum(axis=-1).max()<=2
        f.L[l]=odd,vals
    bound=max(math.factorial(degree)*2**sum(map(len,w))*colours**len(w) for w in words)
    assert bound<2**63 and all(m.nferm(w)==degree for w in words)
    return f.eval_words(words,dtype=np.int64)


def row_scales(right):
    return [math.lcm(*(Fraction(n).denominator for n in r.values())) for r in right]


def evaluate(rows,scales,word_values,width):
    matrix=[]
    for row,scale in zip(rows,scales):
        out=[0]*width
        for w,n in row.items():
            x=n*scale
            assert x==int(x)
            for j,v in enumerate(word_values[w]):
                out[j]+=int(x)*v
        matrix.append(out)
    return matrix


def coordinates(terminal,diagonal,advanced,scales,prior_scales,mu,colours,group):
    dim=dimension(colours,mu,sum(mu),group)
    all_words=sorted({w for block in (terminal,diagonal,advanced) for row in block for w in row})
    if not all_words:
        assert dim==0
        return [[] for _ in diagonal],[[] for _ in advanced],dict(dimension=0)
    degree=mu[0]+mu[1]
    seed=164209+317*degree+11*colours
    # The polynomial bounds apply to raw words. All singlet and quantum combinations use Python integers.
    blocks=[];cols=0
    for count in range(48):
        size=min(8,max(2,math.ceil((dim+8-cols)/(degree+1))))
        blocks.append(sample(all_words,degree,colours,size,seed+4099*count,group))
        values=np.column_stack(blocks)
        cols=values.shape[1]
        if cols<dim:
            continue
        table=dict(zip(all_words,values.tolist()))
        reference=evaluate(terminal,[1]*len(terminal),table,cols)
        if not dim:
            assert all(not any(row) for row in reference)
            piv=[]
            break
        a=nmod_mat(reference,1000003)
        piv=m.pivotcols(a)
        assert len(piv)<=dim
        if len(piv)==dim:
            break
    else:
        raise ValueError(("terminal evaluation failed to span",mu,dim,len(piv)))
    chosen=dict(zip(all_words,values[:,piv].tolist()))
    right=evaluate(advanced,prior_scales,chosen,dim)
    left=evaluate(diagonal,scales,chosen,dim)
    return left,right,dict(dimension=dim,formal=len(terminal),points=sum(v.shape[1] for v in blocks)//(degree+1),
                           seeds=[seed+4099*j for j in range(len(blocks))],columns=piv)


def interval(blocks,left,right):
    ranks={}
    multiplicities={}
    for width in range(1,len(blocks)+1):
        for lo in range(len(blocks)-width+1):
            hi=lo+width-1
            n=sum(blocks[t]["dimension"] for t in range(lo,hi+1))
            count=sum(len(left[t]) for t in range(lo,hi+1))
            mat=fmpz_mat(count,n)
            i=j=0
            for t in range(lo,hi+1):
                for row in range(len(left[t])):
                    for k,x in enumerate(left[t][row]):
                        mat[i,k+j]=x
                    if t<hi:
                        for k,x in enumerate(right[t+1][row]):
                            mat[i,k+j+blocks[t]["dimension"]]=x
                    i+=1
                j+=blocks[t]["dimension"]
            assert i==count and j==n
            # FLINT computes rank over Q using the integral matrix.
            r=mat.rank()
            ranks[lo,hi]=r
            mlt=r-ranks.get((lo+1,hi),0)-ranks.get((lo,hi-1),0)+ranks.get((lo+1,hi-1),0)
            assert mlt>=0
            multiplicities[lo,hi]=mlt
    return [[lo,hi,r,multiplicities[lo,hi]] for (lo,hi),r in sorted(ranks.items())]


def run(flavours,colours,work,group="SO"):
    work.mkdir(parents=True,exist_ok=True)
    assert group in ("SO","Sp")
    if group=="Sp":
        assert all(size%2==0 for size in colours)
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"so7_cartan_target_20261118"/"native"))
    weights=[(t,t)+tuple(n-t for n in flavours) for t in range(min(flavours)+1)]
    plans=[];targets=[]
    for t,w in enumerate(weights):
        grouped=m.prepare(w,sum(w)-1)
        sources=decode_rows(grouped,"sources"); images=decode_rows(grouped,"images")
        terminal=decode_rows(m.prepare(w,sum(w)),"sources")
        assert all(not m.q(row) for row in terminal)
        plans.append((sources,images));targets.append(terminal)
        print("formal",t,len(sources),len(terminal),flush=True)
    for size in colours:
        records=[];left=[];right=[]
        incoming=[];previous=[]
        for t,(w,(src,images),target) in enumerate(zip(weights,plans,targets)):
            quantum=[first(row,size,group) for row in src]
            if t+1==len(weights):
                assert not any(quantum)
            scales=row_scales(quantum)
            a,b,record=coordinates(target,images,incoming,scales,previous,w,size,group)
            left.append(a);right.append(b)
            record.update(weight=w,formal_parents=len(src),minimum_row_multiplier=min(scales,default=1),
                          maximum_row_multiplier=max(scales,default=1))
            records.append(record)
            incoming=quantum;previous=scales
        data=interval(records,left,right)
        totals=[sum(x[3] for x in data if x[1]-x[0]==r) for r in range(len(weights))]
        matrix_digest=hashlib.sha256(m.compact([left,right])).hexdigest()
        cert=dict(colours=size,parity="orthogonal traces",flavours=flavours,targets=records,
                  interval_ranks_and_pairings=data,ranks_by_order=totals,integer_matrix_sha256=matrix_digest)
        if group=="Sp":
            cert["parity"]="full symplectic invariants"
        m.write(work/("terminal_"+"_".join(map(str,flavours))+"_"+group.lower()+str(size)+".json"),cert)
        print(size,flavours,totals,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flavours",nargs=3,required=True,type=int)
    p.add_argument("--colours",nargs="+",required=True,type=int)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--group",choices=("SO","Sp"),default="SO")
    a=p.parse_args()
    assert __debug__ and a.flavours==sorted(a.flavours,reverse=True) and min(a.colours)>=3
    run(a.flavours,a.colours,a.work,a.group)
