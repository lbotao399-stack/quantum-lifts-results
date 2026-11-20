"""Project the physical F4 first action on every ordinary trace singlet to the fifth gaugino power."""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np
from flint import fmpq_mat,nmod_mat

import mixed_incoming as m
from f4_colour import basis
from f4_first_source import G, multiply, singlet
from f4_formal_cycles import coefficient_matrix,probe
from f4_singlet_characters import singlets
from mixed33_exact import read
from rank_four_s5 import decode_rows
from q1d_exact import pair_terms
from fortuity.full.q1 import _reorder_sign
from fortuity.full.letters import is_ferm


@lru_cache(None)
def colours():
    gen,_=basis()
    a=gen[0]
    j=next(j for j,b in enumerate(gen[1:],1) if not (a@b-b@a).any())
    h=np.stack((a,gen[j]))
    gram=np.einsum("iab,jba->ij",gen,gen)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    den=int(np.lcm.reduce(abs(gram.diagonal())))
    dual=gen*(den//gram.diagonal())[:,None,None]
    cartan=np.einsum("iab,jba->ij",h,h)
    assert int(cartan[0,0])*int(cartan[1,1])-int(cartan[0,1])**2
    return h,gen,dual,den,j


def top(word, pair=None, exchanged=None):
    h,gen,dual,_,_=colours()
    lengths=list(map(len,word))
    labels=[l for tr in word for l in tr]
    if pair is not None:
        for idx,l in zip(pair,exchanged):
            labels[idx]=l
    assert all(l in G for l in labels) and len(labels)==4
    width=1 if pair is None else len(gen)
    total=0
    for coeff in product(range(2),repeat=4):
        exterior=[2*k+G.index(l) for k,l in zip(coeff,labels)]
        if len(set(exterior))!=4:
            continue
        sign=(-1)**sum(a>b for i,a in enumerate(exterior) for b in exterior[i+1:])
        product_val=np.ones(width,np.int64)
        offset=0
        for length in lengths:
            mat=None
            for pos in range(offset,offset+length):
                x=h[coeff[pos]]
                if pair is not None and pos in pair:
                    t=gen if pos==pair[0] else dual
                    x=t@x-x@t
                mat=x if mat is None else mat@x
            product_val*=mat.trace(axis1=-2,axis2=-1)
            offset+=length
        total+=sign*sum(map(int,product_val))
    return total


@lru_cache(None)
def reduced_action(word):
    letters=[l for tr in word for l in tr]
    i,j=[k for k,l in enumerate(letters) if l[2]]
    old=list(map(is_ferm,letters))
    move=(-1)**(old[i]*sum(old[:i])+old[j]*(sum(old[:j])-old[i]))
    ref=[i,j]+[k for k in range(4) if k not in (i,j)]
    result=Fraction(0)
    for c,left,right in pair_terms(letters[i],letters[j]):
        post=old.copy()
        post[i],post[j]=is_ferm(left),is_ferm(right)
        phase=_reorder_sign(list(range(4)),ref,lambda k:post[k])
        result+=c*move*phase*top(word,(i,j),(left,right))
    return result/Fraction(colours()[3]*top((G,G)))


def action(word):
    removable=[]
    for k,tr in enumerate(word):
        if tr==G or tr==G[::-1]:
            removable.append((k,1 if tr==G else -1))
    if len(removable)<3:
        return Fraction(0)
    removed=dict(removable[:3])
    reduced=tuple(tr for i,tr in enumerate(word) if i not in removed)
    if not any(l[2] for tr in reduced for l in tr):
        return Fraction(0)
    return (int(np.prod(list(removed.values()))))*reduced_action(reduced)


def run(plan,previous,work,workers):
    work.mkdir(exist_ok=True,parents=True)
    z,old=read(plan),read(previous)
    p=old["prime"]
    weight=(2,2,0,0,0)
    lower=m.prepare(weight,3)
    image_words=list(map(m.decode,lower["target_dictionary"]))
    evaluated,_=probe(lower["images"],image_words,singlets("F4",weight,4),p,m.library(work,27,p),workers,750022)
    minor=m.minor(evaluated,p)
    assert singlets("F4",weight,4)[4]-minor["rank"]==1 and top((G,G))
    words=list(map(m.decode,z["source_dictionary"]))
    coeff=[action(w) if m.nferm(w)==9 else Fraction(0) for w in words]
    quantum=[sum((c*coeff[j] for j,c in row),Fraction(0)) for row in z["sources"]]
    assert sum(c*action(w) for w,c in multiply(singlet(),4).items())==18
    cols=old["formal_image_word_columns"]
    lookup={j:i for i,j in enumerate(cols)}
    n=len(cols)
    formal=[[(lookup[j],c) for j,c in row if j in lookup] for row in z["images"]]
    modular=coefficient_matrix(formal,n,p).toarray()
    rows=m.pivotcols(nmod_mat(modular.T.tolist(),p))
    assert len(rows)==n==old["formal_terminal_singlets"]
    exact=fmpq_mat(n,n)
    for i,r in enumerate(rows):
        for j,c in formal[r]:
            exact[i,j]=c
    sol=exact.solve(fmpq_mat([[str(quantum[r])] for r in rows]))
    values=[Fraction(str(v)) for v in sol.entries()]
    assert [sum((c*values[j] for j,c in row),Fraction(0)) for row in formal]==quantum
    num_nonzero=sum(bool(v) for v in quantum)
    print("four-letter rank",minor["rank"],"pure S2",top((G,G)),"first-action rows",num_nonzero,
          "formal identities",len(quantum),flush=True)
    record=dict(prime_for_pivot_selection=p,integer_cartan_indices=[0,colours()[-1]],trace_dual_denominator=colours()[3],
                S2_top_coefficient=top((G,G)),four_letter_character=singlets("F4",weight,4)[4],
                four_letter_boundary=minor,source_words=len(words),nine_fermion_words=sum(m.nferm(w)==9 for w in words),
                formal_singlets=len(quantum),nonzero_rows=[[i,str(c)] for i,c in enumerate(quantum) if c],
                formal_rank=n,formal_columns=cols,independent_source_rows=rows,
                covector=[[cols[i],str(c)] for i,c in enumerate(values) if c],rational_residual_nonzero=0)
    m.write(work/"first_functional.json",record)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("plan","previous","work"):
        p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--workers",type=int,default=6)
    assert __debug__
    run(**vars(p.parse_args()))
