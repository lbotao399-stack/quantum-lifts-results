"""Reconstruct the two rational first-page detectors and the return at level 50."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
import json

import descendants as d
import first_slice as f
from su2_bmn_character50 import characters


def add_variation(word,operator,i):
    out=defaultdict(int)
    for j,l in enumerate(word):
        value=change(l,operator,i)
        if value is None:
            continue
        n,new=value
        exponent=sum(f.parity(s) for _,_,s in word[:j])
        out[word[:j]+(new,)+word[j+1:]]+=n*(-1)**exponent
    return out


def change(letter,operator,i):
    p,q,s=letter
    sg=(-1)**sum(k>i for k in f.components(s))
    if operator=='D' and not s&(1<<i):
        return sg,(p,q,s|(1<<i))
    degree=p+q+s.bit_count()
    if operator=='S' and s&(1<<i) and degree>1:
        return -(degree-1)*sg,(p,q,s^(1<<i))


def native_tree(letter):
    p,q,s=letter
    out=defaultdict(int)
    for a,b,t in product(range(p+1),range(q+1),f.flavours):
        if t&~s:
            continue
        u=s^t
        one,two=(a,b,t),(p-a,q-b,u)
        if one==(0,0,0) or two==(0,0,0):
            continue
        n=(-1)**(s.bit_count()+f.parity(t)*u.bit_count())*f.sign(f.components(t)+f.components(u))
        out[one,two]+=n
    return d.clean(out)


def native_pair(word):
    x,y=word
    assert x[:2]==y[:2]==(0,0)
    out=defaultdict(int)
    for n,s,t,plus in f.pair_rule(x[2],y[2]):
        a,b=((1,0),(0,1)) if plus else ((0,1),(1,0))
        out[(*a,s),(*b,t)]+=n
    return d.clean(out)


def verify_ward():
    tree_checks=pair_checks=0
    for p,q,s,op,i in product(range(2),range(2),f.flavours,('D','S'),d.axes):
        l=(p,q,s)
        if l==(0,0,0):
            continue
        combined=defaultdict(int)
        changed=change(l,op,i)
        if changed:
            n,new=changed
            d.add_to(combined,native_tree(new),n)
        for word,n in native_tree(l).items():
            d.add_to(combined,add_variation(word,op,i),n)
        inner=defaultdict(int)
        if op=='D':
            scalar=(0,0,1<<i)
            inner[l,scalar]+=1
            inner[scalar,l]-=1
        assert d.clean(combined)==d.clean(inner)
        tree_checks+=1
    for s,t,op,i in product(range(1,8),range(1,8),('D','S'),d.axes):
        w=((0,0,s),(0,0,t))
        combined=defaultdict(int)
        for word,n in native_pair(w).items():
            d.add_to(combined,add_variation(word,op,i),n)
        for word,n in add_variation(w,op,i).items():
            d.add_to(combined,native_pair(word),n)
        assert not d.clean(combined)
        pair_checks+=1
    return {'tree':tree_checks,'ordered_first_pairs':pair_checks}


def fraction_remainder(basis,vec):
    vals={k:Fraction(n) for k,n in vec.items()}
    for col,pivot in sorted(basis.rows.items()):
        c=vals.get(col,0)
        if c:
            tmp=defaultdict(int,vals)
            d.add_to(tmp,pivot,-c/pivot[col])
            vals=d.clean(tmp)
    return vals


def validate():
    check=verify_ward()
    A=d.core(2)
    B=d.dressed(1)
    for i in d.axes:
        assert not d.act('Q',0,d.hair(i))
        assert all(not d.act('D',k,d.hair(i)) for k in d.axes)
    assert all(not d.act('Q',0,z) for z in [d.core(0),d.core(1),A,B])
    weights=[(8,8,8),(9,8,8),(9,9,7),(8,8,7)]
    dims=characters(weights)
    ranks={}
    for w in weights:
        data=d.frame.run(w,0)
        assert all(v['dimension']==dims[','.join(map(str,w))].get(str(y),0) for y,v in data.items())
        assert all(int(k) in data for k in dims[','.join(map(str,w))])
        ranks[','.join(map(str,w))]=[
            [y,v['dimension'],v['incoming'],v['out'],v['cohomology_upper']]
            for y,v in sorted(data.items())]
    assert [[r[0],r[4]] for r in ranks['8,8,8'] if r[4]]==[[11,2]]
    assert [[r[0],r[4]] for r in ranks['9,8,8'] if r[4]]==[[11,2]]
    assert not any(r[4] for r in ranks['9,9,7'])
    source=d.boundaries((8,8,8),10)
    b=fraction_remainder(source,d.slice_vector(B))
    assert b
    lower=d.boundaries((8,8,7),10)
    x=fraction_remainder(lower,d.covariant('S',A))
    y=fraction_remainder(lower,d.covariant('S',B))
    assert y and {k:102*c for k,c in x.items()}==y
    coeffs=(Fraction(-1,3),Fraction(-34))
    for z,h in zip((A,B),coeffs):
        t=defaultdict(int)
        u=defaultdict(int)
        for i in d.axes:
            d.add_to(t,d.act('D',i,d.act('S',i,z)))
            d.add_to(u,d.act('S',i,d.act('D',i,z)))
        assert fraction_remainder(source,d.slice_vector(t))=={k:h*c for k,c in b.items()}
        d.add_to(u,t)
        d.add_to(u,z,48)
        assert not d.clean(u)

    target=sorted(f.basis((7,7,7),11,None))
    prev=list(f.basis((7,7,7),10,None))
    assert (len(prev),len(target))==(21,12)
    # This covector is independently checked on every predecessor below.
    weights=(1,Fraction(1,2),1)
    dual={entry:weights[i%3]*(1 if i//3 in (0,2,3) else -1)
          for i,entry in enumerate(target)}
    def score(z):
        assert set(z)<=set(target)
        return sum(dual[k]*v for k,v in z.items())
    trees=d.frame.Rank(0)
    for key in prev:
        z=f.tree(key,None)
        assert score(z)==0
        trees.insert(z)
    assert len(trees)==11
    pair=[score(f.quantum(z,None)) for z in (A,B)]
    ward=[c*pair[1] for c in coeffs]
    determinant=pair[0]*ward[1]-pair[1]*ward[0]
    assert pair==[-9600,-172800] and ward==[57600,5875200] and determinant==-46448640000
    # Tabulate the two independent classes and their D descendants modulo all
    # rational weight-zero predecessors on the scalar frame.
    lifted=d.joined((9,8,8),{'D_O2':d.covariant('D',A),'D_dressing':d.covariant('D',B)})
    assert lifted['additional_rank']==2
    return dict(ward_checks=check,cochains_incoming_outgoing_cohomology=ranks,
                classical_lowering_ratio='1/102',classical_DS_coefficients=list(map(str,coeffs)),
                boundary_matrix_shape=[len(target),len(prev)],boundary_rank=len(trees),
                dual=list(map(str,dual.values())),first_pairing=[list(map(int,pair)),list(map(int,ward))],
                determinant=int(determinant),SD_return_eigenvalues=[-48,-14],
                first_rank_level48=2,first_rank_level50=2)


if __name__=='__main__':
    if not __debug__:
        raise ValueError('Enable assertions')
    print(json.dumps(validate(),indent=2))
