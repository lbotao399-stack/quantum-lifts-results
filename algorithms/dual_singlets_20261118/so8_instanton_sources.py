"""Protected products in the five possible incoming SO(8) level-18 singlets."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from math import gcd, lcm, prod
from pathlib import Path
import random

from flint import fmpz_mat
import mixed_incoming as m
from mixed11_products import factor, partitions, sign
from fortuity.full.letters import is_ferm, Q_letter
from full_extremal import action as raise_action


def cleaned(row):
    return {w: n for w, n in row.items() if n}


def outside(expr, axis, odd):
    out = Counter()
    for (word,), coefficient in expr.items():
        phase = 1
        for i, (a, b, theta) in enumerate(word):
            if not odd or axis not in theta:
                if odd:
                    c = phase * (-1)**sum(j > axis for j in theta)
                    letter = (a, b, tuple(sorted((*theta, axis))))
                else:
                    c = (a, b)[axis] + 1
                    letter = (a + int(axis == 0), b + int(axis == 1), theta)
                orientation, target = m.so_canon((word[:i] + (letter,) + word[i+1:],))
                if orientation:
                    out[target] += orientation * c * coefficient
            if odd and is_ferm((a,b,theta)):
                phase = -phase
    return cleaned(out)


@lru_cache(None)
def component(spec, flavours, spins):
    length, inner_spin, outer_theta, outer_spin = spec
    n = length - inner_spin
    assert len(flavours) == n + outer_theta and len(spins) == inner_spin + outer_spin
    row = factor((n, inner_spin, 0), flavours[:n], spins[:inner_spin])
    for i in flavours[n:]:
        row = outside(row, i, True)
    for i in spins[inner_spin:]:
        row = outside(row, i, False)
    assert not m.q(row)
    return row


def families(outer=2, spins=2):
    choices = [(size,h,d,z) for size in (2,4,6) for h in range(3)
               for d in range(outer+1) for z in range(outer+1-d) if h <= size]
    for count in (3,2,1):
        for family in combinations_with_replacement(choices,count):
            if (sum(x[0] for x in family)==6 and sum(x[2]+x[3] for x in family)==outer
                    and sum(x[1]+x[3] for x in family)==spins):
                yield family


def contracted(family, spin_values, flavour_values):
    out = Counter()
    for spin,spin_sign in spin_values:
        for flavours,flavour_sign in flavour_values:
            i = j = 0
            pieces = []
            for spec in family:
                size,h,d,z = spec
                fi, sj = i+size-h+d, j+h+z
                pieces.append(component(spec,tuple(flavours[i:fi]),tuple(spin[j:sj])))
                i,j = fi,sj
            assert i == len(flavours) and j == len(spin)
            for terms in product(*(e.items() for e in pieces)):
                words = tuple(w[0] for w,c in terms)
                s, w = m.so_canon(words)
                if s:
                    out[w] += spin_sign*flavour_sign*s*prod(c for _,c in terms)
    return cleaned(out)


def singlet(family, triples):
    assignments=[]
    for left,right in product(permutations(range(3)),repeat=2):
        flavours=[None]*6
        for positions,indices in zip(triples,(left,right)):
            for i,j in zip(positions,indices):
                flavours[i]=j
        assignments.append((flavours,sign(left)*sign(right)))
    return contracted(family,[((0,1),1),((1,0),-1)],assignments)


def fundamental_labels(count, free_slot):
    records=[]
    for perm in permutations(range(count)):
        labels=list(perm)
        labels.insert(free_slot,0)
        records.append((labels,sign(perm)))
    return records


def cartan(expr):
    result = Counter()
    for w,c in expr.items():
        parts = []
        for tr in w:
            odd = [v for v in tr if is_ferm(v)]
            if len(tr)%2 or len(odd)!=len(set(odd)):
                break
            c *= sign(odd)
            parts.append(tuple(sorted(tr)))
        else:
            s,key = m.fw.canon_multi(parts)
            if s:
                result[key] += c*s
    return cleaned(result)


def primitive(expr):
    denominator = lcm(*(Fraction(c).denominator for c in expr.values()))
    scaled = {w:int(c*denominator) for w,c in expr.items()}
    common = gcd(*scaled.values())
    return {w:n//common for w,n in scaled.items()}, Fraction(denominator,common)


def matrix(rows, columns):
    return fmpz_mat([[r.get(c,0) for c in columns] for r in rows])


def cartan_probe(row, assignments):
    # On one colour plane H=[[0,1],[-1,0]]; Tr H^m=2*(-1)^(m/2).
    # Each odd letter is a linear combination of four exterior generators.
    def multiply(a,b):
        value=Counter()
        for i,x in a.items():
            for j,y in b.items():
                if not i&j:
                    degree_swaps=sum((j&((1<<k)-1)).bit_count() for k in range(4) if (i>>k)&1)
                    value[i|j]+=(-1)**degree_swaps*x*y
        return cleaned(value)
    @lru_cache(None)
    def trace(word):
        out=Counter()
        for values in assignments:
            term={0:2*(-1)**(len(word)//2)}
            for l in word:
                term=multiply(term,values[l])
            out.update(term)
        return cleaned(out)
    result=Counter()
    for multi,c in row.items():
        term={0:c}
        for tr in multi:
            term=multiply(term,trace(tr))
        result.update(term)
    return cleaned(result)


def run():
    return certify_candidates(singlet_candidates(), 5, (1,1,2,2,2))


def singlet_candidates():
    for family in families():
        for triples in partitions(tuple(range(6))):
            yield dict(factors=family,flavour_triples=triples),singlet(family,triples)


def lower_candidates():
    for family in families(outer=1,spins=3):
        for spin,flavour in product(range(3),range(4)):
            expr=contracted(family,fundamental_labels(2,spin),fundamental_labels(3,flavour))
            yield dict(factors=family,free_spin_slot=spin,free_flavour_slot=flavour),expr


def middle_candidates():
    flavour=[(p,sign(p)) for p in permutations(range(3))]
    for family in families(outer=1,spins=4):
        for partner in (1,2,3):
            other=[i for i in range(1,4) if i!=partner]
            spins=[]
            for left,right in product(((0,1),(1,0)),repeat=2):
                labels=[0]*4
                labels[0],labels[partner]=left
                labels[other[0]],labels[other[1]]=right
                spins.append((labels,sign(left)*sign(right)))
            yield dict(factors=family,first_spin_partner=partner),contracted(family,spins,flavour)


def certify_candidates(candidates, dimension, weight):
    selected, rows, columns = [], [], []
    for recipe,expr in candidates:
        pic = cartan(expr)
        if not pic:
            continue
        pic, scale = primitive(pic)
        trial = rows+[pic]
        support = sorted(set().union(*(r.keys() for r in trial)))
        mat = matrix(trial,support)
        if mat.rank() == len(rows):
            continue
        assert mat.rank()==len(trial) and not m.q(expr)
        for root in ((0,1),(1,2),("spin",1)):
            assert not raise_action(expr,*root)
        for traces in expr:
            letters=[v for tr in traces for v in tr]
            assert len(letters)==6 and tuple(sum(x[k] for x in letters) for k in (0,1)) + tuple(
                sum(i in x[2] for x in letters) for i in range(3))==weight
        rows.append(pic)
        selected.append(dict(**recipe,scale=str(scale),full_trace_terms=len(expr),cartan_terms=len(pic)))
        if len(rows)==dimension:
            break
    assert len(rows)==dimension
    support = sorted(set().union(*(r.keys() for r in rows)))
    allcols = matrix(rows,support)
    for j in range(len(support)):
        if matrix(rows,columns+[support[j]]).rank()>len(columns):
            columns.append(support[j])
    coeff = matrix(rows,columns)
    assert len(rows)==len(columns) and coeff.det()
    generators = sorted({l for row in rows for multi in row for tr in multi for l in tr})
    rng = random.Random(831799)
    probes = []
    for g in (3,4):
        measured,labels=[],[]
        for i in range(12):
            assignments=[{l:({1<<j:rng.randrange(-2,3) for j in range(4)} if is_ferm(l)
                              else {0:rng.randrange(-2,3)}) for l in generators} for _ in range(g)]
            vals=[cartan_probe(r,assignments) for r in rows]
            for mask in range(16):
                col=[v.get(mask,0) for v in vals]
                if fmpz_mat(measured+[col]).rank()>len(measured):
                    measured.append(col)
                    labels.append([i,mask])
            if len(measured)==len(rows):
                break
        mat=fmpz_mat(measured)
        assert mat.nrows()==mat.ncols() and mat.det()
        probes.append(dict(rank=g,point_and_exterior_masks=labels,
                           matrix=[list(map(int,r)) for r in mat.tolist()],determinant=str(mat.det())))
    encoded = [[m.encode(c),[r.get(c,0) for r in rows]] for c in columns]
    return dict(weight=weight,length=6,doubled_charges=m.sector(weight,6),
                products=selected,power_sum_pivots=encoded,
                power_sum_determinant=str(coeff.det()),cartan_probes=probes,
                probe_seed=831799,protected_source_span=len(rows))

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    assert __debug__
    m.write(args.out,dict(level18=run(),level17=certify_candidates(lower_candidates(),9,(2,1,2,1,1)),
                          level18_middle=certify_candidates(middle_candidates(),2,(2,2,1,1,1))))
