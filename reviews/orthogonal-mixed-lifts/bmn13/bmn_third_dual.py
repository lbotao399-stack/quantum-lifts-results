"""Pull the corrected SO(13) second dual back to the derivative-free first residual."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product,permutations
from math import factorial,prod,lcm
from pathlib import Path
import json,sys

def ordered_dual(data,rows):
    from mixed_relations import check
    from orthogonal_glueball_first_pair import orientations, GHOSTS,PHI
    from rank_four_lowered_replay import terminal
    from mixed_incoming import so_canon
    from orthogonal_glueball_third_pair import encode
    from terminal_commutators import words
    from mixed_ward import trace_to_planes
    from flint import fmpq_mat
    dual=Counter()
    for w,row in data['second']:
        k=sum((6 if i==0 else -1 if i==1 else 0)*Fraction(v) for i,v in row)
        if k:dual[tuple(map(tuple,w))]=k
    coords=fmpq_mat([r['planes'] for r in rows['rows']])
    values=fmpq_mat([[str(6*Fraction(r['scores'][0])-Fraction(r['scores'][1]))] for r in rows['rows']])
    def piv(m):
        red,d=m.rref();return [next(j for j in range(red.ncols()) if red[i,j]) for i in range(d)]
    ri,ci=piv(coords.transpose()),piv(coords)
    sol=fmpq_mat([[coords[i,j] for j in ci] for i in ri]).solve(fmpq_mat([[values[i,0]] for i in ri]))
    assert all(sum(coords[i,j]*sol[k,0] for k,j in enumerate(ci))==values[i,0] for i in range(coords.nrows()))
    basis=words((4,4,4,2,2))
    # At most one merge of two sites is possible for the five selected columns.
    selected=[basis[rows['source_plane_indices'][j]] for j in ci]
    from mixed_ward import normalize
    trvar={}
    for word in basis:
        if len(word)<7:continue
        terms=[([ (k,) for k in range(len(word)) ],word,1)]
        if len(word)==8:
            for i in range(8):
                for j in range(i+1,8):
                    g=tuple(x+y for x,y in zip(word[i],word[j]))
                    if max(g[:2])>1:continue
                    partition=[(k,) for k in range(8) if k not in (i,j)]+[(i,j)]
                    blocks=[tuple(sum(word[k][n] for k in B) for n in range(5)) for B in partition]
                    seq=sorted(range(len(blocks)),key=lambda h:(sum(blocks[h]),blocks[h]))
                    slot={p:x for x,h in enumerate(seq) for p in partition[h]}
                    odd=[2*slot[p]+n for p,tr in enumerate(word) for n in (0,1) if tr[n]]
                    sign=(-1)**sum(x>y for s,x in enumerate(odd) for y in odd[s+1:])
                    terms.append((partition,tuple(blocks[h] for h in seq),sign))
        val=Fraction()
        for _,v,s in terms:
            if v in selected:
                val+=s*2**len(word)*prod(factorial(n) for n in Counter(v).values())*Fraction(str(sol[selected.index(v),0]))
        if val:trvar[word]=val
    seen=set()
    for monomial,k in trvar.items():
        opts=[orientations(tuple(l for l,n in zip(GHOSTS+PHI,t) for _ in range(n))) for t in monomial]
        for ts in product(*opts):
            s,w=so_canon(ts)
            if not s or w in seen:continue
            seen.add(w);sign,mon=terminal(w);assert sign and mon==monomial
            dual[encode(w)]-=k*sign
    return {w:n for w,n in dual.items() if n},[str(v) for v in sol.entries()]

def run(a):
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    data=json.loads(a.dual.read_text());rows=json.loads(a.rows.read_text());dual,c=ordered_dual(data,rows)
    den=lcm(*(v.denominator for v in dual.values()))
    with a.out.open('w') as out:
        for w,n in sorted(dual.items()):
            print(int(n*den),'/'.join(''.join(str(x) for x in t) for t in w),file=out)
    print('corrected second covector',len(dual),'denominator',den,'plane corrections',c,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('repo','dual','rows','out'):p.add_argument('--'+k,type=Path,required=True)
    run(p.parse_args())
