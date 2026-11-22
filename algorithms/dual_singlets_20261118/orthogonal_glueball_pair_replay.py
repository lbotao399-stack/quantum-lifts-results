"""Polarize disjoint colour cycles in a literal Pfaffian and in the adjoint sum."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np
from orthogonal_glueball_first_pair import ENTRIES, GHOSTS, PHI, PfaffianCoefficient, highest_product
from orthogonal_odd_trace_primitive import primitive
from fortuity.full.letters import is_ferm
from fortuity.full.q1 import pair_rule


def permutation_sign(sequence):
    return (-1)**sum(i>j for k,i in enumerate(sequence) for j in sequence[k+1:])


def edge_orientations(length):
    return [(i,(i+1)%length) for i in range(length)]


def cycle_pin(word, double_count, projector):
    """One independent polarizing edge per occurrence on a connected cycle."""
    n=len(word)
    edges=edge_orientations(n)
    odd_mask=sum(1<<i for i,v in enumerate(word) if is_ferm(v))
    # Key consists of the exterior mask (six markers, then odd probe edges),
    # all occupied edges (even probe edges are also multilinearized), the
    # scalar marker monomial and the number of double paths.
    ground=(0,0,0,0)
    entries=[[Counter() for _ in range(2*n)] for _ in range(2*n)]
    def put(i,j,key,value):
        if i<j:
            entries[i][j][key]+=value
    for i in range(n):
        put(i,n+i,ground,1)
    for k,(u,v) in enumerate(edges):
        field=word[k]
        for i,j,d in ((u,v,1),(v,u,-1)):
            for s,matrix in enumerate(ENTRIES):
                if field in GHOSTS:
                    bit=(1<<(3*GHOSTS.index(field)+s)) | (1<<(6+k))
                    weight=0
                else:
                    bit=0
                    weight=1<<(3*(3*PHI.index(field)+s))
                for a,b,c in matrix:
                    # J times the marker matrix.
                    put((1-a)*n+i,b*n+j,(bit,1<<k,weight,0),d*c*(1 if a else -1))
    for k,l in permutations(range(n),2):
        if (word[k],word[l]) not in (GHOSTS,GHOSTS[::-1]):
            continue
        orient=1 if word[k]==GHOSTS[0] else -1
        for i,h,d in ((*edges[k],1),(*edges[k][::-1],-1)):
            for h2,j,e in ((*edges[l],1),(*edges[l][::-1],-1)):
                if h==h2:
                    used=(1<<k)|(1<<l)
                    put(i,n+j,(used<<6,used,0,1),orient*d*e*(1 if k<l else -1))
    @lru_cache(None)
    def pf(indices):
        if not indices:
            return {ground:1}
        out=Counter()
        i=indices[0]
        for jidx in range(1,len(indices)):
            j=indices[jidx]
            if not entries[i][j]:
                continue
            minor=pf(indices[1:jidx]+indices[jidx+1:])
            for (m,g,x,z),c in entries[i][j].items():
                if not c:
                    continue
                for (o,h,y,t),v in minor.items():
                    if m&o or g&h or z+t>double_count or x+y not in projector.divisors:
                        continue
                    odd_swaps=sum((m >> (s+1)).bit_count() for s in range(o.bit_length()) if o&(1<<s))
                    out[m|o,g|h,x+y,z+t]+=(-1)**(jidx-1+odd_swaps)*c*v
        return {k:v for k,v in out.items() if v}
    # Evaluate the candidate trace itself on these colour edges.
    factor=0
    for order in permutations(range(n)):
        if tuple(word[j] for j in order)!=word:
            continue
        for start in edges[order[0]]:
            position=start
            weight=1
            for j in order:
                u,v=edges[j]
                if position==u:
                    position=v
                elif position==v:
                    position=u
                    weight=-weight
                else:
                    break
            else:
                if position==start:
                    factor+=weight*permutation_sign([i for i in order if odd_mask&(1<<i)])
    assert factor
    result={}
    for (m,g,x,z),v in pf(tuple(range(2*n))).items():
        if g != (1<<n)-1:
            continue
        assert m>>6==odd_mask
        mask=m&63
        # Remove the sign of Pf(J) and undo the marker/letter reordering
        # omitted from the individual log factors until the global product.
        number=mask.bit_count()
        s=(-1)**(n*(n-1)//2+number*(number-1)//2)
        result[mask,x,z]=Fraction(s*v,factor)
    return result


def adjoint_coordinate(colours, active, output_pairs):
    rotations=[]
    for k in range(len(output_pairs)):
        matrix=np.zeros((colours,colours),dtype=np.int64)
        matrix[2*k,2*k+1]=1
        matrix[2*k+1,2*k]=-1
        rotations.append(matrix)
    field_values={}
    count=0
    for i,pair in enumerate(output_pairs):
        for field in set(pair):
            label=count if is_ferm(field) else None
            count+=int(label is not None)
            field_values.setdefault(field,[]).append((label,rotations[i]))
    def assignments(letters):
        for values in product(*(field_values.get(x,[]) for x in letters)):
            indices=[label for label,_ in values if label is not None]
            if sorted(indices)==list(range(count)):
                yield permutation_sign(indices),[m for _,m in values]
    denominator=0
    target=tuple(x for pair in output_pairs for x in pair)
    for orient,m in assignments(target):
        denominator+=orient*prod_traces(m)
    assert denominator
    value=Fraction()
    for key, coefficient in primitive((active,)).items():
        assert len(key)==1
        w=key[0]
        for i,j in combinations(range(len(w)),2):
            if w[i][:2]!=(0,0) or w[j][:2]!=(0,0):
                continue
            front=is_ferm(w[i])*sum(map(is_ferm,w[:i]))+is_ferm(w[j])*sum(is_ferm(x) for k,x in enumerate(w[:j]) if k!=i)
            rest=[k for k in range(len(w)) if k not in (i,j)]
            for c,x,y in pair_rule(w[i][2],w[j][2]):
                for orientation,matrices in assignments((x,y)+tuple(w[k] for k in rest)):
                    total=0
                    xmat,ymat=matrices[:2]
                    for u,v in combinations(range(colours),2):
                        generator=np.zeros_like(xmat)
                        generator[u,v]=1
                        generator[v,u]=-1
                        left=generator@xmat-xmat@generator
                        right=generator@ymat-ymat@generator
                        if not left.any() or not right.any():
                            continue
                        positions=dict(zip([i,j]+rest,[left,right]+matrices[2:]))
                        m=positions[0]
                        for k in range(1,len(w)):
                            m=m@positions[k]
                        # The completeness in the trace kernel is the negative
                        # sum over the real skew elementary matrices.
                        total-=int(m.trace())
                    value+=coefficient*((-1)**front)*c*orientation*total
    return value/denominator


def prod_traces(matrices):
    value=1
    for i in range(0,len(matrices),2):
        value*=int((matrices[i]@matrices[i+1]).trace())
    return value


def verify(colours,power,order):
    record=highest_product(colours,power,order)
    pf=PfaffianCoefficient(order,power-order-3)
    checked=set()
    # Check every contributing active trace, all bilinear spectator types,
    # and the complete ordered alphabets of the two input lengths.
    from orthogonal_glueball_first_pair import orientations, remove_once, single_quantum
    from rank_four_lowered_replay import terminal
    letters=PHI+GHOSTS
    target=tuple(tuple(letters[j] for j in w) for w in record['target'])
    checked.update(target)
    adjoint_count=0
    for size in (1,2):
        for subset in sorted(set(combinations(target,size))):
            leaves=tuple(x for w in subset for x in w)
            if not all(g in leaves for g in GHOSTS):
                continue
            for active in orientations(remove_once(leaves,GHOSTS)+PHI):
                checked.add(active)
                s, wanted=terminal(subset)
                assert adjoint_coordinate(colours,active,subset)==s*single_quantum(colours,active).get(wanted,0)
                adjoint_count+=1
    for w in sorted(checked):
        actual=cycle_pin(w,pf.doubles,pf)
        assert actual==pf.trace(w), (w,actual,pf.trace(w))
    return dict(colours=colours,power=power,order=order,pfaffian_cycles=len(checked),
                direct_adjoint_coordinates=adjoint_count,coefficient=record['coefficient'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',nargs=3,type=int,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error('Run with assertions enabled')
    result=verify(*a.case)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(result)
