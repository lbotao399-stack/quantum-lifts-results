"""Pull the lower-page rows and the third pairing back to determinant coefficients."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

import mixed_incoming as traces
from orthogonal_glueball_first_pair import GHOSTS, PHI, orientations
from orthogonal_glueball_second_pair import predecessors, tree_split
from rank_four_lowered_replay import Q1d_multi, terminal
from terminal_commutators import words


def encode(w):
    return tuple(tuple((PHI+GHOSTS).index(l) for l in t) for t in w)


def decode(w):
    return tuple(tuple((PHI+GHOSTS)[i] for i in t) for t in w)


def targets(power,dual):
    basis=words((power-2,power-2,2,2,2))
    lookup={}
    for i,m in enumerate(basis):
        choices=[]
        for w in m:
            leaves=tuple(x for x,n in zip(GHOSTS+PHI,w) for _ in range(n))
            choices.append(orientations(leaves))
        for ts in product(*choices):
            s,v=traces.so_canon(ts)
            if s and v not in lookup:
                sign, check=terminal(v)
                assert sign and check==m
                lookup[v]=(i,sign)
    score={decode(w):Fraction(c) for w,c in dual['scores']}
    return basis,lookup,score


def rows(colours,selected,lookup,score):
    result=[]
    for w in selected:
        h,section=tree_split(w)
        assert not section
        coordinates=Counter()
        third=Fraction()
        for pre,c in h.items():
            if not any(len(x[2])==2 for t in pre for x in t):
                continue
            for raw,n in Q1d_multi(pre,colours,'SO_odd').items():
                s,v=traces.so_canon(raw)
                if not s:
                    continue
                if v in lookup:
                    i,k=lookup[v]
                    coordinates[i]+=c*n*s*k
                if v in score:
                    third-=c*n*s*score[v]
        coeff=[[i,str(n)] for i,n in sorted(coordinates.items()) if n]
        if coeff or third:
            result.append([encode(w),coeff,str(third)])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dual',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--list',action='store_true')
    p.add_argument('--range',nargs=2,type=int)
    a=p.parse_args()
    if not __debug__:
        p.error('Run with assertions enabled')
    d=json.loads(a.dual.read_bytes())
    basis,lookup,scores=targets(d['power'],d)
    support=predecessors(set(lookup)|set(scores))
    print('commuting, ordered, dual and predecessor counts:',len(basis),len(lookup),len(scores),len(support),flush=True)
    if a.list:
        a.out.write_text(json.dumps([encode(w) for w in support],separators=(',',':'))+'\n')
    else:
        selected=support[slice(*a.range)] if a.range else support
        data=rows(d['colours'],selected,lookup,scores)
        a.out.write_text(json.dumps(data,separators=(',',':'))+'\n')
