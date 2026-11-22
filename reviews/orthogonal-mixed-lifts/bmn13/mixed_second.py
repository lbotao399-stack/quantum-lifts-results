"""Pull the second obstruction back from a mixed terminal orthogonal weight."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import product
import json
from multiprocessing import get_context
from pathlib import Path
import sys


def ordered_coordinates(weight):
    from orthogonal_glueball_first_pair import GHOSTS, PHI, orientations
    from mixed_incoming import so_canon
    from rank_four_lowered_replay import terminal
    from terminal_commutators import words
    basis=words(weight)
    rows={}
    for i,monomial in enumerate(basis):
        choices=[orientations(tuple(a for a,n in zip(GHOSTS+PHI,t) for _ in range(n))) for t in monomial]
        for traces in product(*choices):
            s,w=so_canon(traces)
            if s and w not in rows:
                sign,check=terminal(w)
                assert sign and check==monomial
                rows[w]=((i,sign),)
    return basis,rows


def pull(word):
    from mixed_incoming import q,so_canon
    from orthogonal_glueball_second_pair import tree_split
    from rank_four_lowered_replay import Q1d_multi
    primitive,section=tree_split(word)
    expected=Counter({word:1})
    expected.subtract(section)
    assert q(primitive)=={w:v for w,v in expected.items() if v}
    values=Counter()
    for parent,a in primitive.items():
        if not any(len(l[2])==2 for trace in parent for l in trace):
            continue
        for output,b in Q1d_multi(parent,COLOURS,'SO_odd').items():
            sign,w=so_canon(output)
            if sign:
                for i,v in DUAL.get(w,()):
                    values[i]+=SIGN*a*b*v*sign
    return word,tuple(sorted((i,v) for i,v in values.items() if v))


def stage(colours,dual,sign,workers):
    global COLOURS,DUAL,SIGN
    from orthogonal_glueball_second_pair import predecessors
    COLOURS,DUAL,SIGN=colours,dual,sign
    support=predecessors(dual)
    print('possible ordered parents',len(support),flush=True)
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        result={w:row for w,row in pool.map(pull,support,chunksize=8) if row}
    print('nonzero ordered parents',len(result),flush=True)
    return result,len(support)


def encode(dual):
    from orthogonal_glueball_third_pair import encode as encode_word
    return [[encode_word(w),[[i,str(v)] for i,v in row]] for w,row in sorted(dual.items())]


def decode(rows):
    from orthogonal_glueball_third_pair import decode as decode_word
    return {decode_word(w):tuple((i,Fraction(v)) for i,v in row) for w,row in rows}


def run(colours,power,scalars,workers,out,functional=None):
    weight=(power,power,*scalars)
    basis,target=ordered_coordinates(weight)
    if functional is not None:
        if len(functional)!=len(basis):
            raise ValueError('functional length differs from the stable terminal dimension')
        target={w:((0,functional[i]*value),) for w,((i,value),) in target.items() if functional[i]}
    cache=out.with_suffix('.first.json')
    if cache.exists():
        data=json.loads(cache.read_text())
        assert [colours,list(weight)]==[data['colours'],data['weight']]
        assert data.get('functional')==functional
        first=decode(data['first'])
    else:
        first,n=stage(colours,target,1,workers)
        data=dict(colours=colours,weight=weight,target=basis,target_ordered=len(target),
                  first_candidates=n,first=encode(first))
        if functional is not None:
            data['functional']=functional
        cache.write_text(json.dumps(data,separators=(',',':'))+'\n')
    second,n=stage(colours,first,-1,workers)
    data.update(second_candidates=n,second=encode(second))
    out.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print('target dimension',len(basis),'ordered',len(target),'first and second',len(first),len(second),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--colours',type=int,required=True)
    parser.add_argument('--power',type=int,required=True)
    parser.add_argument('--scalars',nargs=3,type=int,default=[2,0,0])
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--functional',type=int,nargs='+')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('Use assertions and a positive worker count')
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    run(args.colours,args.power,args.scalars,args.workers,args.out,args.functional)
