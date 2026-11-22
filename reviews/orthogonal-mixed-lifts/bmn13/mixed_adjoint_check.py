"""Check the two native adjoint equations for all three mixed terminal coordinates."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import json
from pathlib import Path
import sys


def check(column):
    from mixed_incoming import q,so_canon
    from rank_four_lowered_replay import Q1d_multi
    difference=Counter()
    for raw,a in Q1d_multi(column,COLOURS,'SO_odd').items():
        s,w=so_canon(raw)
        if s:
            for i,v in BEFORE.get(w,()):
                difference[i]+=SIGN*a*s*v
    for w,a in q({column:1}).items():
        for i,v in AFTER.get(w,()):
            difference[i]-=a*v
    assert not any(difference.values()),(column,difference)
    return 1


def stage(colours,before,after,sign,workers):
    from orthogonal_glueball_second_pair import tree_columns,quantum_columns
    global COLOURS,BEFORE,AFTER,SIGN
    COLOURS,BEFORE,AFTER,SIGN=colours,before,after,sign
    tree=tree_columns(after)
    one=quantum_columns(before)
    columns=sorted(tree|one)
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        assert sum(pool.map(check,columns,chunksize=32))==len(columns)
    return dict(tree=len(tree),one_loop=len(one),union=len(columns))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',required=True,type=Path)
    p.add_argument('--dual',required=True,type=Path)
    p.add_argument('--workers',default=1,type=int)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    assert __debug__ and a.workers>0
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from mixed_second import decode,ordered_coordinates
    d=json.loads(a.dual.read_text())
    _,target=ordered_coordinates(tuple(d['weight']))
    if 'functional' in d:
        f=d['functional']
        target={w:((0,f[i]*v),) for w,((i,v),) in target.items() if f[i]}
    first,second=decode(d['first']),decode(d['second'])
    out=dict(colours=d['colours'],weight=d['weight'],first=stage(d['colours'],target,first,1,a.workers),
             second=stage(d['colours'],first,second,-1,a.workers))
    a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(out,flush=True)
