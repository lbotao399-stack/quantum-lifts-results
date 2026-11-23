"""Classical endpoint evaluations with three ordered Albert exchanges at level 42."""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import product
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from seventh_crossings import FIVE_PAIRS, PAIRS, action, apply, balanced, colour_tensors, half5, native26
from global_source_boundaries import minor
from global_source_scan import square_determinant


def evaluate(n,p,free=False):
    native=native26().reshape(52,-1)
    colours=(np.random.default_rng(532961+n).integers(-3,4,size=(15,52)) @ native).reshape(15,26,26)
    x=colours[:5]
    marks=np.array([balanced(x[i]@x[j]-x[j]@x[i],p) for i,j in FIVE_PAIRS])
    low,high,ex,et=colour_tensors()
    actual=colours[5:]
    def arm(side,y):
        if side:
            return half5(high,x,y,(0,2,3,1),(1,2,3),p)
        return half5(low,x.transpose(0,2,1),y.transpose(0,2,1),(1,0,3,2),(0,1,3),p)
    def ordered(left,right):
        l=np.stack([left.ravel(),*[balanced(apply(et,left,ab),p).ravel() for ab in PAIRS]])
        one=[balanced(apply(ex,right,ab),p) for ab in PAIRS]
        r0=right.ravel()
        def pair(r):
            with np.errstate(all='ignore'):
                scores=l@r.T
            assert 26**4*(p//2)**2<2**53 and np.array_equal(scores,scores.astype(np.int64))
            return scores.astype(np.int64)%p
        first=pair(r0[None,:])[:,0]
        out=[first[:1],first[1:]]
        order2=pair(np.stack([v.ravel() for v in one]))
        assert np.array_equal(first[1:],order2[0])
        out.append(order2[1:].T.ravel())
        for v in one:
            rr=np.stack([balanced(apply(ex,v,cd),p).ravel() for cd in PAIRS])
            out.append(pair(rr)[1:].T.ravel())
        return np.concatenate(out).astype(np.uint16)
    left,right=arm(0,marks),arm(1,marks)
    if free:
        return np.stack([ordered(arm(0,actual),right),ordered(left,arm(1,actual))])
    return ordered(left,right)


def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--prime',type=int,default=1009);ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--fixed',type=Path)
    args=ap.parse_args()
    p=args.prime
    fixed=json.loads(args.fixed.read_text()) if args.fixed else None
    ids=fixed['points'] if fixed else list(range(260))
    # Free rows test for classical syzygies beyond relations of the colour insertion.
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        rows=np.stack(list(pool.map(partial(evaluate,p=p),ids,chunksize=1)))
        extra=np.concatenate(list(pool.map(partial(evaluate,p=p,free=True),range(400,418),chunksize=1)))
    if fixed:
        rr=list(range(len(ids))); cc=fixed['columns']
        d=square_determinant(rows[:,cc],p)
    else:
        with np.errstate(all='ignore'):
            rr,cc,d=minor(rows,p)
    from mixed_source_search import solve_left
    coords,det=solve_left(rows[np.ix_(rr,cc)],extra[:,cc],p)
    assert d==det
    remainder=(extra.astype(np.int64)-coords.astype(np.int64) @ rows[rr].astype(np.int64))%p
    with np.errstate(all='ignore'):
        nr,nc,dd=minor(remainder,p)
    out=dict(prime=p,columns=list(cc),points=[ids[i] for i in rr],rank=len(rr),determinant=int(d),
             number_of_recipes=259,extra_free_rank=len(nr),extra_free_minor=int(dd),
             extra_rows=nr,extra_columns=nc)
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    np.savez_compressed(args.output.with_suffix('.npz'),closed=rows[rr].astype(np.uint16),free=extra)
    print('rank',len(rr),'minor',d,'extra free',len(nr),'minor',dd,flush=True)

if __name__=='__main__':
    main()
