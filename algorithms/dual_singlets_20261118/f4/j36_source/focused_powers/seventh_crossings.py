"""Full classical kernel in a family of genuine F4 level-42 Albert sources."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from fractions import Fraction
from itertools import combinations, permutations
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'global_boundaries'))
from albert_boundaries import PAIRS, SIGNS, action, apply, balanced, colour_tensors, features, native26, pure_column
from global_source_boundaries import minor
from global_source_scan import square_determinant
from integer_section import FIRST_PAGE, tensors

FIVE_PAIRS = tuple(combinations(range(5), 2))


def sign(sequence):
    return (-1) ** sum(x > y for j,x in enumerate(sequence) for y in sequence[j+1:])


def half5(corner, x, marks, axis, middle, p):
    xs = [action(v)[0] for v in x]
    total = np.zeros_like(corner)
    for z, pair in enumerate(FIVE_PAIRS):
        initial = balanced(apply(action(marks[z])[0],corner,(axis[0],)),p)
        for tail in permutations(i for i in range(5) if i not in pair):
            t=initial
            for k,leg in zip(tail,axis[1:]):
                t=balanced(apply(xs[k],t,(leg,)),p)
            total+=sign((*pair,*tail))*t
    return pure_column(balanced(total,p),xs,middle,p)


def point(number,p,free=False):
    native=native26().reshape(52,-1)
    a=(np.random.default_rng(532961+number).integers(-3,4,size=(15,52)) @ native).reshape(15,26,26)
    x=a[:5]
    b=np.array([balanced(x[i]@x[j]-x[j]@x[i],p) for i,j in FIVE_PAIRS])
    low,high,tr,tt=colour_tensors()
    def arms(side,mark):
        if side==0:
            return half5(low,x.transpose(0,2,1),mark.transpose(0,2,1),(1,0,3,2),(0,1,3),p)
        return half5(high,x,mark,(0,2,3,1),(1,2,3),p)
    def pair(left,right):
        l=features([left],tt,p)[0]
        r=features([right],tr,p)[0]
        with np.errstate(all='ignore'):
            m=l@r.T
        assert 26**4*(p//2)**2 < 2**53 and np.array_equal(m,m.astype(np.int64))
        # The first pair acts on the right tensor. Zero, one, then two exchanges.
        rounded=m.astype(np.int64)%p
        assert np.array_equal(rounded[1:,0],rounded[0,1:])
        return np.r_[rounded[0,0],rounded[0,1:],rounded[1:,1:].T.ravel()]
    lb,rb=arms(0,b),arms(1,b)
    if not free:
        return pair(lb,rb)
    left=pair(arms(0,a[5:]),rb)
    right=pair(lb,arms(1,a[5:]))
    return np.r_[left,right]


def check_relations(values):
    for i,j in combinations(range(6),2):
        if set(PAIRS[i]).isdisjoint(PAIRS[j]):
            assert np.array_equal(values[...,7+6*i+j],values[...,7+6*j+i])


def colour_relations():
    # A symmetric traceless pair of open 26 indices sees only eigenvalues -36 and 6.
    exchange=tensors()[2]
    native=exchange.astype(np.int64)
    assert np.array_equal(native,exchange)
    identity=np.eye(676,dtype=np.int64)
    assert 676*int(np.abs(native).max())**2+30*int(np.abs(native).max())+216<2**63
    operator=(native@native+30*native-216*identity).reshape(676,26,26)
    raw=json.loads(gzip.decompress((FIRST_PAGE/'representation26.json.gz').read_bytes()))
    metric=np.zeros((26,26),dtype=np.int64)
    for i,j,v in raw['form']:
        q=Fraction(v)
        assert q.denominator==1
        metric[i,j]=q.numerator
    inverse_times_three=3*metric.copy()
    inverse_times_three[:2,:2]=[[2,1],[1,2]]
    assert np.array_equal(metric@inverse_times_three,3*np.eye(26,dtype=np.int64))
    symmetric=(operator+operator.transpose(0,2,1)).reshape(676,676)
    assert np.array_equal(symmetric,240*np.outer(inverse_times_three.ravel(),metric.ravel()))
    # Tensor products in the fundamental-weight coordinates of the retained Weyl calculation.
    sys.path.insert(0,str(HERE.parents[2]))
    from f4_singlet_characters import decomposition, multiply
    product=decomposition('F4',(1,1))
    doubled=decomposition('F4',(2,))
    differences={w:product.get(w,0)-doubled.get(w,0) for w in set(product)|set(doubled)}
    assert all(n%2==0 for n in differences.values())
    exterior={w:n//2 for w,n in differences.items() if n}
    mixed=dict(multiply('F4',(0,0,0,1),1))
    assert exterior=={(0,1,0,0):1,(1,0,0,0):1}
    assert mixed=={(0,0,0,1):1,(0,0,1,0):1,(1,0,0,1):1}
    assert not (set(exterior)&set(mixed))
    return True


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prime',type=int,choices=(1009,1013),default=1009)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--fixed',type=Path)
    args=ap.parse_args()
    p=args.prime
    open_values=point(3,p,True).reshape(2,43)
    check_relations(open_values)
    assert colour_relations()
    assert not np.any((open_values[:,35]+30*open_values[:,5]-216*open_values[:,0])%p)
    assert np.any((open_values[0]-open_values[1])%p)
    fixed=json.loads(args.fixed.read_bytes()) if args.fixed else None
    ids=fixed['points'] if fixed else list(range(64))
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        a=np.stack(list(pool.map(partial(point,p=p),ids,chunksize=1)))
    check_relations(a)
    assert not np.any((a[:,35]+30*a[:,5]-216*a[:,0])%p)
    # Both marks map to the same classical tensor by direct substitution.
    if fixed:
        col=fixed['columns']
        d=square_determinant(a[:,col],p)
        rows=list(range(len(ids)))
        assert d
    else:
        with np.errstate(all='ignore'):
            rows,col,d=minor(a,p)
    out=dict(prime=p,points=[ids[i] for i in rows],columns=col,rank=len(rows),determinant=int(d),
             double_exchange_columns=36,disjoint_commutations=3,
             left_minus_right_nonzero_coordinates=int(np.count_nonzero((open_values[0]-open_values[1])%p)))
    print(json.dumps(out),flush=True)
    if args.output:
        args.output.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
