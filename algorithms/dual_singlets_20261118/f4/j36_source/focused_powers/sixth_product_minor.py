"""Rebuild the fixed 3017-column classical kernel minor for the sixth F4 glueball."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import permutations
import json
from multiprocessing import get_context
from pathlib import Path

import numpy as np

from albert_product_cycles import AXES, PAIRS, OMITTED, action, alternating, balanced, colour_tensors, highest, point, samples, scalar_polynomials
from sixth_boundary_minor import central, evaluate, native_old, previous, record, solve_left


HERE = Path(__file__).resolve().parent

def product():
    cols, coeff = highest()
    active = np.flatnonzero(np.any(coeff, axis=0))
    assert len(cols) == 96 and len(coeff) == 1 and len(active) == 87
    return [cols[i] for i in active], coeff[:,active]


def controls(prime):
    cols, coeff = product()
    x, _, marks = samples(32,prime)
    lower, upper = colour_tensors()[:2]
    weights = (1,-2,3,1)
    matrices = [action(t)[0] for t in x]
    value = upper
    for _ in range(3):
        value = balanced(sum(w*alternating(value,matrices,i,prime) for i,w in enumerate(weights)),prime)
    direct = int(lower.astype(np.int64).ravel() @ value.astype(np.int64).ravel()) % prime
    scalars = scalar_polynomials(x,prime,OMITTED)
    expected = sum(len(set(permutations(t)))*int(np.prod([weights[i] for i in t]))*s for t,s in scalars.items())
    assert 6*direct%prime == expected%prime
    def contract(xs, ys):
        terms = scalar_polynomials(xs,prime,{t[0] for t in cols})
        gram = xs.reshape(4,-1) @ xs.transpose(0,2,1).reshape(4,-1).T % prime
        marked = ys.reshape(6,-1) @ xs.transpose(0,2,1).reshape(4,-1).T % prime
        return int(sum(int(a)*terms[t]*int(gram[g])*int(marked[PAIRS.index(pair),k])
                       for a,(t,g,pair,k) in zip(coeff[0],cols))%prime)
    original = int(point(32,prime,cols,coeff)[0])
    assert original and original == contract(x,marks)
    for i in range(3):
        xs,ys = x.copy(),marks.copy()
        xs[i+1] += 2*x[i]
        for k,pair in enumerate(PAIRS):
            if i+1 in pair and i not in pair:
                other = tuple(i if a==i+1 else a for a in pair)
                ys[k] += 2*marks[PAIRS.index(other)]
        assert contract(xs,ys) == original
    assert not point(32,prime,cols,coeff,closed=True).any()
    return dict(prime=prime,free_value=original,split_value=0,
                highest_weight_shears=3,direct_cubic_polarization=direct)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prime',type=int,choices=(1009,1013),required=True)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--controls',action='store_true')
    args=ap.parse_args()
    p=args.prime
    if args.controls:
        print(json.dumps(controls(p),indent=2))
        return
    fixed=json.loads((HERE/f'sixth-product-{p}.json').read_bytes())
    ids=previous()['rows']+central()['new_rows']+fixed['new_points']
    assert len(ids)==len(set(ids))==3017
    cols,v=product()
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        old=np.stack(list(pool.map(partial(native_old,prime=p),ids,chunksize=2)))
        weights,base=solve_left(old[:3002],old[3002:],p)
        assert np.array_equal(weights.astype(np.int64)@old[:3002].astype(np.int64)%p,old[3002:])
        assert base==fixed['previous_minor']
        reference=np.stack(list(pool.map(partial(evaluate,prime=p),ids,chunksize=2)))
        new=np.stack(list(pool.map(partial(point,p=p,columns=cols,coeff=v),ids,chunksize=2)))
    assert reference.shape==(3017,49) and new.shape==(3017,1)
    checked=record(weights,np.column_stack((reference,new)),ids,int(base),p,fixed)
    assert checked==fixed
    print(json.dumps(checked,indent=2))

if __name__=='__main__':
    main()
