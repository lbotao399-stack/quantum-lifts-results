"""Finite-colour path relations and second scores for one scalar bilinear at the target."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from math import factorial,prod
from multiprocessing import get_context
from pathlib import Path
import sys


def mixed_projector(axes,order=2):
    from orthogonal_glueball_first_pair import PfaffianCoefficient
    projection=PfaffianCoefficient(order,0)
    weights={}
    divisors=set()
    for degree,value in projection.weights.items():
        counts=[(degree >> (3*i)) & 7 for i in range(9)]
        for i in axes:
            counts[i]+=1
            value*=counts[i]
        code=sum(x<<(3*i) for i,x in enumerate(counts))
        weights[code]=value
        subset={0}
        for i,n in enumerate(counts):
            for _ in range(n):
                subset.update([v+(1 << (3*i)) for v in subset])
        divisors.update(subset)
    projection.weights=weights
    projection.divisors=divisors
    return projection


def charge(spec,axes):
    weights=(0,1,-1,0)
    single=sum(weights[i] for row in spec['single'] for i in row)
    scalar=sum(weights[i] for i in axes)
    scalar-=sum(weights[v.index(1)] for part in spec['scalar_directions'] for v in part)
    return single+scalar+sum(weights[m.index(1)] for m in spec['matrices'])


def evaluate(task):
    from orthogonal_glueball_relation_span import determinant
    index,spec,axes=task
    relation=determinant(POWER-len(spec['closed']),spec)
    relation.project=mixed_projector(axes)
    relation.pf_degree+=2
    assert relation.pf_degree>COLOURS
    relation.closed=bool(spec['closed'])
    coords=[relation.cartan(w) for _,w in PLANES]
    scores=[Fraction() for _ in range(TARGET_DIM)]
    for w,entries in DUAL:
        coefficient=relation.coefficient(w)
        for i,v in entries:
            scores[i]+=coefficient*v
    relation.clear()
    return dict(index=index,axes=axes,spec=spec,planes=list(map(str,coords)),scores=list(map(str,scores)))


def check(rows):
    from flint import fmpq_mat
    cols=fmpq_mat([r['planes'] for r in rows])
    aug=fmpq_mat([r['planes']+r['scores'] for r in rows])
    rank=cols.rank()
    result=dict(relation_plane_rank=rank,augmented_rank=aug.rank())
    if result['augmented_rank']>rank:
        rref,rank=cols.transpose().rref()
        pivots=[next(j for j in range(rref.ncols()) if rref[i,j]) for i in range(rank)]
        for j,row in enumerate(rows):
            if j in pivots:
                continue
            c=[Fraction(0) for _ in rows]
            c[j]=1
            for i,k in enumerate(pivots):
                c[k]=-Fraction(str(rref[i,j]))
            assert all(sum(x*Fraction(r['planes'][z]) for x,r in zip(c,rows))==0 for z in range(cols.ncols()))
            v=[sum(x*Fraction(r['scores'][z]) for x,r in zip(c,rows)) for z in range(TARGET_DIM)]
            if any(v):
                result['witness']={str(r['index']):str(x) for x,r in zip(c,rows) if x}
                result['image']=list(map(str,v))
                break
        assert 'witness' in result
    return result


def main():
    global COLOURS,POWER,DUAL,PLANES,TARGET_DIM
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--dual',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('Use assertions and a positive worker count')
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from second_path_quotient import polarizations
    from terminal_commutators import words
    d=json.loads(args.dual.read_text())
    COLOURS,POWER=d['colours'],d['weight'][0]
    assert d['weight']==[POWER,POWER,2,0,0]
    assert COLOURS>=POWER+7 and POWER>=6
    DUAL=[(tuple(map(tuple,w)),[(i,Fraction(v)) for i,v in row]) for w,row in d['second']]
    TARGET_DIM=len(d['target'])
    basis=words((POWER-2,POWER-2,4,2,2))
    PLANES=[(i,w) for i,w in enumerate(basis) if len(w)>COLOURS//2]
    tasks=[(i,spec,axes) for i,(spec,axes) in enumerate(
        (spec,axes) for spec in polarizations(COLOURS-2,POWER)
        for axes in combinations_with_replacement(range(3),2)) if not charge(spec,axes)]
    data=dict(colours=COLOURS,weight=d['weight'],target=d['target'],source_plane_indices=[i for i,_ in PLANES],
              all_source_coordinates=len(basis),polarizations=len(tasks),rows=[])
    if args.out.exists():
        other=json.loads(args.out.read_text())
        assert all(other[k]==v for k,v in data.items() if k!='rows')
        data=other
    print('nonzero-weight-allowed polarizations',len(tasks),'source coordinates',len(basis),'above rank',len(PLANES),flush=True)
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        n=len(data['rows'])
        for start in range(n,len(tasks),4*args.workers):
            data['rows'].extend(pool.map(evaluate,tasks[start:start+4*args.workers]))
            data.update(check(data['rows']))
            tmp=args.out.with_suffix('.part')
            tmp.write_text(json.dumps(data,separators=(',',':'))+'\n')
            tmp.replace(args.out)
            print('rows',len(data['rows']),'ranks',data['relation_plane_rank'],data['augmented_rank'],flush=True)
            if 'witness' in data:
                break

if __name__=='__main__':
    main()
