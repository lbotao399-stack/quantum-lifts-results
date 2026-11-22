"""Finite orthogonal relations modulo internal tree commutators at pure glueball weight.

Pass --repo pointing to the quantum-lifts-results checkout at ee9e064.
The supported cases have at most one closed spectator, the two-gaugino trace.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations, product
import json
from multiprocessing import get_context
from pathlib import Path
import sys


def basis_vector(length, axis):
    return [int(i == axis) for i in range(length)]


def path_types(n, power):
    """All multisets of non-single paths after commutation inside each path."""
    limit = 2*power+1-n
    catalogue = []
    for p, m, s in product(range(2), range(2), range(7)):
        length = p+m+s
        if length < 2:
            continue
        name = 'm'*m + 'p'*p + 's'*s
        for closed in (False, True):
            cost = length if closed else length-1
            if cost <= limit and (not closed or length % 2 == 0):
                catalogue.append((closed, name, cost, p, m, s))
    catalogue.sort()
    out = []
    def visit(start, cost, ps, ms, ss, paths):
        if min(ps,ms,ss) < 0 or max(0,ps-3)+max(0,ms-3) > 2*cost:
            return
        if max(ps,ms) <= 3:
            out.append((tuple(paths), ps, ms))
        for i in range(start,len(catalogue)):
            closed,name,used,p,m,s=catalogue[i]
            if used<=cost:
                visit(i,cost-used,ps-p,ms-m,ss-s,paths+[(closed,name)])
    visit(0,limit,power-2,power-2,6,[])
    return out


def polarizations(n, power):
    for paths, plus, minus in path_types(n,power):
        closed = [p for c,p in paths if c]
        if closed not in ([], ['mp']):
            raise ValueError(f'Unsupported closed spectator {closed}')
        opens = [p for c,p in paths if not c]
        for axes in product(*(range(3) if len(p)%2 else [3] for p in opens)):
            for scalar in product(range(3),repeat=sum(p.count('s') for p in opens)):
                it=iter(scalar)
                directions=[[basis_vector(3,next(it)) for _ in range(p.count('s'))] for p in opens]
                for left,right in product(combinations(range(3),plus),combinations(range(3),minus)):
                    yield {'closed':closed, 'paths':opens, 'matrices':[basis_vector(4,i) for i in axes],
                           'single':[left,right], 'scalar_directions':directions}


def evaluate(spec):
    from orthogonal_glueball_relation_span import determinant
    inside = POWER-len(spec['closed'])
    relation=determinant(inside,spec)
    assert relation.pf_degree > COLOURS
    # This flag adds the single canonical spectator only after constructing the
    # marked open paths. It leaves their multivariate coefficient unchanged.
    relation.closed = bool(spec['closed'])
    second=sum(c*relation.coefficient(w) for w,c in SCORES)
    coords=[relation.cartan(w) if len(w)>COLOURS//2 else 0 for w in BASIS]
    relation.clear()
    return dict(spec=spec,cartan=[str(x) for x in coords],second_score=str(second))


def main():
    global COLOURS, POWER, SCORES, BASIS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--colours',type=int,default=13)
    parser.add_argument('--power',type=int,default=8)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('Run with assertions and positive worker count')
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from flint import fmpq, fmpq_mat
    from orthogonal_glueball_relation_span import relation_dimension
    from orthogonal_glueball_second_pair import second_dual
    from orthogonal_glueball_second_sources import ranks
    from orthogonal_glueball_sources import physical_image_relations
    from terminal_commutators import words, audit_rules
    COLOURS,POWER=args.colours,args.power
    assert len(audit_rules())==12
    lower=[relation_dimension(COLOURS,POWER,k,False) for k in (0,1)]
    assert lower==[0,0]
    dual=second_dual(COLOURS,POWER)
    SCORES=[(tuple(map(tuple,w)),Fraction(v)) for w,v in dual['scores']]
    BASIS=words((POWER-2,POWER-2,2,2,2))
    specs=list(polarizations(COLOURS,POWER))
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        rows=list(pool.map(evaluate,specs))
    result=ranks(rows)
    matrix=fmpq_mat([row['cartan'] for row in rows])
    reduced,rank=matrix.rref()
    columns=[next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]
    reduced_rows,row_rank=matrix.transpose().rref()
    indices=[next(j for j in range(reduced_rows.ncols()) if reduced_rows[i,j]) for i in range(row_rank)]
    minor=fmpq_mat([[matrix[i,j] for j in columns] for i in indices])
    if result['augmented_rank']==rank:
        coefficients=minor.solve(fmpq_mat([[rows[i]['second_score']] for i in indices]))
        assert all(sum(matrix[i,j]*coefficients[k,0] for k,j in enumerate(columns))==fmpq(row['second_score'])
                   for i,row in enumerate(rows))
        result['functional']=[str(coefficients[k,0]) for k in range(rank)]
    result.update(group=f'SO({COLOURS})',power=POWER,lower=lower,balanced_dimension=relation_dimension(COLOURS,POWER,2),
                  singlet_dimension=physical_image_relations(POWER,2,COLOURS),types=len(path_types(COLOURS,POWER)),
                  polarizations=len(specs),stable_plane_coordinates=len(BASIS),
                  factorization_columns=[dual['first_factorization_columns'],dual['second_factorization_columns']],
                  pivot_rows=indices,pivot_columns=columns,minor=str(minor.det()),
                  rows=[dict(spec=r['spec'],nonzero_planes={str(j):c for j,c in enumerate(r['cartan']) if c!='0'},
                             score=r['second_score']) for r in rows])
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='rows'},flush=True)


if __name__=='__main__':
    main()
