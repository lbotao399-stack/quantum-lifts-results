"""Independent combinatorial and trace-coordinate checks for second_path_quotient.py."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import permutations, product
import json
from multiprocessing import get_context
from pathlib import Path
import sys


def raw_graph_types(colours, power):
    """Enumerate multisets of open paths and closed cycles before internal sorting."""
    budget = 2*power+1-colours
    catalogue = []
    for length in range(2, budget+2):
        for word in product("mps", repeat=length):
            word = ''.join(word)
            if word <= word[::-1]:
                catalogue.append((False, word, length-1))
            cyclic = {s[k:]+s[:k] for s in (word, word[::-1]) for k in range(length)}
            if word == min(cyclic) and length <= budget:
                catalogue.append((True, word, length))
    records = Counter()
    reduced = set()
    def enumerate_graphs(start, remaining, counts, graphs):
        m, p, s = counts
        if m <= 3 and p <= 3:
            records['raw'] += 1
            reasons = set()
            for closed, word in graphs:
                if max(word.count('m'), word.count('p')) > 1:
                    reasons.add('repeated_odd')
                if closed and len(word)%2:
                    reasons.add('odd_cycle')
            if reasons:
                records['discarded'] += 1
            else:
                reduced.add((tuple(sorted((closed,''.join(sorted(w))) for closed,w in graphs)),p,m))
        for i in range(start,len(catalogue)):
            closed, word, cost = catalogue[i]
            after = (m-word.count('m'), p-word.count('p'), s-word.count('s'))
            if cost <= remaining and min(after) >= 0:
                enumerate_graphs(i, remaining-cost, after, graphs+[(closed,word)])
    enumerate_graphs(0,budget,(power-2,power-2,6),[])
    return records,reduced


def flavour_projector():
    from flint import fmpz_mat
    from orthogonal_glueball_first_pair import PfaffianCoefficient, scalar_code
    projector = PfaffianCoefficient(2,0)
    balanced = sorted(set(permutations((0,0,1,1,2,2))))
    matrix = fmpz_mat([[projector.weights.get(scalar_code([3*i+j for i,j in zip(a,b)]),0)
                        for b in balanced] for a in balanced])
    assert matrix == matrix.transpose() and matrix*matrix == 144*matrix
    assert matrix.rank() == 5
    return dict(balanced_words=len(balanced),rank=5,square_multiplier=144)


def trace_check(row):
    from flint import fmpq_mat
    from orthogonal_glueball_relation_span import determinant
    spec = row['spec']
    relation = determinant(POWER-len(spec['closed']),spec)
    assert relation.pf_degree > COLOURS
    relation.closed = bool(spec['closed'])
    traces = fmpq_mat([[str(relation.commuting_coefficient(w)) for w in BASIS]])
    expected = fmpq_mat([[row['nonzero_planes'].get(str(i),'0') for i in range(len(BASIS))]])
    assert traces*CHANGE == expected
    relation.clear()
    return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--uncompressed',type=Path)
    parser.add_argument('--rebuild-uncompressed',action='store_true')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if not __debug__ or args.workers<1 or (args.rebuild_uncompressed and not args.uncompressed):
        parser.error('Use assertions and positive worker count')
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from orthogonal_glueball_cartan_basis import matrix
    from second_path_quotient import path_types
    from terminal_commutators import words
    case = json.loads(args.case.read_text())
    global POWER,COLOURS,BASIS,CHANGE
    POWER,COLOURS=case['power'],int(case['group'][3:-1])
    counts,reduced=raw_graph_types(COLOURS,POWER)
    assert reduced == set(path_types(COLOURS,POWER))
    assert set(w for paths,_,_ in reduced for closed,w in paths if closed) <= {'mp'}
    assert max(sum(w.count('s') for _,w in paths) for paths,_,_ in reduced) == 2
    projection = flavour_projector()
    BASIS=words((POWER-2,POWER-2,2,2,2))
    CHANGE=matrix(POWER)
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        verified=sum(pool.map(trace_check,case['rows']))
    out=dict(group=case['group'],power=POWER,graph_counts=dict(counts),reduced_types=len(reduced),
             flavour_projection=projection,reconstructed_trace_rows=verified,planes_per_row=len(BASIS))
    if args.uncompressed:
        uncompressed=json.loads(args.uncompressed.read_text())
        assert (uncompressed['colours'],uncompressed['power'])==(COLOURS,POWER)
        if args.rebuild_uncompressed:
            from orthogonal_glueball_even_colours import evaluate
            from orthogonal_glueball_second_pair import second_dual
            dual=second_dual(COLOURS,POWER)
            assert dual==uncompressed['dual']
            scores=[(tuple(map(tuple,w)),Fraction(v)) for w,v in dual['scores']]
            tasks=[(r['index'],r['spec'],COLOURS,POWER,scores) for r in uncompressed['rows']]
            with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
                for rebuilt,expected in zip(pool.map(evaluate,tasks),uncompressed['rows'],strict=True):
                    assert rebuilt==expected
        functional=list(zip(case['pivot_columns'],map(Fraction,case['functional'])))
        for r in uncompressed['rows']:
            assert sum(Fraction(r['cartan'][i])*v for i,v in functional)==Fraction(r['second_score'])
        out.update(uncompressed_rows=len(uncompressed['rows']),nonzero_uncompressed_scores=
                   sum(bool(Fraction(r['second_score'])) for r in uncompressed['rows']))
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(out,flush=True)


if __name__=='__main__':
    main()
