"""Finite-rank highest weights, certified by exact character dimensions.

Project canonical trace words with the SU(3) extremal projector, and retain
independent evaluations until the complete exact dimension is reached. This
spans physical cochains, not the universal closed formal-trace space.
"""
import hw_extremal as projection
import sp_census as census
import argparse
import hashlib
import json
import pathlib
import pickle
import random
import signal
import time
from math import ceil, comb
from collections import defaultdict
base = census.base
np = base.np


def evaluate(expressions, frame):
    words = sorted({w for e in expressions for w in e})
    if not words:
        return np.zeros((len(expressions), frame.ncols), dtype=np.int64)
    index = {word:i for i,word in enumerate(words)}
    rr, cc, vv = [], [], []
    for i, expr in enumerate(expressions):
        for word, coefficient in expr.items():
            rr.append(i)
            cc.append(index[word])
            vv.append(int(coefficient)%frame.p)
    assert max(map(len,expressions), default=0)*(frame.p-1)**2 < 2**63
    matrix = base.csr_matrix((np.array(vv, dtype=np.int64),(rr,cc)), shape=(len(expressions),len(words)))
    return np.asarray(matrix@frame.values(words))%frame.p


def independent_rows(values, prime):
    matrix = base.nmod_mat(values.shape[1], values.shape[0], prime)
    for j, row in enumerate(values):
        for i, coefficient in enumerate(row.tolist()):
            matrix[i,j] = coefficient
    reduced, rank = matrix.rref()
    return [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]


def build(sector, rank, prime, output=None, batch_size=32):
    start = time.monotonic()
    census.configure(rank, prime)
    exact = base.exact_hw(sector,0)
    if not exact:
        return [], dict(sector=sector, dimension=0, exact_dimensions={}, seconds=time.monotonic()-start)
    words = list(base.sector_basis(sector,0))
    random.Random(38711).shuffle(words)
    frames, chosen, values, pending, seen = {}, {}, {}, {}, {}
    for degree, dimension in exact.items():
        k = degree+2 if degree else 0
        frames[degree] = base.Frame(degree, ceil((dimension+16)/comb(k,degree)),
                                    125131+7919*degree, prime, k=k)
        chosen[degree], pending[degree], seen[degree] = [], [], set()
        values[degree] = np.zeros((0, frames[degree].ncols), dtype=np.int64)
    counts = dict(words_examined=0, nonzero_projections=0, distinct_projections=0)

    def flush(degree):
        if not pending[degree]:
            return
        new = evaluate(pending[degree],frames[degree])
        joint = np.concatenate([values[degree],new],axis=0)
        ids = independent_rows(joint,prime)
        expressions = chosen[degree]+pending[degree]
        chosen[degree] = [expressions[i] for i in ids]
        values[degree] = joint[ids]
        pending[degree] = []
        assert len(ids) <= exact[degree]

    last = start
    for word in words:
        degree = sum(base.count(word)[3:6])
        if degree not in exact or len(chosen[degree]) == exact[degree]:
            continue
        counts['words_examined'] += 1
        expr = projection.project({word:1},sector,prime)
        if not expr:
            continue
        counts['nonzero_projections'] += 1
        inverse = pow(expr[min(expr)],-1,prime)
        expr = {w:c*inverse%prime for w,c in expr.items()}
        fingerprint = tuple(sorted(expr.items()))
        if fingerprint in seen[degree]:
            continue
        seen[degree].add(fingerprint)
        counts['distinct_projections'] += 1
        assert not projection.action(expr,0,1,prime)
        assert not projection.action(expr,1,2,prime)
        pending[degree].append(expr)
        if len(pending[degree]) >= batch_size:
            flush(degree)
        if time.monotonic()-last > 15:
            print(json.dumps(dict(stage='projected_basis',sector=sector,
                                  dimensions={d:len(v) for d,v in chosen.items()},
                                  exact_dimensions=exact,**counts,
                                  seconds=time.monotonic()-start)),flush=True)
            last = time.monotonic()
        if all(len(chosen[d]) == dim for d,dim in exact.items()):
            break
    for degree in exact:
        flush(degree)
        assert len(chosen[degree]) == exact[degree], ('Incomplete projected basis',degree,len(chosen[degree]),exact[degree])
    result = [expr for degree in sorted(chosen) for expr in chosen[degree]]
    record = dict(sector=sector, group='Sp(%d)'%rank, prime=prime,
                  dimension=len(result), exact_dimensions=exact,
                  evaluation_columns={d:fr.ncols for d,fr in frames.items()},
                  canonical_weight_words=len(words), **counts,
                  seconds=time.monotonic()-start,
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                 (pathlib.Path(__file__),pathlib.Path(projection.__file__))},
                  scope='Complete finite-rank cochain basis over the stated field; universal formal cohomology is not computed.')
    if output:
        output = pathlib.Path(output)
        output.mkdir(parents=True,exist_ok=True)
        data=pickle.dumps(result,protocol=5)
        (output/'basis.pkl').write_bytes(data)
        record['basis_sha256']=hashlib.sha256(data).hexdigest()
        (output/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    return result,record


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--sector',required=True)
    parser.add_argument('--rank',type=int,required=True)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--seconds',type=int,default=120)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('projected basis wall cap')))
    signal.alarm(args.seconds)
    build(tuple(map(int,args.sector.split(','))),args.rank,args.prime,args.output)
