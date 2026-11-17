"""Checkpointed finite-rank cochains, with independent dimension saturation.

Each committed snapshot includes its evaluation row space and selected trace
expressions. It is hashed before an atomic pointer is replaced. Incomplete
ranks are not cohomology results. Both BMN and full-letter adapters are used.
"""
import sp_environment
import argparse
import contextlib
import gc
import hashlib
import io
import json
import os
import pathlib
import pickle
import random
import signal
import time
from math import ceil, comb
import numpy as np
from flint import nmod_mat
from sp_projected_basis import independent_rows


class CheckpointStop(Exception):
    pass


def atomic(path, data):
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('wb') as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def json_bytes(data):
    return (json.dumps(data, indent=2)+'\n').encode()


def extend_rows(previous, pivots, values, prime):
    """Keep the pivot block equal to the identity; select original new rows."""
    assert max(previous.shape[0], values.shape[0])*(prime-1)**2 < 2**63
    reduced = values.copy()
    if pivots:
        reduced = (reduced-reduced[:,pivots]@previous) % prime
    ids = independent_rows(reduced, prime)
    if not ids:
        return previous, pivots, ids
    matrix, rank = nmod_mat(reduced[ids].tolist(), prime).rref()
    assert int(rank) == len(ids)
    new = np.array(matrix.tolist(), dtype=np.int64)
    new_pivots = [int(np.flatnonzero(row)[0]) for row in new]
    assert not set(pivots).intersection(new_pivots)
    if len(previous):
        previous = (previous-previous[:,new_pivots]@new) % prime
    result = np.concatenate((previous,new),axis=0)
    all_pivots = pivots+new_pivots
    assert np.array_equal(result[:,all_pivots], np.eye(len(result),dtype=np.int64))
    return result, all_pivots, ids


def build(mode, sector, rank, prime, output, letters=None, derivatives=None,
          batch_size=24, stop_after_batches=0):
    start = time.monotonic()
    output = pathlib.Path(output)
    output.mkdir(parents=True,exist_ok=True)
    if mode == 'bmn':
        import sp_census as census
        import hw_extremal as projection
        from sp_projected_basis import evaluate
        census.configure(rank,prime)
        exact = census.base.exact_hw(sector,0)
        enumerate_words = lambda: list(census.base.sector_basis(sector,0))
        degree_of = lambda word: sum(census.base.count(word)[3:6])
        def make_frame(degree, dimension, words):
            k = degree+2 if degree else 0
            return census.base.Frame(degree,ceil((dimension+16)/comb(k,degree)),
                                     125131+7919*degree,prime,k=k)
        def check(expr):
            assert not projection.action(expr,0,1,prime)
            assert not projection.action(expr,1,2,prime)
        additional = [pathlib.Path(census.__file__),
                      pathlib.Path(census.__file__).with_name('sp_dimensions.py')]
        checks_per_projection = 2
    elif mode == 'full':
        import full_extremal as projection
        import full_hw_so as full
        from sp_full_dimensions import highest_dimensions
        from sp_full_projected_basis import evaluate
        from fortuity.full.evalg import Frame,nferm
        exact = highest_dimensions(sector,letters,rank)
        def enumerate_words():
            with contextlib.redirect_stdout(io.StringIO()):
                shapes = full.enumerate_shapes(sector,letters,derivatives)
            return sorted(set().union(*shapes.values())) if shapes else []
        degree_of = nferm
        def make_frame(degree,dimension,words):
            alphabet = sorted({letter for word in words if nferm(word)==degree
                               for trace in word for letter in trace})
            k = degree+2 if degree else 0
            return Frame(rank,ceil((dimension+16)/comb(k,degree)),k,[degree],alphabet,
                         seed=257371+7919*degree,p=prime,group='Sp')
        def check(expr):
            for raising in (1,2,3):
                assert not full.raise_expr(expr,raising,prime)
        additional = [pathlib.Path(full.__file__),
                      pathlib.Path(__file__).with_name('sp_full_dimensions.py'),
                      pathlib.Path(__file__).with_name('sp_dimensions.py')]
        checks_per_projection = 3
    else:
        raise ValueError(mode)
    sources = [pathlib.Path(__file__),pathlib.Path(projection.__file__),
               pathlib.Path(sp_environment.__file__),pathlib.Path(__file__).with_name('sp_projected_basis.py')]+additional
    manifest = dict(mode=mode,sector=list(sector),rank=rank,prime=prime,
                    letters=letters,derivatives=derivatives,batch_size=batch_size,
                    exact_dimensions={str(k):v for k,v in exact.items()},
                    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    marker = output/'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text()) == manifest, 'Checkpoint manifest mismatch'
    else:
        atomic(marker,json_bytes(manifest))
    manifest_hash = hashlib.sha256(marker.read_bytes()).hexdigest()
    if (output/'result.json').exists():
        result = json.loads((output/'result.json').read_text())
        data = (output/'basis.pkl').read_bytes()
        assert result['manifest_sha256'] == manifest_hash
        assert hashlib.sha256(data).hexdigest() == result['basis_sha256']
        return pickle.loads(data),result
    if (output/'words.json').exists():
        record = json.loads((output/'words.json').read_text())
        data = (output/'words.pkl').read_bytes()
        assert record['manifest_sha256'] == manifest_hash
        assert hashlib.sha256(data).hexdigest() == record['sha256']
        words = pickle.loads(data)
    else:
        words = enumerate_words()
        random.Random(59623 if mode=='full' else 38711).shuffle(words)
        data = pickle.dumps(words,protocol=5)
        atomic(output/'words.pkl',data)
        atomic(output/'words.json',json_bytes(dict(sha256=hashlib.sha256(data).hexdigest(),
                                                  manifest_sha256=manifest_hash,count=len(words))))
    print(json.dumps(dict(stage='checkpoint_words',mode=mode,sector=sector,
                          words=len(words),exact_dimensions=exact)),flush=True)
    frames = {d:make_frame(d,dimension,words) for d,dimension in exact.items()}
    pointer = output/'checkpoint.json'
    if pointer.exists():
        receipt = json.loads(pointer.read_text())
        assert receipt['manifest_sha256'] == manifest_hash
        data = (output/receipt['file']).read_bytes()
        assert hashlib.sha256(data).hexdigest() == receipt['sha256'], 'Checkpoint hash mismatch'
        state = pickle.loads(data)
        reused_batches = state['batches']
    else:
        state = dict(cursor=0,examined=0,distinct=0,batches=0,
                     chosen={d:[] for d in exact},seen={d:set() for d in exact},
                     pivots={d:[] for d in exact},
                     rows={d:np.zeros((0,frames[d].ncols),dtype=np.int64) for d in exact})
        reused_batches = 0
    pending = {d:[] for d in exact}
    def commit():
        for degree in exact:
            if not pending[degree]:
                continue
            new = evaluate(pending[degree],frames[degree])
            rows,pivots,ids = extend_rows(state['rows'][degree],state['pivots'][degree],new,prime)
            state['rows'][degree],state['pivots'][degree] = rows,pivots
            state['chosen'][degree].extend(pending[degree][i] for i in ids)
            pending[degree] = []
            assert len(pivots) == len(state['chosen'][degree]) <= exact[degree]
        state['batches'] += 1
        data = pickle.dumps(state,protocol=5)
        filename = 'state_%06d.pkl'%state['batches']
        atomic(output/filename,data)
        old_receipt = json.loads(pointer.read_text()) if pointer.exists() else None
        receipt = dict(file=filename,sha256=hashlib.sha256(data).hexdigest(),
                       manifest_sha256=manifest_hash,cursor=state['cursor'],
                       batches=state['batches'],examined=state['examined'],
                       dimensions={d:len(v) for d,v in state['chosen'].items()},
                       exact_dimensions=exact,seconds_this_run=time.monotonic()-start)
        atomic(pointer,json_bytes(receipt))
        if old_receipt and old_receipt['file'] != filename:
            (output/old_receipt['file']).unlink()
        print(json.dumps(dict(stage='projected_checkpoint',**receipt)),flush=True)
        projection.word_action.cache_clear()
        gc.collect()
        if stop_after_batches and state['batches']-reused_batches >= stop_after_batches:
            raise CheckpointStop()
    for position in range(state['cursor'],len(words)):
        word = words[position]
        degree = degree_of(word)
        if degree in exact and len(state['chosen'][degree]) < exact[degree]:
            state['examined'] += 1
            expr = projection.project({word:1},sector,prime)
            if expr:
                inverse = pow(expr[min(expr)],-1,prime)
                expr = {w:c*inverse%prime for w,c in expr.items()}
                fingerprint = hashlib.sha256(pickle.dumps(tuple(sorted(expr.items())),protocol=5)).digest()
                if fingerprint not in state['seen'][degree]:
                    check(expr)
                    state['seen'][degree].add(fingerprint)
                    state['distinct'] += 1
                    pending[degree].append(expr)
        state['cursor'] = position+1
        if sum(map(len,pending.values())) >= batch_size:
            commit()
        if all(len(state['chosen'][d]) == dimension for d,dimension in exact.items()):
            break
    if any(pending.values()):
        commit()
    assert all(len(state['chosen'][d]) == dimension for d,dimension in exact.items()), 'Incomplete projection rank'
    basis = [expr for d in sorted(exact) for expr in state['chosen'][d]]
    data = pickle.dumps(basis,protocol=5)
    atomic(output/'basis.pkl',data)
    result = dict(group='Sp(%d)'%rank,mode=mode,sector=sector,prime=prime,
                  letters=letters,derivatives=derivatives,dimension=len(basis),exact_dimensions=exact,
                  canonical_weight_words=len(words),examined=state['examined'],
                  distinct_projections=state['distinct'],batches=state['batches'],reused_batches=reused_batches,
                  highest_weight_checks=checks_per_projection*state['distinct'],
                  basis_sha256=hashlib.sha256(data).hexdigest(),manifest_sha256=manifest_hash,
                  seconds_this_run=time.monotonic()-start,
                  scope='Complete finite-rank cochain basis over the stated field. '
                        'No universal-closed subtraction or quantum obstruction is asserted.')
    atomic(output/'result.json',json_bytes(result))
    print(json.dumps(result),flush=True)
    return basis,result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode',choices=('bmn','full'),required=True)
    parser.add_argument('--sector',required=True)
    parser.add_argument('--rank',type=int,required=True)
    parser.add_argument('--letters',type=int)
    parser.add_argument('--derivatives',type=int)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--seconds',type=int,default=600)
    parser.add_argument('--batch-size',type=int,default=24)
    parser.add_argument('--stop-after-batches',type=int,default=0)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('projection time cap')))
    signal.alarm(args.seconds)
    try:
        build(args.mode,tuple(map(int,args.sector.split(','))),args.rank,args.prime,args.output,
              args.letters,args.derivatives,args.batch_size,args.stop_after_batches)
    except (CheckpointStop,TimeoutError) as error:
        print(json.dumps(dict(status=type(error).__name__,checkpoint=str(pathlib.Path(args.output)/'checkpoint.json'))),flush=True)
