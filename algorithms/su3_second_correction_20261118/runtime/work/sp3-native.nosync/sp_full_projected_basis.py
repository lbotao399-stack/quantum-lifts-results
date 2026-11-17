"""Complete finite-rank full-letter cochains from projected trace words.

Every retained row passes all three raising equations. Evaluation rank must
reach the independently computed integer dimension in every fermion degree.
This does not calculate a universal formal-trace cohomology space.
"""
import full_extremal as projection
import full_hw_so as full
from sp_full_dimensions import highest_dimensions
from sp_projected_basis import independent_rows
from fortuity.full.evalg import Frame,nferm
import argparse
import contextlib
import gc
import hashlib
import io
import json
import pathlib
import pickle
import random
import signal
import time
from math import ceil,comb
import numpy as np
from scipy.sparse import csr_matrix


def evaluate(expressions,frame):
    words=sorted({w for e in expressions for w in e})
    if not words:
        return np.zeros((len(expressions),frame.ncols),dtype=np.int64)
    index={w:i for i,w in enumerate(words)}
    rr,cc,vv=[],[],[]
    for i,expr in enumerate(expressions):
        for word,c in expr.items():
            rr.append(i);cc.append(index[word]);vv.append(int(c)%frame.p)
    assert max(map(len,expressions),default=0)*(frame.p-1)**2<2**63
    matrix=csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(expressions),len(words)))
    return np.asarray(matrix@frame.eval_words(words,dtype=np.int64))%frame.p


def build(sector,letters,derivatives,rank,prime,output,batch_size=16):
    start=time.monotonic()
    output=pathlib.Path(output)
    output.mkdir(parents=True,exist_ok=True)
    manifest=dict(sector=list(sector),letters=letters,derivatives=derivatives,rank=rank,prime=prime,
                  sources={pathlib.Path(p).name:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
                           for p in [__file__,projection.__file__,full.__file__,
                                     pathlib.Path(__file__).with_name('sp_full_dimensions.py'),
                                     pathlib.Path(__file__).with_name('sp_dimensions.py')]})
    marker=output/'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text())==manifest
    else:
        marker.write_text(json.dumps(manifest,indent=2)+'\n')
    if (output/'result.json').exists():
        result=json.loads((output/'result.json').read_text())
        data=(output/'basis.pkl').read_bytes()
        assert hashlib.sha256(data).hexdigest()==result['basis_sha256']
        return pickle.loads(data),result
    exact=highest_dimensions(sector,letters,rank)
    with contextlib.redirect_stdout(io.StringIO()):
        shapes=full.enumerate_shapes(sector,letters,derivatives)
    words=sorted(set().union(*shapes.values())) if shapes else []
    del shapes
    print(json.dumps(dict(stage='full_weight_words',sector=sector,letters=letters,derivatives=derivatives,
                          weight_words=len(words),exact_dimensions=exact,seconds=time.monotonic()-start)),flush=True)
    random.Random(59623).shuffle(words)
    frames,chosen,pending,values,seen={},{},{},{},{}
    for degree,dimension in exact.items():
        alphabet=sorted({L for w in words if nferm(w)==degree for trace in w for L in trace})
        k=degree+2 if degree else 0
        frames[degree]=Frame(rank,ceil((dimension+16)/comb(k,degree)),k,[degree],alphabet,
                             seed=257371+7919*degree,p=prime,group='Sp')
        chosen[degree],pending[degree],seen[degree]=[],[],set()
        values[degree]=np.zeros((0,frames[degree].ncols),dtype=np.int64)
    examined=distinct=0

    def flush(degree):
        if not pending[degree]:
            return
        new=evaluate(pending[degree],frames[degree])
        joint=np.concatenate([values[degree],new],axis=0)
        ids=independent_rows(joint,prime)
        expressions=chosen[degree]+pending[degree]
        chosen[degree]=[expressions[i] for i in ids]
        values[degree]=joint[ids]
        pending[degree]=[]
        assert len(ids)<=exact[degree]

    last=start
    for word in words:
        degree=nferm(word)
        if degree not in exact or len(chosen[degree])==exact[degree]:
            continue
        examined+=1
        expr=projection.project({word:1},sector,prime)
        if not expr:
            continue
        inverse=pow(expr[min(expr)],-1,prime)
        expr={w:c*inverse%prime for w,c in expr.items()}
        fingerprint=tuple(sorted(expr.items()))
        if fingerprint in seen[degree]:
            continue
        seen[degree].add(fingerprint)
        distinct+=1
        for raising in (1,2,3):
            assert not full.raise_expr(expr,raising,prime)
        pending[degree].append(expr)
        if len(pending[degree])>=min(batch_size,exact[degree]-len(chosen[degree])+8):
            flush(degree)
        if time.monotonic()-last>15:
            print(json.dumps(dict(stage='full_projected_basis',dimensions={d:len(v) for d,v in chosen.items()},
                                  exact_dimensions=exact,examined=examined,seconds=time.monotonic()-start)),flush=True)
            projection.word_action.cache_clear()
            gc.collect()
            last=time.monotonic()
        if all(len(chosen[d])==dimension for d,dimension in exact.items()):
            break
    for degree in exact:
        flush(degree)
        assert len(chosen[degree])==exact[degree],('Incomplete full basis',degree,len(chosen[degree]),exact[degree])
    basis=[expr for degree in sorted(chosen) for expr in chosen[degree]]
    data=pickle.dumps(basis,protocol=5)
    (output/'basis.pkl').write_bytes(data)
    result=dict(group='Sp(%d)'%rank,sector=sector,letters=letters,derivatives=derivatives,prime=prime,
                exact_dimensions=exact,dimension=len(basis),canonical_weight_words=len(words),
                examined=examined,distinct_projections=distinct,
                evaluation_columns={d:fr.ncols for d,fr in frames.items()},
                highest_weight_checks=3*distinct,basis_sha256=hashlib.sha256(data).hexdigest(),
                manifest_sha256=hashlib.sha256(marker.read_bytes()).hexdigest(),seconds=time.monotonic()-start,
                scope='Complete finite-rank full-letter cochain basis over the stated prime. No universal-closed subtraction or quantum obstruction was calculated.')
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return basis,result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--sector',required=True)
    parser.add_argument('--letters',type=int,required=True)
    parser.add_argument('--derivatives',type=int,required=True)
    parser.add_argument('--rank',type=int,required=True)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--seconds',type=int,default=300)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('full projected basis cap')))
    signal.alarm(args.seconds)
    build(tuple(map(int,args.sector.split(','))),args.letters,args.derivatives,args.rank,args.prime,args.output)
