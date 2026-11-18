"""Exact evaluation system for rational graviton boundary primitives.

Full target evaluation rank is checked modulo the saved good prime against
the rational character dimension; this proves injectivity over Q. Vanishing
in these integer coordinates then certifies an actual polynomial identity.
"""
import sp_environment
import integer_sp_evaluation as exact
import hw_extremal as projection
from rational_j30_inputs import HERE, OUT, ROOT, PRIME, digest
from exact_trace_eval import rational_reconstruct
from projected_checkpoint import atomic,json_bytes
from fractions import Fraction
from math import ceil,comb,lcm
from flint import nmod_mat
from pathlib import Path
import argparse
import gc
import hashlib
import io
import json
import pickle
import time
import numpy as np

OUTPUT=HERE/'sp2_j30_integer_system'


def save_array(path,value):
    stream=io.BytesIO()
    np.save(stream,value,allow_pickle=False)
    atomic(path,stream.getvalue())


def source_inputs():
    record=json.loads((OUT/'primitive_screen.json').read_text())
    source=pickle.loads((OUT/'basis6.pkl').read_bytes())['basis']
    source=[source[i] for i in record['source_rows']]
    saved=pickle.loads((OUT/'gravitons.pkl').read_bytes())
    gravitons=[];denominators=[]
    for i in record['graviton_rows']:
        e=saved['expressions'][i]
        denominator=lcm(*(c.denominator for c in e.values()))
        gravitons.append({w:int(c*denominator) for w,c in e.items()})
        denominators.append(denominator)
    return source,gravitons,denominators,record


def target_inputs():
    destination=OUTPUT/'rational_basis5.pkl'
    if destination.exists():
        metadata=json.loads(destination.with_suffix('.json').read_text())
        assert digest(destination)==metadata['sha256']
        return pickle.loads(destination.read_bytes()),metadata
    basis=pickle.loads((ROOT/'basis5/basis.pkl').read_bytes())
    metadata=json.loads((ROOT/'basis5/result.json').read_text())
    assert digest(ROOT/'basis5/basis.pkl')==metadata['basis_sha256']
    reconstruction={c:rational_reconstruct(c,PRIME) for c in {c for e in basis for c in e.values()}}
    result=[];denominators=[]
    for i,e in enumerate(basis):
        den=lcm(*(reconstruction[c].denominator for c in e.values()))
        value={w:int(reconstruction[c]*den) for w,c in e.items()}
        assert all((value[w]-c*den)%PRIME==0 for w,c in e.items())
        assert not projection.action(value,0,1)
        assert not projection.action(value,1,2)
        result.append(value);denominators.append(den)
        if i%300==0:
            print(json.dumps(dict(stage='target_rational_reconstruction',rows=i+1)),flush=True)
    atomic(destination,pickle.dumps(result,protocol=5))
    metadata=dict(status='passed',sha256=digest(destination),rows=len(result),
        exact_raising_checks=2*len(result),maximum_denominator=max(denominators),
        modular_basis_sha256=digest(ROOT/'basis5/basis.pkl'),denominators=denominators)
    atomic(destination.with_suffix('.json'),json_bytes(metadata))
    projection.word_action.cache_clear()
    return result,metadata


def main(only_degree=None):
    start=time.monotonic()
    OUTPUT.mkdir(exist_ok=True)
    validated=json.loads((HERE/'integer_sp_evaluation_validation.json').read_text())
    assert validated['status']=='passed' and validated['source_sha256']['integer_sp_evaluation.py']==digest(Path(exact.__file__))
    source,gravitons,denominators,screen=source_inputs()
    basis,target_meta=target_inputs()
    manifest=dict(group='Sp(2)',sector=[5,5,5,5],prime=PRIME,
        source_indices=screen['source_rows'],graviton_rows=screen['graviton_rows'],
        graviton_denominators=denominators,rational_basis5_sha256=target_meta['sha256'],
        rational_basis6_sha256=digest(OUT/'basis6.pkl'),
        graviton_source_sha256=digest(OUT/'gravitons.pkl'),
        target_dimensions={'1':312,'3':1537,'5':1028},
        evaluator_sha256=digest(Path(exact.__file__)),seed_formula='817001+7919*degree')
    marker=OUTPUT/'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text())==manifest
    else:
        atomic(marker,json_bytes(manifest))
    marker_hash=digest(marker)
    controls=[]
    for degree,dimension in ((1,312),(3,1537),(5,1028)):
        if only_degree is not None and degree!=only_degree:
            continue
        receipt=OUTPUT/('degree%d.json'%degree)
        if receipt.exists():
            old=json.loads(receipt.read_text())
            assert old['manifest_sha256']==marker_hash
            assert all(digest(OUTPUT/name)==value for name,value in old['sha256'].items())
            controls.append(old)
            continue
        npts=ceil((dimension+16)/comb(degree+2,degree))
        frame=exact.IntegerSymplecticFrame(2,npts,degree+2,degree,817001+7919*degree)
        selected=[i for i,e in enumerate(basis) if sum(exact.fermions(t) for t in next(iter(e))[1])==degree]
        assert len(selected)==dimension
        def progress(stage):
            return lambda end,total,bound:print(json.dumps(dict(stage=stage,degree=degree,
                words=end,total=total,overflow_bound=bound,seconds=time.monotonic()-start)),flush=True)
        source_traces={t for e in source for w in e if abs(sum(exact.fermions(t) for t in w[1])-degree)==1 for t in w[1]}
        ordinary,images={},{}
        traces=sorted(source_traces)
        for offset in range(0,len(traces),4096):
            ordinary1,images1=exact.trace_values(traces[offset:offset+4096],frame,derivative=True)
            ordinary.update(ordinary1);images.update(images1)
            if offset%(4*4096)==0:
                print(json.dumps(dict(stage='integer_source_traces',degree=degree,
                    traces=min(offset+4096,len(traces)),total=len(traces),seconds=time.monotonic()-start)),flush=True)
        boundaries=exact.evaluate(source,frame,ordinary,images,derivative=True,progress=progress('integer_boundaries'))
        save_array(OUTPUT/('boundaries%d.npy'%degree),boundaries)
        del images,ordinary1,images1
        gc.collect()
        target=[basis[i] for i in selected]
        needed={t for e in target+gravitons for w in e if sum(exact.fermions(t) for t in w[1])==degree for t in w[1]}-ordinary.keys()
        traces=sorted(needed)
        for offset in range(0,len(traces),4096):
            ordinary1,_=exact.trace_values(traces[offset:offset+4096],frame)
            ordinary.update(ordinary1)
            if offset%(4*4096)==0:
                print(json.dumps(dict(stage='integer_target_traces',degree=degree,
                    traces=min(offset+4096,len(traces)),total=len(traces),seconds=time.monotonic()-start)),flush=True)
        target_values=exact.evaluate(target,frame,ordinary,progress=progress('integer_target_basis'))
        graviton_values=exact.evaluate(gravitons,frame,ordinary,progress=progress('integer_gravitons'))
        save_array(OUTPUT/('target%d.npy'%degree),target_values)
        save_array(OUTPUT/('gravitons%d.npy'%degree),graviton_values)
        rank=int(nmod_mat((target_values%PRIME).tolist(),PRIME).rank())
        assert rank==dimension,('target evaluation not injective',degree,rank,dimension)
        # Compare three actual source vectors to a separately expanded symbolic Q0.
        nonzero=[i for i in range(len(source)) if np.any(boundaries[i])]
        checks=0
        for i in (nonzero[0],nonzero[len(nonzero)//2],nonzero[-1]):
            expr=exact.census.base.action_expr(source[i])
            expected=exact.evaluate([expr],frame)[0]
            assert np.array_equal(expected,boundaries[i]),(degree,i)
            checks+=frame.ncols
        row=dict(status='passed',degree=degree,source_rows=len(source),target_dimension=dimension,
            target_rank_mod_prime=rank,columns=frame.ncols,independent_symbolic_entries=checks,
            target_basis_indices=selected,manifest_sha256=marker_hash,
            sha256={name:digest(OUTPUT/name) for name in ('boundaries%d.npy'%degree,
                'target%d.npy'%degree,'gravitons%d.npy'%degree)},seconds=time.monotonic()-start)
        atomic(receipt,json_bytes(row));controls.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('target_basis_indices','sha256')}),flush=True)
        del ordinary,target_values,graviton_values,boundaries,target,traces,ordinary1,needed
        exact.census.base.clear_caches();gc.collect()
    if all((OUTPUT/('degree%d.json'%d)).exists() for d in (1,3,5)):
        controls=[json.loads((OUTPUT/('degree%d.json'%d)).read_text()) for d in (1,3,5)]
        record=dict(status='passed',field='Z',group='Sp(2)',sector=[5,5,5,5],
            target_cochain_dimension=2877,total_evaluation_columns=sum(c['columns'] for c in controls),
            injective_over_Q=True,modular_prime=PRIME,graviton_denominators=denominators,
            blocks=controls,manifest_sha256=marker_hash,seconds=time.monotonic()-start,
            scope='Exact integer Q0 and graviton coordinates in an injective rational evaluation. Rational primitive solve is separate.')
        atomic(OUTPUT/'result.json',json_bytes(record))
        print(json.dumps({k:v for k,v in record.items() if k!='blocks'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--only-degree',type=int)
    args=parser.parse_args()
    main(args.only_degree)
