#!/usr/bin/env python3
"""Portable entry point for the exact classical fortuity certificates."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import importlib
import json
import lzma
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
UPSTREAM=HERE/'upstream/algorithms'


def environment():
    os.environ['N4_WORKSPACE_ROOT']=str(HERE)
    os.environ['N4_DUAL_SOURCE']=str(UPSTREAM/'dual_singlets_20261118')
    os.environ['N4_TRACE_RUNTIME']=str(UPSTREAM/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync')
    os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    sys.dont_write_bytecode=True
    sys.path.insert(0,os.environ['N4_DUAL_SOURCE'])
    sys.path.insert(0,os.environ['N4_TRACE_RUNTIME'])


def hydrate(rank):
    folder=HERE/f'symplectic/sp{rank}'
    for path in folder.glob('*.pkl.xz'):
        target=path.with_suffix('')
        if not target.exists():
            target.write_bytes(lzma.decompress(path.read_bytes()))
    return folder


def symplectic(rank):
    start=time.monotonic();folder=hydrate(rank);sys.path.insert(0,str(folder))
    import mixed_incoming as m
    receipt=json.loads((folder/f'sp{rank}_discovery_certificate.json').read_text())
    data=json.loads(gzip.decompress((folder/f'sp{rank}_new_fortuitous_operator.json.gz').read_bytes()))
    polynomial={m.decode(w):Fraction(c) for c,w in data['terms']}
    denominator=math.lcm(*(v.denominator for v in polynomial.values()))
    integral={w:int(v*denominator) for w,v in polynomial.items()}
    variation=m.q(integral)
    assert not m.q(variation)
    source_length=2*rank+3
    source,_=pickle.loads((folder/f'exact_relations_sp{rank}_L{source_length}.pkl').read_bytes())
    target_name='target_snapshot_181.pkl' if rank==4 else 'target_sp5_snapshot_181.pkl'
    target,_=pickle.loads((folder/target_name).read_bytes())
    assert len(source)==188 and len(target)==181
    coefficients=list(map(Fraction,receipt['defect_identity_coefficients']))
    common=math.lcm(denominator,*(v.denominator for v in coefficients))
    combination=Counter()
    for c,row in zip(coefficients,target):
        for w,v in row.items():combination[w]+=int(c*common)*v
    assert {w:c for w,c in combination.items() if c}=={w:c*(common//denominator) for w,c in variation.items()}
    for raising in (1,2,3):
        raised=Counter()
        for w,c in integral.items():
            for v,d in m.raise_word(w,raising).items():raised[v]+=c*d
        assert not any(raised.values())
    dual_record=json.loads((folder/'rational_trace_dual.json').read_text()) if rank==4 else receipt
    dual=[(m.decode(w),Fraction(c)) for c,w in dual_record['dual_word_coefficients']]
    for i,row in enumerate(source):
        image=m.q(row)
        assert sum(c*image.get(w,0) for w,c in dual)==0,('boundary',i)
    pairing=sum(c*variation.get(w,0) for w,c in dual)/denominator
    assert pairing==1
    result={'passed':True,'group':f'Sp({rank})','native_weight':[3,3,2*rank-1,0,0],
        'length':source_length,'j':4*rank+16,'operator_terms':len(integral),
        'least_trace_coefficient_denominator':denominator,'variation_terms':len(variation),
        'exact_closure':True,'native_highest_weights':True,'source_Q0_images_recomputed':188,
        'dual_field':'Q','dual_nonzero_trace_coefficients':len(dual),'dual_boundary_pairings':'all zero',
        'dual_defect_pairing':str(pairing),'seconds':time.monotonic()-start}
    out=HERE/'replay_outputs';out.mkdir(exist_ok=True)
    (out/f'sp{rank}_witness.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


def identities(rank,seconds):
    folder=hydrate(rank);sys.path.insert(0,str(folder))
    import pfaffian_relations
    pfaffian_relations.identity=lru_cache(maxsize=768)(pfaffian_relations.identity)
    import path_relation_completion
    path_relation_completion.expand_paths=lru_cache(maxsize=768)(path_relation_completion.expand_paths)
    out=HERE/'replay_outputs';out.mkdir(exist_ok=True)
    if rank==4:
        import replay_identities
        replay_identities.main(out/'sp4_identities.json')
    else:
        import replay_sp5_identities
        replay_sp5_identities.main(out/'sp5_identities.json',out/'sp5_identity_checkpoint.json',seconds)


def dependencies():
    import mixed_incoming
    from trace_descent import dimensions
    assert sum(dimensions((3,3,7,0,0),11,8,'Sp').values())==11439
    for record in json.loads((HERE/'upstream/PROVENANCE.json').read_text())['files']:
        path=HERE/record['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
    foreign=[]
    for module in sys.modules.values():
        path=getattr(module,'__file__',None)
        if path and ('quantum_lifts_frontier_source' in path or 'so11_pin' in path):foreign.append(path)
    assert not foreign,foreign
    print(json.dumps({'passed':True,'pinned_source_hashes_verified':48,'original_workspace_imports':foreign},indent=2))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('part',choices=['scalar','theory','symplectic','identities','dependencies','all']);ap.add_argument('--rank',type=int,choices=[4,5],default=4);ap.add_argument('--seconds',type=int,default=300);a=ap.parse_args();environment()
    outputs=HERE/'replay_outputs';outputs.mkdir(exist_ok=True)
    if a.part=='scalar':
        subprocess.run([sys.executable,str(HERE/'scalar/verify_rank_birth.py'),'--rebuild-lower','--out',str(outputs/'scalar_rank_births.json')],check=True)
    elif a.part=='theory':
        subprocess.run([sys.executable,str(HERE/'structure/universal_auxiliary_identity.py'),'--out',str(outputs/'universal_auxiliary_identity.json')],check=True)
        subprocess.run([sys.executable,str(HERE/'structure/closed_dual.py'),'--verify','--out',str(outputs/'closed_dual.json')],check=True)
    elif a.part=='symplectic':symplectic(a.rank)
    elif a.part=='identities':identities(a.rank,a.seconds)
    elif a.part=='dependencies':dependencies()
    else:
        for part in ['dependencies','scalar','theory']:
            subprocess.run([sys.executable,str(HERE/'replay.py'),part],check=True)
        for rank in (4,5):
            subprocess.run([sys.executable,str(HERE/'replay.py'),'symplectic','--rank',str(rank)],check=True)


if __name__=='__main__':main()
