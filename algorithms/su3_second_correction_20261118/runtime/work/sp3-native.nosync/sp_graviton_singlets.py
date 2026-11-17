"""A specified closed product family; no completeness claim for the ring."""
import sp_environment
import sp_census as census
from fortuity.bmn.words import graviton_generators
from protected import graph_expr,expr_sec,monomials
from hw_products import Raising
from projected_checkpoint import atomic,json_bytes
from collections import defaultdict
import argparse
import hashlib
import json
import os
import pathlib
import pickle
import signal
import time
import numpy as np
from scipy.sparse import csr_matrix


def prepare(sector,rank,prime,output):
    started=time.monotonic()
    output=pathlib.Path(output)
    output.mkdir(parents=True,exist_ok=True)
    census.configure(rank,prime)
    generators=[]
    for name,expression in graviton_generators(4,1,qmax=sector):
        expression=graph_expr(expression)
        if not expression:
            continue
        assert not census.base.action_expr(expression),('Generator not formally closed',name)
        degrees={sum(census.base.count(w)[3:6])%2 for w in expression}
        assert len(degrees)==1
        generators.append(dict(name=name,expr=expression,sector=expr_sec(expression),
                               odd=degrees.pop(),reflection=0,family='trace'))
    products=monomials(generators,sector,0)
    print(json.dumps(dict(stage='graviton_products',generators=len(generators),weight_products=len(products))),flush=True)
    raising=Raising(generators,prime)
    images=[{(i,m):c for i in (1,2) for m,c in raising.product(product,i).items()} for product in products]
    kernel,raising_rank=census.base.formal_closed_kernel(images,prime)
    checks=0
    for row in kernel:
        check=defaultdict(int)
        for index in np.flatnonzero(row):
            for word,coefficient in images[index].items():
                check[word]+=int(row[index])*coefficient
        assert not any(value%prime for value in check.values())
        checks+=2
    record=dict(group='Sp(%d)'%rank,sector=sector,prime=prime,
                specified_generator_degrees=[2,4],generators=len(generators),
                formally_closed_generator_checks=len(generators),weight_products=len(products),
                raising_rank=raising_rank,formal_highest_weight_products=len(kernel),
                raising_checks=checks,seconds=time.monotonic()-started,
                scope='Only a known closed product span. No full-ring or cohomology completeness is asserted.')
    data=pickle.dumps(dict(generators=generators,products=products,kernel=kernel),protocol=5)
    atomic(output/'family.pkl',data)
    record['family_sha256']=hashlib.sha256(data).hexdigest()
    record['source_sha256']=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
    atomic(output/'family.json',json_bytes(record))
    print(json.dumps(record),flush=True)
    return record


def product_values(generators,products,frame,degree):
    """Keep the mixed fermion degrees of w until the full product is formed."""
    indices=sorted({index for product in products for index in product})
    necklaces=set()
    for index in indices:
        for paths,traces in generators[index]['expr']:
            assert not paths and len(traces)==1
            necklaces.add(traces[0])
    traced=frame.necklace_values(necklaces)
    values={}
    for index in indices:
        parts={}
        for (_,traces),coefficient in generators[index]['expr'].items():
            raw,current=traced[traces[0]]
            if current>degree:
                continue
            part=raw*(coefficient%frame.p)%frame.p
            parts[current]=(parts[current]+part)%frame.p if current in parts else part
        values[index]=parts
    result=np.zeros((len(products),frame.ncols),dtype=np.int64)
    for index,product in enumerate(products):
        parts={0:np.ones((frame.P,1),dtype=np.int64)}
        for factor in product:
            following={}
            for left_degree,left in parts.items():
                for right_degree,right in values[factor].items():
                    if left_degree+right_degree>degree:
                        continue
                    raw,total=frame.scal_mul(left,left_degree,right,right_degree)
                    following[total]=(following[total]+raw)%frame.p if total in following else raw
            parts=following
        if degree in parts:
            result[index]=parts[degree].reshape(-1)
    return result


def evaluate_family(family_folder,basis_folder,incoming,output):
    import projected_differential as differential
    from flint import nmod_mat
    started=time.monotonic()
    family_folder=pathlib.Path(family_folder)
    record=json.loads((family_folder/'family.json').read_text())
    data=(family_folder/'family.pkl').read_bytes()
    assert hashlib.sha256(data).hexdigest()==record['family_sha256']
    family=pickle.loads(data)
    _,basis_record,manifest,dimensions=differential.load_basis(basis_folder)
    assert manifest['sector']==record['sector'] and manifest['prime']==record['prime']
    prime=record['prime']
    differential.initialize(manifest['rank'],prime,dimensions)
    coefficients=csr_matrix(family['kernel'])
    assert int(max(np.diff(coefficients.indptr),default=0))*(prime-1)**2<2**63
    blocks=[]
    for degree in sorted(dimensions):
        raw=product_values(family['generators'],family['products'],differential._frames[degree],degree)
        blocks.append(np.asarray(coefficients@raw)%prime)
        print(json.dumps(dict(stage='graviton_values',degree=degree,seconds=time.monotonic()-started)),flush=True)
    values=np.concatenate(blocks,axis=1)
    incoming=pathlib.Path(incoming)
    previous=json.loads((incoming/'manifest.json').read_text())
    previous_result=json.loads((incoming/'result.json').read_text())
    assert previous['target_basis_sha256']==basis_record['basis_sha256']
    data=(incoming/'matrix.npy').read_bytes()
    assert hashlib.sha256(data).hexdigest()==previous_result['matrix_sha256']
    boundaries=np.load(incoming/'matrix.npy',allow_pickle=False)
    incoming_rank=int(nmod_mat(boundaries.tolist(),prime).rank())
    combined_rank=int(nmod_mat(np.concatenate([boundaries,values],axis=0).tolist(),prime).rank())
    assert incoming_rank==previous_result['rank']
    result=dict(group=record['group'],sector=record['sector'],prime=prime,
                dimension=basis_record['dimension'],incoming_rank=incoming_rank,
                known_closed_rank=combined_rank,known_cohomology_lower_bound=combined_rank-incoming_rank,
                source_family_sha256=record['family_sha256'],seconds=time.monotonic()-started,
                scope='Exact rank of this closed product family modulo boundaries. Completeness requires the outgoing map.')
    atomic(pathlib.Path(output),json_bytes(result))
    np.save(pathlib.Path(output).with_suffix('.npy'),values,allow_pickle=False)
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    print('FAMILY_PID',os.getpid(),flush=True)
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    parser.add_argument('--seconds',type=int,default=180)
    args=parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('graviton family pilot cap')))
    signal.alarm(args.seconds)
    prepare((5,5,5,5),2,1000003,args.output)
