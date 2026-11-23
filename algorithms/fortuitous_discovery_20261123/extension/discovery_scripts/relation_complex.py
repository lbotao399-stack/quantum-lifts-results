"""Certify the outgoing rank in the small finite trace-identity complex.

The previous rational witness supplies a 122-dimensional closed subspace:
121 independent exact boundaries and one rational nonboundary. If a 61
minor is nonzero in the 183-dimensional target relation space, its
cohomology therefore has dimension exactly one over Q.
"""
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import signal
import sys
import time

HERE=Path(__file__).resolve().parent
JETS=HERE.parent/'fortuity_discovery_20261123/jets'
sys.path.insert(0,str(JETS))
import pfaffian_relations as pf
from flint import nmod_mat


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rank',type=int,default=4)
    ap.add_argument('--seconds',type=int,default=300)
    ap.add_argument('--prime',type=int,default=1000003)
    args=ap.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('relation-complex cap')))
    signal.alarm(args.seconds)
    start=time.monotonic()
    source_path=HERE/f'complete_target_sp{args.rank}.pkl'
    rows,recipes=pickle.loads(source_path.read_bytes())
    assert len(rows)==183
    cache=HERE/f'outgoing_sp{args.rank}'
    cache.mkdir(exist_ok=True)
    echelon=pf.Echelon(args.prime)
    selected=[]
    selected_indices=[]
    for i,row in enumerate(rows):
        path=cache/f'row_{i}.pkl'
        digest=hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest()
        if path.exists():
            old,image=pickle.loads(path.read_bytes())
            assert old==digest
        else:
            image=pf.m.q(row)
            path.write_bytes(pickle.dumps((digest,image),protocol=4))
        if echelon.add(image):
            selected.append(image)
            selected_indices.append(i)
        if i%10==0 or i>=181:
            print(json.dumps(dict(row=i,rank=len(selected),elapsed_seconds=time.monotonic()-start)),flush=True)
        assert len(selected)<=61
        if len(selected)==61:
            break
    index_to_word={i:w for w,i in echelon.words.items()}
    columns=[index_to_word[i] for i in sorted(echelon.rows)]
    minor=nmod_mat([[r.get(w,0)%args.prime for w in columns] for r in selected],args.prime)
    determinant=int(minor.det())
    assert determinant
    result=dict(rank=args.rank,relation_dimension=len(rows),prime=args.prime,
                outgoing_rank_lower_bound=len(selected),
                selected_relation_indices=selected_indices,
                determinant_mod_prime=determinant,
                pivot_words=[pf.m.encode(w) for w in columns],
                complete_target_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                established_rational_closed_subspace_dimension=122,
                established_rational_boundary_dimension=121,
                established_rational_nonboundary_count=1,
                elapsed_seconds=time.monotonic()-start)
    if len(selected)==61:
        result.update(outgoing_rational_rank=61,relation_cohomology_dimension=183-121-61,
                      exact_fortuitous_quotient_dimension=1,
                      full_finite_cohomology_dimension=1,
                      uses_previous_stable_centre_acyclicity=True)
    (HERE/f'relation_complex_sp{args.rank}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('pivot_words','selected_relation_indices')}),flush=True)


if __name__=='__main__':main()
