"""Evaluate the saved rational -Q1 x in the same target frames as Q0."""
import sp_environment
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import pickle
import gzip
import time
import numpy as np
from sp_full_profiles import profile
from native_su_parallel_v1 import evaluate
from su3_correction_map_v2 import PRIME, target_frame, save, sha

HERE = Path(__file__).resolve().parent


def run(index, source, output, chunk=2048, workers=8):
    started = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plan = json.loads((HERE/'su3_correction_profiles_v1.json').read_text())
    target = plan['targets'][index]
    selected = tuple((tuple(k), n) for k, n in target['profile'])
    source = Path(source)
    declared = json.loads((source.parent/'result.json').read_text())
    assert source.name in ('Q1.pkl','Q1.json.gz')
    assert sha(source) == declared['files'][source.name]
    if source.name == 'Q1.pkl':
        expression = pickle.loads(source.read_bytes())['expression']
    else:
        with gzip.open(source, 'rt') as stream:
            terms = json.load(stream)['terms']
        expression = {tuple(tuple((letter[0],letter[1],tuple(letter[2])) for letter in trace) for trace in word):Fraction(c) for c,word in terms}
    terms = [(w, (-Fraction(c).numerator*pow(Fraction(c).denominator,-1,PRIME)) % PRIME)
             for w,c in sorted(expression.items()) if profile(tuple(L for t in w for L in t)) == selected]
    frame, description = target_frame(target)
    manifest = dict(target=target, frame=description, input_terms=len(terms), all_input_terms=len(expression),
                    source_sha256=sha(source), driver_sha256=sha(__file__), chunk=chunk,
                    map_driver_sha256=sha(HERE/'su3_correction_map_v2.py'))
    marker = output/'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text()) == json.loads(json.dumps(manifest))
    else:
        save(marker, manifest)
    total = np.zeros(frame.ncols, dtype=np.int64)
    for start in range(0, len(terms), chunk):
        path = output/('chunk%06d.npy'%start)
        receipt = path.with_suffix('.json')
        if receipt.exists():
            record = json.loads(receipt.read_text())
            assert sha(path) == record['sha256']
            value = np.load(path).astype(np.int64)
        else:
            before = time.monotonic()
            value = evaluate([dict(terms[start:start+chunk])], frame, workers)[0]
            temporary = path.with_suffix('.tmp.npy')
            np.save(temporary, value.astype(np.uint32))
            temporary.replace(path)
            record = dict(start=start, stop=min(start+chunk,len(terms)), sha256=sha(path),
                          seconds=time.monotonic()-before)
            save(receipt, record)
        total = (total+value) % PRIME
        save(output/'progress.json', dict(target=index, completed_terms=min(start+chunk,len(terms)),
             required_terms=len(terms), seconds=time.monotonic()-started))
    path = output/'rhs.npy'
    np.save(path, total.astype(np.uint32))
    result = dict(status='complete_negative_Q1_evaluations', target=index, input_terms=len(terms),
                  columns=frame.ncols, nonzero_entries=int(np.count_nonzero(total)), sha256=sha(path),
                  seconds=time.monotonic()-started, source_sha256=sha(source), correction_solved=False,
                  scope='Right-hand side evaluations in the fixed coupled target frames; no exactness conclusion.')
    save(output/'result.json', result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--chunk', type=int, default=2048)
    parser.add_argument('--workers', type=int, default=8)
    run(**vars(parser.parse_args()))
