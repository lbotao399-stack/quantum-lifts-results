"""Match Dot's R-007 source to the native finite-field correction coordinates."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import pickle
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT/'algorithms/su3_independent_prime_20261118'
PACKET = '1f2d7e10b7e25ad5cd2f0213f577566f636e92343f25f32334a1eada58e38fa4'
sys.path.insert(0, str(PACKAGE/'runtime/work/sp3-native.nosync'))
import sp_environment  # Registers the vendored full-jet trace runtime.
from su_trace_adapter_v1 import canonical


def compare(path):
    assert sha256(path.read_bytes()).hexdigest() == PACKET
    with ZipFile(path) as packet:
        bare = packet.read('q1/input/dmk_o330_hw.json')
        assert bare == packet.read('q0-certificate/source/dmk_o330_hw.json')
        image = packet.read('q1/input/Q1.json.gz')
    public = ROOT/'algorithms/dual_singlets_20261118/taylor_hooks/su3_bmn_sources/dmk_o330_hw.json.gz'
    terms = json.loads(bare)['terms']
    assert json.loads(gzip.decompress(public.read_bytes()))['terms'] == terms
    # Translate the public BMN fermions before canonicalizing the graded
    # cyclic words. Only the fermion at public index 4 has a negative sign.
    letters = [(0,0,(i,)) for i in range(3)] + [(0,0,(1,2)), (0,0,(0,2)), (0,0,(0,1))]
    translated = defaultdict(Fraction)
    for coefficient, traces in terms:
        word = tuple(tuple(letters[i] for i in trace) for trace in traces)
        fermions_at_four = sum(i == 4 for trace in traces for i in trace)
        translated[word] += (-1)**fermions_at_four * Fraction(coefficient)
    expression = canonical(translated)
    assert len(expression) == len(terms) == 6139
    # Recreate the exact bytes consumed by the production solver. No pickle is loaded.
    encoded = pickle.dumps({'expression':expression}, protocol=4)
    manifest = json.loads((PACKAGE/'quantum_source/result.json').read_bytes())
    assert sha256(encoded).hexdigest() == manifest['source_sha256']
    assert image == (PACKAGE/'quantum_source/Q1.json.gz').read_bytes()
    assert sha256(image).hexdigest() == manifest['files']['Q1.json.gz']
    assert manifest['image_terms']['Q1'] == 213104
    return dict(public_terms=len(terms), public_sha256=sha256(bare).hexdigest(),
                native_source_bytes=len(encoded), native_source_sha256=sha256(encoded).hexdigest(),
                q1_terms=manifest['image_terms']['Q1'], compressed_q1_bytes=len(image),
                compressed_q1_sha256=sha256(image).hexdigest(), packet_sha256=PACKET)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    print(json.dumps(compare(args.packet), indent=2))
