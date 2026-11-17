"""Replay all saved G2 determinants and selected exact symbolic Q0 probes.

This verifies the production minor determinants, but reconstructs only the
archived probe entries. It is not a reconstruction of every minor entry.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import time

import g2_bmn_complex_v1 as base
import numpy as np
from flint import nmod_mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    hashes = json.loads((args.evidence / 'SHA256SUMS.json').read_text())
    for name, digest in hashes.items():
        assert base.sha(args.evidence / name) == digest
    manifest = json.loads((args.evidence / 'frame_manifest.json').read_text())
    base.setup(tuple(manifest['weights']), tuple(range(1, 10)), manifest['prime'])
    assert json.loads(json.dumps(base.CONTEXT['dimensions'])) == manifest['dimensions']
    records = []
    for path in sorted(args.evidence.glob('s*.pkl')):
        payload = pickle.loads(path.read_bytes())
        spin, prime = payload['spin'], payload['prime']
        minor = payload['minor']
        rank = payload['result']['outgoing_rank_lower_bound']
        assert minor.shape == (rank, rank)
        determinant = int(nmod_mat(minor.tolist(), prime).det())
        assert determinant
        images = []
        for word in payload['probe_recipes']:
            expression = base.hw_extremal.project({word: 1}, (spin,) + tuple(payload['weights']))
            assert not base.hw_extremal.action(expression, 0, 1)
            assert not base.hw_extremal.action(expression, 1, 2)
            image = base.native.sp_census.base.action_expr(expression)
            assert not base.native.sp_census.base.action_expr(image)
            images.append({word: Fraction(value).numerator * pow(Fraction(value).denominator, -1, prime) % prime
                           for word, value in image.items() if value})
        coordinates, nonzero = 0, 0
        for probe in payload['probes']:
            frame = base.native.Frame(probe['degree'], probe['points'], 1, prime, k=probe['k'])
            frame.L = probe['letters']
            expressions = [{word: value for word, value in image.items()
                            if sum(3 <= letter <= 5 for trace in word[1] for letter in trace) == probe['degree']}
                           for image in images]
            actual = base.original_evaluate(expressions, frame)
            assert np.array_equal(actual, probe['values']), (spin, probe['degree'])
            coordinates += actual.size
            nonzero += int(np.count_nonzero(actual))
        records.append(dict(spin=spin, rank_lower_bound=rank, minor_determinant=determinant,
            independent_probe_recipes=len(images), coordinates=coordinates, nonzero_coordinates=nonzero))
        print(json.dumps(records[-1]), flush=True)
    ranks = {row['spin']: row['rank_lower_bound'] for row in records}
    ranks.update({0: 0, 10: 0})
    dimensions = {spin: sum(parts.values()) for spin, parts in base.CONTEXT['dimensions'].items()}
    excluded = []
    for spin in range(10):
        if spin in ranks and spin + 1 in ranks:
            upper = dimensions[spin] - ranks[spin] - ranks[spin + 1]
            assert upper >= 0
            if upper == 0:
                excluded.append(spin)
    assert excluded == [0, 1, 6, 7, 8, 9]
    result = dict(status='passed_all_saved_minors_exact_dimensions_and_independent_symbolic_probes',
        records=records, payload_hashes=len(hashes), excluded_spins=excluded,
        total_probe_coordinates=sum(row['coordinates'] for row in records),
        driver_sha256=base.sha(__file__), seconds=time.monotonic()-started,
        all_minor_entries_independently_reconstructed=False,
        scope='Every saved minor determinant and exact character dimension; selected rational projections, exact Q0 nilpotence and dense-matrix production probes. Full minor entries are retained production data.')
    (args.evidence / 'replay.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key:value for key,value in result.items() if key != 'records'}), flush=True)


if __name__ == '__main__':
    main()
