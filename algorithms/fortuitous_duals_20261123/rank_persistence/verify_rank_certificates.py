#!/usr/bin/env python3
"""Replay the saved finite-field minors and the exact rank-saturation arithmetic.

This does not regenerate native matrix entries. The separate production
drivers do that. No upstream imports, compiler, or pickle cache are required.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import time

import numpy as np
from flint import nmod_mat

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=HERE)
    ap.add_argument('--positive', type=Path, default=HERE.parent / 'positive_sector')
    ap.add_argument('--output', type=Path, default=HERE / 'rank_certificate_replay.json')
    args = ap.parse_args()
    start = time.monotonic()
    root = args.data
    lower = json.loads((root / 'lower_rank_scan.json').read_text())
    upper = json.loads((root / 'sp4_outgoing_rank.json').read_text())
    dims = json.loads((root / 'dimensions.json').read_text())
    assert lower['weight'] == upper['weight'] == [3, 3, 5, 0, 0]
    assert lower['length'] == upper['length'] == 9
    assert lower['prime'] == upper['prime'] == 1000003
    prime = lower['prime']
    checks = []
    table = []
    inputs = ['lower_rank_scan.json', 'sp4_outgoing_rank.json', 'dimensions.json', 'exact_rank_control.json']
    for group in lower['groups']:
        size = group['defining_dimension']
        name = 'sp%d_rank_minors.npz' % (size // 2)
        path = root / name
        assert sha(path) == group['minors_sha256']
        data = np.load(path)
        inputs.append(name)
        for key, certificate_key in [('incoming', 'incoming'), ('outgoing', 'outgoing'), ('augmented', 'incoming_plus_omega')]:
            certificate = group[certificate_key]
            minor = data[key]
            assert minor.shape == (certificate['rank'], certificate['rank'])
            determinant = int(nmod_mat(minor.tolist(), prime).det())
            assert determinant == certificate['determinant'] != 0
            checks.append({'rank_N': size // 2, 'map': key, 'size': len(minor), 'determinant': determinant})
        assert np.count_nonzero(data['omega_evaluation']) == group['omega_evaluation_nonzeros']
        dimensions = {int(k): v for k, v in group['character_dimensions'].items()}
        incoming, outgoing = group['incoming']['rank'], group['outgoing']['rank']
        assert dimensions[9] - incoming - outgoing == group['cohomology_upper_bound'] == 0
        table.append({'rank_N': size // 2, 'dimensions_C8_C9_C10': [dimensions[k] for k in (8, 9, 10)],
                      'incoming_rank': incoming, 'outgoing_rank': outgoing, 'H9_dimension': 0})

    path = root / 'sp4_outgoing_minor.npz'
    assert sha(path) == upper['minor_sha256']
    inputs.append(path.name)
    data = np.load(path)
    assert data['minor'].shape == (1198, 1198)
    assert data['source_rows'].tolist() == upper['selected_stable_outgoing_rows']
    determinant = int(nmod_mat(data['minor'].tolist(), prime).det())
    assert determinant == upper['outgoing_minor_determinant'] == 435941
    checks.append({'rank_N': 4, 'map': 'outgoing', 'size': 1198, 'determinant': determinant})

    dimension_map = {(row['defining_dimension'] // 2, row['length']): row['dimension'] for row in dims['records']}
    assert [dimension_map[(3, k)] for k in (8, 9, 10)] == [2559, 2689, 1281]
    assert [dimension_map[(4, k)] for k in (8, 9, 10)] == [2624, 2870, 1451]
    assert [dimension_map[(5, k)] for k in (8, 9, 10)] == [2624, 2870, 1455]
    assert upper['stable_incoming_rank'] == 1672
    assert upper['stable_outgoing_rank'] == 1198
    assert 2870 - 1672 - 1198 == upper['cohomology_upper_bound'] == 0

    companion = {}
    if args.positive.is_dir() and (args.positive / 'stable_rank_result.json').is_file():
        stable = json.loads((args.positive / 'stable_rank_result.json').read_text())
        native = json.loads((args.positive / 'native_frames_result.json').read_text())
        rational = json.loads((args.positive / 'rational_omega_result.json').read_text())
        relations = json.loads((args.positive / 'rational_relations_result.json').read_text())
        relation_image = json.loads((args.positive / 'exact_relation_image_result.json').read_text())
        assert stable['formal_centre_dimension'] == 2870
        assert stable['incoming']['rank'] + stable['outgoing']['rank'] == 2870
        assert native['passed'] and rational['passed']
        omega_json = gzip.decompress((args.positive / 'omega_rational.json.gz').read_bytes())
        assert hashlib.sha256(omega_json).hexdigest() == rational['omega_json_sha256']
        assert relations['passed'] and relation_image['passed']
        assert relations['source_basis_dimension'] == 2870
        assert relations['finite_character_dimension'] == 2689
        assert relations['rational_relation_dimension'] == 181
        assert relations['crt_modulus'] > 2 * relations['integer_evaluation_absolute_bound']
        assert relation_image['coefficient_field'] == 'Q'
        assert relation_image['image_rank'] == 116
        assert relation_image['nullity'] == 65
        assert relation_image['all65_exact_null_residuals_zero']
        assert 1672 - (181 - 116) == 1607
        assert native['incoming_minor_rank'] == 1607
        assert native['outgoing_minor_rank'] == 1081
        assert native['omega_pairing'] == 1 and native['boundary_pairings_nonzero'] == 0
        companion = {name: sha(args.positive / name) for name in [
            'stable_rank_result.json', 'native_frames_result.json', 'rational_omega_result.json',
            'rational_relations_result.json', 'exact_relation_image_result.json']}
    table.append({'rank_N': 3, 'dimensions_C8_C9_C10': [2559, 2689, 1281],
                  'incoming_rank': 1607, 'outgoing_rank': 1081, 'H9_dimension': 1,
                  'positive_sector_receipts_checked': bool(companion)})
    for rank in (4, 5):
        table.append({'rank_N': rank, 'dimensions_C8_C9_C10': [dimension_map[(rank, k)] for k in (8, 9, 10)],
                      'incoming_rank': 1672, 'outgoing_rank': 1198, 'H9_dimension': 0})

    result = {'passed': True, 'weight': [3, 3, 5, 0, 0], 'length': 9,
        'highest_weight_only': True, 'prime': prime, 'minor_checks': checks,
        'dimension_and_rank_table': sorted(table, key=lambda row: row['rank_N']),
        'Sp3_rational_rank_proof': {
            'stable_incoming_image_equals_stable_cycle_space': True,
            'incoming_image_intersect_relation_dimension': '181 - 116 = 65',
            'finite_incoming_rank': '1672 - 65 = 1607',
            'omega_nonboundary': 'the modular functional kills every boundary and pairs 1 with omega; the boundary rank is already exact over Q',
            'finite_outgoing_rank': 'at least 1081 by the minor, at most 2689 - 1607 - 1 = 1081 by exact closure and non-exactness',
            'H9_dimension': '2689 - 1607 - 1081 = 1',
            'companion_receipts_checked': bool(companion)},
        'all_N_at_least_4_inference': {
            'source_trace_injectivity_at_N4': 'stable C9 dimension = finite C9 dimension = 2870; specialization is onto',
            'nested_kernel_argument': 'symplectic block restriction gives I_N subset I_4 for N>=4; hence I_N at C9 is zero',
            'incoming_rank': 1672,
            'outgoing_rank': 1198,
            'H9_dimension': 0,
            'outgoing_minor_embedding': 'the Sp4 evaluations embed into Sp(N) by a symplectic direct sum with zero matrices',
            'target_degree_injectivity_assumed': False},
        'input_sha256': {name: sha(root / name) for name in inputs},
        'companion_receipt_sha256': companion,
        'seconds': time.monotonic() - start,
        'scope': 'Saved-minor arithmetic and rank-saturation proof only. Original native generation is recorded separately. The all-N conclusion uses defining-trace generation and symplectic block restriction, and applies only to this fixed classical highest-weight degree. It is not a statement about other grades or quantum corrections.'}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
