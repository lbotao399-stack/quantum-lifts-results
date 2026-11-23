#!/usr/bin/env python3
"""Regenerate and prove a 40-term rational fortuity detector for Sp(3) omega.

The search uses one finite field; the reported detector is subsequently
verified over the rationals on the COMPLETE certified relation space.
There is no rational-reconstruction assumption in this final verification.
"""
from collections import defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import pickle
import time

import numpy as np
from scipy.sparse import csr_matrix, load_npz
from flint import fmpz_mat, nmod_mat

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_rank(array):
    return fmpz_mat(np.asarray(array, dtype=np.int64).tolist()).rank()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=HERE.parent / 'positive_sector')
    ap.add_argument('--output', type=Path, default=HERE / 'sp3_intrinsic_dual.json')
    args = ap.parse_args()
    start = time.monotonic()
    root = args.input
    proof_path = root / 'rational_relations_result.json'
    basis_path = root / 'rational_relation_basis.json.gz'
    proof = json.loads(proof_path.read_text())
    assert proof['passed']
    assert proof['crt_modulus'] > 2 * proof['integer_evaluation_absolute_bound']
    rational_bytes = gzip.decompress(basis_path.read_bytes())
    assert hashlib.sha256(rational_bytes).hexdigest() == proof['rational_basis_sha256']
    rational = json.loads(rational_bytes)
    assert rational['complete_relation_dimension'] == 181
    modular = np.load(root / 'relation_image.npz')
    q = load_npz(root / 'stable_q_columns.npz').tocsr()
    assert q.shape == (6020, 2870)
    prime = 1000003

    # This is only a deterministic SEARCH. Its result is checked below over Q.
    rows = modular['minor_rows']
    augmented = np.vstack((modular['images'][rows], modular['omega_defect'])).astype(np.int64)
    reduced, rank = nmod_mat(augmented.tolist(), prime).rref()
    assert rank == 117
    columns = [next(j for j in range(6020) if reduced[i, j]) for i in range(117)]
    rhs = nmod_mat(117, 1, prime)
    rhs[116, 0] = 1
    solution = nmod_mat(augmented[:, columns].tolist(), prime).solve(rhs)
    ell = np.zeros(6020, dtype=np.int64)
    for i, column in enumerate(columns):
        coefficient = int(solution[i, 0])
        ell[column] = coefficient if coefficient <= prime // 2 else coefficient - prime
    assert np.count_nonzero(ell) == 40

    # Clear each certified rational relation's denominators solely for exact
    # integer arithmetic. The source cochains and the detector are not rescaled.
    relation_entries, relation_rows, relation_columns = [], [], []
    denominators = []
    for i, row in enumerate(rational['basis_coefficients']):
        fractions = [(j, Fraction(c)) for j, c in row]
        denominator = lcm(*(c.denominator for _, c in fractions))
        denominators.append(denominator)
        for j, coefficient in fractions:
            relation_rows.append(i)
            relation_columns.append(j)
            relation_entries.append(int(coefficient * denominator))
    relations_integer = csr_matrix((np.asarray(relation_entries, dtype=np.int64),
                                   (relation_rows, relation_columns)), shape=(181, 2870))
    # Prove that the fast signed-int64 products below cannot wrap. The ranks
    # themselves are subsequently computed by arbitrary-precision fmpz.
    maximum_q = max(abs(int(x)) for x in q.data)
    maximum_relation_l1 = max(sum(abs(int(x)) for x in relations_integer.data[
        relations_integer.indptr[i]:relations_integer.indptr[i + 1]]) for i in range(181))
    assert maximum_relation_l1 * maximum_q < 2**63
    relation_images_integer = (relations_integer @ q.T).toarray()
    assert int(np.max(np.abs(relation_images_integer))) * 6020 * int(max(abs(ell))) < 2**63
    rational_pairings_numerators = relation_images_integer @ ell
    assert not np.any(rational_pairings_numerators)

    # The actual original rational cycle, in the original finite-basis order.
    cycle_path = root / 'inputs/algorithms/dual_singlets_20261118/sp3_omega/cycle.json.gz'
    original_cycle = json.loads(gzip.decompress(cycle_path.read_bytes()))
    cycle_coefficients = [(i, Fraction(c)) for i, c in original_cycle['cycle_basis_coefficients']]
    omega_denominator = lcm(*(c.denominator for _, c in cycle_coefficients))
    omega_integer = np.zeros(2870, dtype=np.int64)
    for i, coefficient in cycle_coefficients:
        omega_integer[modular['basis_rows'][i]] = int(coefficient * omega_denominator)
    assert int(max(abs(omega_integer))) * maximum_q * int(max(np.diff(q.indptr))) < 2**63
    omega_defect_integer = np.asarray(q @ omega_integer).reshape(-1)
    assert int(max(abs(omega_defect_integer))) * 6020 * int(max(abs(ell))) < 2**63
    omega_numerator = int(ell @ omega_defect_integer)
    omega_pairing = Fraction(omega_numerator, omega_denominator)
    assert omega_pairing == 1

    # Verify that the original rational formal lift is the independently
    # exported 4264-term omega, without changing its representative or scale.
    _, centre, _, outgoing = pickle.loads((root / 'stable_rows.pkl').read_bytes())
    omega = defaultdict(Fraction)
    for i, coefficient in cycle_coefficients:
        for word, value in centre[modular['basis_rows'][i]].items():
            omega[word] += coefficient * value
    omega = {word: value for word, value in omega.items() if value}
    exported = json.loads(gzip.decompress((root / 'omega_rational.json.gz').read_bytes()))
    def decode(word):
        return tuple(tuple((i, j, tuple(t)) for i, j, t in trace) for trace in word)
    expected = {decode(word): Fraction(value) for value, word in exported['terms']}
    assert omega == expected and len(omega) == 4264

    words = sorted({word for row in outgoing for word in row})
    assert len(words) == 6020
    shapes = defaultdict(list)
    for i, word in enumerate(words):
        shapes[tuple(sorted(map(len, word)))].append(i)
    shape_results = []
    for shape, indices in sorted(shapes.items()):
        restricted = relation_images_integer[:, indices]
        augmented_integer = np.vstack((restricted, omega_defect_integer[indices]))
        image_rank = exact_rank(restricted)
        joint_rank = exact_rank(augmented_integer)
        contribution = Fraction(int(ell[indices] @ omega_defect_integer[indices]), omega_denominator)
        shape_results.append({'trace_lengths': list(shape), 'coordinates': len(indices),
            'relation_image_rank_over_Q': image_rank,
            'relation_image_plus_omega_rank_over_Q': joint_rank,
            'chosen_dual_support': int(np.count_nonzero(ell[indices])),
            'chosen_dual_pairing_contribution': str(contribution)})
    assert exact_rank(relation_images_integer) == 116
    assert exact_rank(np.vstack((relation_images_integer, omega_defect_integer))) == 117
    assert sum(Fraction(row['chosen_dual_pairing_contribution']) for row in shape_results) == 1

    result = {
        'passed': True,
        'group': 'Sp(3)=USp(6)',
        'weight': [3, 3, 5, 0, 0],
        'length': 9,
        'level': 28,
        'coefficient_field': 'Q',
        'source_normalization': exported['coefficient_convention'],
        'search_prime_only': prime,
        'detector_integer_support': 40,
        'maximum_absolute_detector_coefficient': int(max(abs(ell))),
        'full_relation_basis_dimension': 181,
        'full_relation_image_rank_over_Q': 116,
        'augmented_relation_image_rank_over_Q': 117,
        'all_relation_pairings_exact_zero': True,
        'omega_pairing': str(omega_pairing),
        'cleared_omega_denominator_internal_only': omega_denominator,
        'cleared_omega_pairing_numerator': omega_numerator,
        'original_omega_trace_terms': len(omega),
        'formal_source_count': 2870,
        'formal_output_coordinate_count': 6020,
        'detector': [{'index': int(i), 'integer_coefficient': int(ell[i]),
                      'ordered_native_traces': words[i]} for i in np.flatnonzero(ell)],
        'trace_partition_results': shape_results,
        'certified_rational_relation_basis': {
            'compressed_sha256': digest(basis_path),
            'decoded_sha256': hashlib.sha256(rational_bytes).hexdigest(),
            'complete_finite_frame_dimension': proof['finite_character_dimension'],
            'frame_minor_determinant': proof['full_finite_frame_minor_determinant'],
            'crt_modulus': proof['crt_modulus'],
            'integer_evaluation_absolute_bound': proof['integer_evaluation_absolute_bound']},
        'input_sha256': {name: digest(root / name) for name in [
            'rational_relations_result.json', 'rational_relation_basis.json.gz',
            'stable_q_columns.npz', 'relation_image.npz', 'stable_rows.pkl',
            'omega_rational.json.gz',
            'inputs/algorithms/dual_singlets_20261118/sp3_omega/cycle.json.gz']},
        'seconds': time.monotonic() - start,
        'scope': 'Exact rational annihilation of Q0 on every finite-colour trace identity in the complete source degree, and exact pairing1 with the unchanged rational omega lift. Certifies a nonzero fortuitous component independently of choosing a graviton/harmonic basis. This is a new certificate for the known omega class, not a new BPS class, physical observable, quantum lift or minimal-support theorem.'
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('detector', 'input_sha256')}, indent=2))


if __name__ == '__main__':
    main()
