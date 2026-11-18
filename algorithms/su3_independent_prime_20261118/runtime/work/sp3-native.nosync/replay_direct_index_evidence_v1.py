"""Replay Cartan minors and archived rows with the dense matrix evaluator.

This checks every archived determinant and independent production probes.
It does not reconstruct every entry of every production minor.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import time

import sp_census
import numpy as np
from flint import nmod_mat
from hw_products import Raising
from sp_graviton_singlets import product_values
from g2_cartan_product_reuse_v1 import torus_certificate


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--full-minors', action='store_true')
    args = parser.parse_args()
    start = time.monotonic()
    manifest = json.loads((args.evidence / 'SHA256SUMS.json').read_text())
    for name, digest in manifest.items():
        assert sha(args.evidence / name) == digest, name
    families, closure_cache, checked_rows, records = {}, set(), set(), []
    full_minor_coordinates = 0
    certificate = json.loads(json.dumps(torus_certificate()))
    for group in ('Sp3', 'G2'):
        for probe_path in sorted((args.evidence / group).glob('s*/probes.pkl')):
            folder = probe_path.parent
            result = json.loads((folder / 'result.json').read_text())
            payload = pickle.loads(probe_path.read_bytes())
            assert payload['result_sha256'] == sha(folder / 'result.json')
            family_path = args.evidence / 'families' / (folder.name + '.pkl')
            assert sha(family_path) == result.get('family_sha256', result.get('universal_source_family_sha256'))
            if folder.name not in families:
                families[folder.name] = pickle.loads(family_path.read_bytes())
            family = families[folder.name]
            prime = result['prime']
            sp_census.configure(3, prime)
            for generator in family['generators']:
                fingerprint = repr(sorted(generator['expr'].items()))
                if fingerprint not in closure_cache:
                    assert not sp_census.base.action_expr(generator['expr'])
                    closure_cache.add(fingerprint)
            minor = payload['minor']
            rank = result['rank_lower_bound']
            assert minor.shape == (rank, rank)
            determinant = int(nmod_mat(minor.tolist(), prime).det()) if rank else 1
            assert determinant == result['minor_determinant'] and determinant
            rows = payload['selected_rows']
            coefficients = family['kernel'][rows]
            indices = sorted(set(map(int, coefficients.indices)))
            coefficients = coefficients[:, indices]
            recipes = [family['products'][index] for index in indices]
            raising = Raising(family['generators'], prime)
            for source_index in rows:
                identity = folder.name, source_index
                if identity in checked_rows:
                    continue
                row = family['kernel'][source_index]
                for root in (1, 2):
                    residual = defaultdict(int)
                    for index, value in zip(row.indices, row.data):
                        for word, coefficient in raising.product(family['products'][int(index)], root).items():
                            residual[word] += int(value) * coefficient
                    assert not any(value % prime for value in residual.values()), identity
                checked_rows.add(identity)
            coordinates = 0
            for probe in payload['probes']:
                frame = sp_census.base.Frame(probe['degree'], probe['points'], 1, prime, k=probe['k'])
                frame.L = probe['letters']
                if group == 'G2':
                    assert result['torus_certificate'] == certificate
                raw = product_values(family['generators'], recipes, frame, probe['degree'])
                actual = np.asarray(coefficients @ raw) % prime
                assert np.array_equal(actual, probe['values']), (group, folder.name, probe['degree'])
                coordinates += actual.size
            if args.full_minors and rank:
                full_coefficients = family['kernel'][result['minor_rows']]
                needed = sorted(set(map(int, full_coefficients.indices)))
                full_coefficients = full_coefficients[:, needed]
                full_recipes = [family['products'][index] for index in needed]
                blocks = []
                for metadata in result['frames']:
                    frame = sp_census.base.Frame(metadata['degree'], metadata['points'], metadata['seed'],
                                                prime, k=metadata['k'])
                    for letter, (odd, matrix) in list(frame.L.items()):
                        diagonal = np.zeros_like(matrix)
                        for index in range(3):
                            diagonal[..., index, index] = matrix[..., index, index]
                        if group == 'G2':
                            diagonal[..., 2, 2] = (-diagonal[..., 0, 0] - diagonal[..., 1, 1]) % prime
                        for index in range(3):
                            diagonal[..., index + 3, index + 3] = -diagonal[..., index, index] % prime
                        frame.L[letter] = odd, diagonal
                    raw = product_values(family['generators'], full_recipes, frame, metadata['degree'])
                    blocks.append(np.asarray(full_coefficients @ raw) % prime)
                actual_minor = np.concatenate(blocks, axis=1)[:, result['minor_columns']]
                assert np.array_equal(actual_minor, minor), (group, folder.name, 'full_minor')
                full_minor_coordinates += actual_minor.size
            records.append(dict(group=group, sector=result['sector'], rank_lower_bound=rank,
                                minor_determinant=determinant, independent_rows=len(rows),
                                dense_matrix_probe_coordinates=int(coordinates)))
            if len(records) % 10 == 0:
                print(json.dumps(dict(stage='replay', completed=len(records), seconds=time.monotonic() - start)), flush=True)
    result = dict(status='passed_all_payload_hashes_minors_and_dense_matrix_probes',
                  payload_hashes=len(manifest), records=records, rational_G2_torus_certificate=certificate,
                  exact_generator_closure_checks=len(closure_cache),
                  modular_raising_rows=len(checked_rows),
                  independent_dense_matrix_coordinates=sum(row['dense_matrix_probe_coordinates'] for row in records),
                  full_minor_coordinates_reconstructed=int(full_minor_coordinates),
                  all_nonzero_minor_entries_reconstructed=args.full_minors,
                  nonempty_minors=sum(row['rank_lower_bound'] > 0 for row in records),
                  driver_sha256=sha(__file__), seconds=time.monotonic() - start,
                  scope=('Every payload hash, every nonzero minor entry and determinant independently reconstructed with the dense matrix trace evaluator, production probes, exact generator closure, and selected modular raising equations.' if args.full_minors else
                         'Every payload hash and nonzero minor determinant, selected production rows with the independent dense matrix evaluator, exact generator closure, and selected modular raising equations. Full production minor entries were not all independently reconstructed.'))
    (args.evidence / 'replay.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('records', 'rational_G2_torus_certificate')}), flush=True)


if __name__ == '__main__':
    main()
