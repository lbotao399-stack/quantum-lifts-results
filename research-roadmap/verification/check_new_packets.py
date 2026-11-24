#!/usr/bin/env python3
"""Bounded read-only checks of the new physical packet and accepted ledger."""
import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def zero(matrix):
    return all(matrix[i,j] == 0 for i in range(matrix.nrows()) for j in range(matrix.ncols()))


def physical(output):
    from flint import fmpq_mat
    directory = ROOT / 'algorithms/physical_block_20261124'
    sys.path.insert(0, str(directory))
    import independent_checks
    records = []
    for length, formal, finite, harmonic_count in [(5,35,35,0),(6,427,379,21),(7,816,605,1),(8,261,149,1)]:
        data = read(directory / f'finite_metric_Y{length}.json.gz')
        harmonic = read(directory / f'harmonic_Y{length}.json.gz')
        full = read(directory / f'full_gram_Y{length}.json.gz')['matrix']
        form = fmpq_mat([row[:formal] for row in full[:formal]])
        metric = fmpq_mat(data['physical_gram'])
        embedding = fmpq_mat(data['embedding_of_all_formal_rows'])
        nulls = fmpq_mat(data['null_relations'])
        indices = data['independent_basis_indices']
        assert data['full_formal_dimension'] == formal
        assert data['finite_dimension'] == finite == len(indices)
        assert form == embedding.transpose() * metric * embedding
        assert zero(embedding * nulls)
        assert zero(form * nulls)
        assert len(data['null_relation_indices']) == formal - finite
        for i, row in enumerate(data['null_relation_indices']):
            assert all(nulls[row,j] == (i == j) for j in range(formal-finite))
        assert all(embedding[i,j] == (i == k) for i in range(finite) for k,j in enumerate(indices))
        coefficients = fmpq_mat(harmonic['harmonic_coefficients'])
        outgoing = fmpq_mat(harmonic['outgoing_gram'])
        incoming = fmpq_mat(harmonic['incoming_pairing'])
        assert harmonic['harmonic_dimension'] == harmonic_count == coefficients.ncols()
        assert harmonic['finite_basis_indices'] == indices
        assert incoming == fmpq_mat([[row[j] for j in indices] for row in full[formal:]])
        assert zero(outgoing * coefficients)
        assert zero(incoming * coefficients)
        assert coefficients.transpose() * metric * coefficients == fmpq_mat(harmonic['harmonic_gram'])
        records.append(dict(length=length, formal=formal, finite=finite,
                            all_null_relations=formal-finite, harmonic=harmonic_count))
    full_metric = read(directory / 'physical_metric.json.gz')
    small = full_metric['highest_weight_multiplicity_gram']
    whole = full_metric['complete_flavour_irreducible_representation_gram']
    assert len(small) == 23 and len(whole) == 69
    for i in range(69):
        for j in range(69):
            assert whole[i][j] == (small[i//3][j//3] if i%3 == j%3 else '0')
    # Change only the receipt destination. Native inputs and numerical code
    # remain those pinned by the original47-file physical manifest.
    independent_checks.write = lambda path, value: (output / 'physical-independent.json').write_text(json.dumps(value, indent=2)+'\n')
    independent_checks.main()
    spec = importlib.util.spec_from_file_location('physical_complete_channel', directory / 'f03/channel.py')
    channel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(channel)
    channel.write = lambda path, value: (output / 'physical-channel.json').write_text(json.dumps(value, indent=2)+'\n')
    channel.main()
    result = dict(passed=True, full_formal_equations=records, complete_physical_metric_entries=69*69,
                  complete_source_channel=3, complete_target_channel=2,
                  scope='Leading interacting physical form and complete first quantum-order small channel. Saved large Wick matrices are authenticated and checked, not freshly regenerated here.')
    (output / 'physical-packet.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


def ledger(output):
    directory = ROOT / 'research-roadmap/claim-ledger'
    receipt = read(directory / 'STRICT_GATE_RECEIPT.json')
    assert receipt['A07_complete'] and receipt['ledger_consistent']
    assert receipt['errors'] == receipt['gate_blockers'] == []
    assert len(receipt['negative_controls']) == 9
    assert all(row['rejected'] for row in receipt['negative_controls'])
    for name, expected in receipt['artifact_sha256'].items():
        assert hashlib.sha256((directory / name).read_bytes()).hexdigest() == expected, name
    claims = read(directory / 'claims.json')
    assert len(claims['claims']) == receipt['canonical_claim_count'] == 247
    releases = [read(directory / 'releases' / (name+'.json')) for name in ['b05','physical','f4']]
    assert sum(len(release['documents']) for release in releases) == 12
    assert all(document['review_complete'] and not document['metadata_gaps']
               for release in releases for document in release['documents'])
    coverage = read(directory / 'coverage/coverage_manifest.json')
    for name, expected in coverage['sha256'].items():
        assert hashlib.sha256((directory / 'coverage' / name).read_bytes()).hexdigest() == expected
    assert coverage['expanded_native_group_slots'] == 1376452
    assert coverage['expanded_physical_group_blocks'] == 675081
    result = dict(passed=True, canonical_claim_count=247, reviewed_release_texts=12,
                  accepted_artifacts_checked=len(receipt['artifact_sha256']),
                  inherited_negative_controls=9,
                  fresh_git_reference_audit=False,
                  scope='Authenticate the unchanged complete semantic/reference audit and its exact data. The separate strict checker reruns reference validation when all pinned Git objects are present.')
    (output / 'ledger-artifacts.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet', choices=['physical','ledger'])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (physical if args.packet == 'physical' else ledger)(args.output)
