#!/usr/bin/env python3
"""Run the declared public verification tier, with a real corruption control.

This verifies bounded public certificates and solver controls. It does not
claim a fresh reconstruction of an entire production census. Receipts and
temporary files live outside the checkout, which must start and end clean.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*arguments):
    return subprocess.check_output(['git', '-C', str(ROOT), *arguments], text=True).strip()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_manifest(path):
    manifest = json.loads(path.read_text())
    entries = manifest.get('files', manifest)
    if isinstance(entries, dict):
        entries = [{'path': name, 'sha256': value} for name, value in entries.items()]
    for entry in entries:
        item = path.parent / entry['path']
        require(item.is_file(), 'Missing manifest input: ' + str(item.relative_to(ROOT)))
        require(digest(item) == entry['sha256'], 'Corrupt manifest input: ' + str(item.relative_to(ROOT)))
    return {'path': str(path.relative_to(ROOT)), 'sha256': digest(path), 'entries_checked': len(entries)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--allow-dirty', action='store_true', help='Development only; cannot produce a clean-checkout receipt.')
    args = parser.parse_args()
    require(__debug__, 'Python assertions must remain enabled.')
    output = args.output.resolve()
    require(not output.is_relative_to(ROOT), 'Put receipts outside the immutable checkout.')
    output.mkdir(parents=True, exist_ok=True)
    initial_status = git('status', '--porcelain', '--untracked-files=all')
    require(args.allow_dirty or not initial_status, 'A clean checkout is required.')
    started = time.monotonic()
    receipt = {
        'schema_version': 1, 'tier': 'public-bounded-verification-v2',
        'commit': git('rev-parse', 'HEAD'), 'tree': git('rev-parse', 'HEAD^{tree}'),
        'clean_at_start': not bool(initial_status), 'development_run': args.allow_dirty,
        'python': sys.version, 'platform': platform.platform(),
        'packages': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'sympy', 'python-flint')},
        'manifests': [], 'checks': [],
    }
    receipt.update(passed=False, status='running')
    receipt_path = output / 'receipt.json'
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    for relative in (
        'algorithms/fortuitous_duals_20261123/MANIFEST.json',
        'algorithms/su3_second_correction_20261118/SOURCE_MANIFEST.json',
        'algorithms/su3_independent_prime_20261118/SOURCE_MANIFEST.json',
    ):
        receipt['manifests'].append(check_manifest(ROOT / relative))

    environment = dict(os.environ)
    environment.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    environment.pop('PYTHONOPTIMIZE', None)

    def run(identifier, arguments, scope, expected_exit=0, required_text=()):
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        begin = time.monotonic()
        completed = subprocess.run(arguments, cwd=ROOT, env=environment, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        log = output / (identifier + '.log')
        log.write_text(completed.stdout)
        require((completed.returncode == 0) if expected_exit == 0 else (completed.returncode != 0),
                identifier + ' had an unexpected exit; see ' + str(log))
        for expected in required_text:
            require(expected in completed.stdout, identifier + ' omitted expected result: ' + expected)
        record = {
            'id': identifier, 'command': arguments, 'scope': scope,
            'exit_code': completed.returncode, 'wall_seconds': time.monotonic() - begin,
            'total_processor_seconds': (after.ru_utime + after.ru_stime) - (before.ru_utime + before.ru_stime),
            'largest_child_peak_kib_so_far': after.ru_maxrss,
            'log': log.name, 'log_sha256': digest(log), 'required_text': list(required_text),
        }
        receipt['checks'].append(record)
        print(json.dumps({'id': identifier, 'passed': True, 'seconds': record['wall_seconds']}), flush=True)

    python = sys.executable
    run('roadmap', [python, '-B', str(ROOT / 'research-roadmap/score.py')],
        'Weights, dependencies, evidence pins and generated documents; no physics inferred from this validator.')
    run('sp3-coefficient-dual', [python, '-B', str(ROOT / 'algorithms/fortuitous_duals_20261123/transpose_detector/score_candidate.py')],
        'Exact coefficient pairing on the unchanged published omega; closure/completeness are inherited certificate evidence.',
        required_text=('"exact_coefficient_pairing": "1"', '"source_coordinates": 186'))

    with tempfile.TemporaryDirectory(prefix='public-verification-') as temporary:
        work = Path(temporary)
        source = ROOT / 'algorithms/fortuitous_duals_20261123/so11_independent_audit/inputs/orthogonal_glueball_first_pair.py'
        audit = ROOT / 'algorithms/fortuitous_duals_20261123/so11_independent_audit/global_return_audit.py'
        (work / 'so11').mkdir()
        run('so11-global-return', [python, '-B', str(audit), '--pfaffian-source', str(source), '--out-dir', str(work / 'so11')],
            'Independent complete 1000-marker return identity and original coefficient oracle; no new source census.')
        for path in (work / 'so11').glob('*.json'):
            shutil.copy2(path, output / ('so11-' + path.name))
        damaged = work / 'damaged_pfaffian_source.py'
        damaged.write_bytes(source.read_bytes() + b'\n# Deliberate public-verification corruption.\n')
        run('corrupted-source-rejected', [python, '-B', str(audit), '--pfaffian-source', str(damaged), '--out-dir', str(work / 'so11')],
            'A copied input has different bytes; the unchanged native verifier must reject its pinned hash.',
            expected_exit=1, required_text=('The Pfaffian source differs from the pinned, audited file',))
        receipt['corruption_control'] = {'original_sha256': digest(source), 'damaged_sha256': digest(damaged),
                                         'original_unchanged': True, 'rejected': True}

        # The historical control writes a receipt beside its script. Execute an
        # unchanged private working copy so no tracked certificate is rewritten.
        old = ROOT / 'algorithms/su3_second_correction_20261118/runtime'
        copied = work / 'runtime'
        shutil.copytree(old, copied, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        runtime = copied / 'work/sp3-native.nosync'
        run('su3-complete-affine-control', [python, '-B', str(runtime / 'validate_su3_second_coupled_solve_v1.py')],
            'Complete 30-equation synthetic control, all five homogeneous directions, restart, incompatible-frame and coverage rejection.',
            required_text=('"homogeneous_dimension": 5', '"rank": 11', 'genuine_inconsistency_agrees_with_independent_dense_rank'))
        shutil.copy2(runtime / 'validate_su3_second_coupled_solve_v1.json', output / 'su3-affine-receipt.json')
        package = ROOT / 'algorithms/su3_independent_prime_20261118'
        run('su3-lossless-packing', [python, '-B', str(package / 'validate_compact_su3_maps_v1.py'),
                                   '--runtime', str(runtime), '--output', str(output / 'su3-packing-receipt.json')],
            'Original and packed complete control equations agree; both missing rows and a modified matrix entry are rejected.',
            required_text=('corrupted_payload_rejected', 'missing_rows_rejected'))

    generator = ROOT / 'algorithms/su3_generator_20261124'
    run('su3-generator-assemble', [python, '-B', str(generator / 'assemble_artifacts.py')],
        'Lossless assembly of the three original certificate files; every chunk and full-file hash is verified.',
        required_text=('"verified": true',))
    run('su3-generator-certificate', [python, '-B', str(generator / 'verify.py'), '--rational-ranks',
                                     '--receipt', str(output / 'su3-generator-certificate.json')],
        'All complete certificate ranks over the rationals and every dual pairing; archived finite-frame entries reused.',
        required_text=('"verified": true', '"rational_rank_recalculation": true'))
    run('su3-generator-independent', [python, '-B', str(generator / 'independent_verify.py'),
                                     '--receipt', str(output / 'su3-generator-independent.json')],
        'Separate polynomial differentiation, symmetry checks and evaluation of every residual at all 1588 complete target coordinates.',
        required_text=('"full_target_coordinates_checked": 1588', '"full_polynomial_residual_nonzero": 0'))

    ledger = ROOT / 'research-roadmap/claim-ledger'
    run('ledger-assemble', [python, '-B', str(ledger / 'assemble_artifacts.py')],
        'Restore the exact reviewed ledger bytes from authenticated bounded chunks.',
        required_text=('"verified": true',))
    checker = HERE / 'check_new_packets.py'
    run('ledger-accepted-artifacts', [python, '-B', str(checker), 'ledger', '--output', str(output)],
        'Authenticate the complete strict-audit receipt, its unchanged data and expanded finite coverage. The full Git-reference audit is separately callable.',
        required_text=('"passed": true', '"fresh_git_reference_audit": false'))
    physical = ROOT / 'algorithms/physical_block_20261124'
    run('physical-assemble', [python, '-B', str(physical / 'assemble_artifacts.py')],
        'Restore every exact physical certificate payload and verify all original file hashes.',
        required_text=('"verified": true',))
    receipt['manifests'].append(check_manifest(physical / 'MANIFEST.json'))
    receipt['manifests'].append(check_manifest(physical / 'f03/MANIFEST.json'))
    run('physical-complete-packet', [python, '-B', str(checker), 'physical', '--output', str(output)],
        'Every finite null and physical harmonic equation, all 4761 metric entries, independent normalization/Wick checks and the full 3-to-2 quantum channel. Stored large Wick evaluations are reused.',
        required_text=('"passed": true', '"complete_physical_metric_entries": 4761'))
    f4 = ROOT / 'algorithms/f4_complete_source_20261124'
    run('f4-publication-manifest', [python, '-B', str(f4 / 'publication.py')],
        'Authenticate the exact permitted F4 packet allowlist; stopped-route source is excluded.')
    inputs = output / 'f4-inputs'
    run('f4-public-inputs', [python, '-B', str(f4 / 'survival_inputs.py'), '--download', '--destination', str(inputs)],
        'Fetch only 100 explicitly pinned public input blobs and authenticate every byte.')
    run('f4-complete-small-identities', [python, '-B', str(f4 / 'verify_packet.py'), '--inputs', str(inputs), '--out', str(output / 'f4-packet.json')],
        'Recompute the full rational direct-product and 60-column bilinear identities; authenticate the 150850-entry extension. The historical complete 3017 minor remains attributed.',
        required_text=('"complete_second_map_rank": 0', '"fresh_native_complete_kernel_minor_rebuild": false'))
    final_status = git('status', '--porcelain', '--untracked-files=all')
    require(final_status == initial_status, 'The verification job modified the checkout.')
    receipt.update(clean_at_end=not bool(final_status), passed=True, status='complete',
                   full_production_sources_rebuilt=False, wall_seconds=time.monotonic() - started,
                   runner_sha256=digest(Path(__file__)))
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'checks': len(receipt['checks']), 'clean': receipt['clean_at_start'] and receipt['clean_at_end']}))


if __name__ == '__main__':
    main()
