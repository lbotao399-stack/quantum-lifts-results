#!/usr/bin/env python3
"""Freeze or verify the exact F02 packet; never change the roadmap score."""
from pathlib import Path
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    receipt = json.loads((HERE / 'execution_receipt.json').read_text())
    assert receipt['passed'] and receipt['status'] == 'complete'
    assert len(receipt['steps']) == 9
    assert all(step['returncode'] == 0 for step in receipt['steps'])
    for entry in receipt['generated_file_hashes']:
        path = HERE / entry['path']
        assert path.is_file() and digest(path) == entry['sha256'], entry['path']
    manifest_path = HERE / 'MANIFEST.json'
    if args.write:
        files = []
        for path in sorted(HERE.iterdir()):
            if not path.is_file() or path.name == manifest_path.name:
                continue
            if path.suffix not in ('.py', '.cpp', '.md', '.json', '.gz', '.log'):
                continue
            files.append({'path': path.name, 'bytes': path.stat().st_size,
                          'sha256': digest(path)})
        manifest = {
            'scope': 'Complete F02 packet at its explicitly stated leading weak-coupling order; F03 is separate.',
            'public_base_commit': json.loads((HERE / 'declaration.json').read_text())['base_commit'],
            'command_executed_from_repository_root': 'PYTHONPATH=/workspace/scratch/24bd17e8b09b/research_deps OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algorithms/physical_block_20261124/reproduce.py --full',
            'reproduction_mode': receipt['mode'],
            'files': files,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest['files']:
        path = HERE / entry['path']
        assert path.is_file() and path.stat().st_size == entry['bytes'], entry['path']
        assert digest(path) == entry['sha256'], entry['path']
    print(json.dumps({'passed': True, 'file_count': len(manifest['files']),
                      'manifest_sha256': digest(manifest_path),
                      'reproduction_mode': manifest['reproduction_mode']}))

if __name__ == '__main__':
    main()
