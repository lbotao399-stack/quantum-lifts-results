#!/usr/bin/env python3
"""Restore exact reviewed ledger bytes from bounded, authenticated XZ chunks."""
import hashlib
import json
import lzma
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent


def authenticated(data, size, digest):
    if len(data) != size or hashlib.sha256(data).hexdigest() != digest:
        raise ValueError('Ledger artifact length or SHA-256 mismatch.')
    return data


def main():
    manifest = json.loads((HERE / 'ARTIFACTS.json').read_text())
    if manifest['schema'] != 'exact-concatenated-xz-v1':
        raise ValueError('Unsupported ledger transport schema.')
    packed = b''.join(authenticated((HERE / item['path']).read_bytes(), item['bytes'], item['sha256'])
                      for item in manifest['chunks'])
    authenticated(packed, manifest['xz_bytes'], manifest['xz_sha256'])
    raw = lzma.decompress(packed, format=lzma.FORMAT_XZ, memlimit=256 * 1024 * 1024)
    authenticated(raw, manifest['raw_bytes'], manifest['raw_sha256'])
    offset = 0
    for item in manifest['files']:
        path = Path(item['path'])
        if path.is_absolute() or '..' in path.parts or item['offset'] != offset:
            raise ValueError('Invalid or noncontiguous ledger artifact path/offset.')
        target = HERE / path
        data = authenticated(raw[offset:offset + item['bytes']], item['bytes'], item['sha256'])
        offset += item['bytes']
        if target.exists():
            authenticated(target.read_bytes(), item['bytes'], item['sha256'])
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(target.name + '.assembling')
        temporary.write_bytes(data)
        os.replace(temporary, target)
    if offset != len(raw):
        raise ValueError('Unaccounted bytes in ledger artifact archive.')
    print(json.dumps({'verified': True, 'files': len(manifest['files']),
                      'chunks': len(manifest['chunks']), 'raw_bytes': len(raw)}))


if __name__ == '__main__':
    main()
