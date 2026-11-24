#!/usr/bin/env python3
"""Reconstitute the large, unchanged exact physical-block certificate files.

The GitHub transport stores bounded binary chunks. Both each chunk and the
assembled original have their own SHA-256; mathematical bytes are unchanged.
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def assemble(directory=HERE):
    manifest = json.loads((directory / 'ARTIFACTS.json').read_text())
    results = []
    for item in manifest['artifacts']:
        target = (directory / item['path']).resolve()
        if not target.is_relative_to(directory.resolve()):
            raise ValueError('Artifact path leaves the certificate directory.')
        pieces = []
        for part in item['parts']:
            path = (directory / part['path']).resolve()
            if not path.is_relative_to(directory.resolve()):
                raise ValueError('Chunk path leaves the certificate directory.')
            raw = path.read_bytes()
            if len(raw) != part['bytes'] or hashlib.sha256(raw).hexdigest() != part['sha256']:
                raise ValueError('Corrupted artifact chunk: ' + part['path'])
            pieces.append(raw)
        raw = b''.join(pieces)
        if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise ValueError('Incorrect assembled artifact: ' + item['path'])
        created = not target.exists()
        if not created and target.read_bytes() != raw:
            raise ValueError('Existing artifact differs from the immutable certificate: ' + item['path'])
        if created:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
        results.append({'path': item['path'], 'bytes': len(raw), 'sha256': item['sha256'],
                        'parts_verified': len(pieces), 'created': created})
    return results


if __name__ == '__main__':
    print(json.dumps({'verified': True, 'artifacts': assemble()}, indent=2))
