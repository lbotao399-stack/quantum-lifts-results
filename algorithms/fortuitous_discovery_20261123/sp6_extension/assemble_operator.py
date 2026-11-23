#!/usr/bin/env python3
"""Reassemble the original verified gzip bytes from small Git transport parts."""
import argparse
import hashlib
import json
from pathlib import Path


def operator_bytes(data_directory):
    data_directory=Path(data_directory)
    manifest=json.loads((data_directory/'sp6_operator_parts.json').read_text())
    parts=[]
    for record in manifest['parts']:
        content=(data_directory/record['file']).read_bytes()
        assert len(content)==record['bytes']
        assert hashlib.sha256(content).hexdigest()==record['sha256']
        parts.append(content)
    content=b''.join(parts)
    assert len(content)==manifest['original_bytes']
    assert hashlib.sha256(content).hexdigest()==manifest['original_sha256']
    return content


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    content=operator_bytes(Path(__file__).resolve().parent/'data')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes(content)
    print(json.dumps(dict(output=str(args.output),bytes=len(content),
                         sha256=hashlib.sha256(content).hexdigest())))


if __name__=='__main__':main()
