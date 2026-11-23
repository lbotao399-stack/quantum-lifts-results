"""Exact native characters beyond the two certified symplectic sectors."""
from pathlib import Path
import argparse
import json
import os
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parent
sys.path.insert(0, str(WORKSPACE / 'fortuity_discovery_20261123/jets'))
import pfaffian_relations
from trace_descent import dimensions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ranks', type=int, nargs='+', default=[6, 7])
    parser.add_argument('--seconds', type=int, default=240)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('character cap')))
    signal.alarm(args.seconds)
    start = time.monotonic()
    path = HERE / 'character_extension.json'
    result = json.loads(path.read_text()) if path.exists() else {'records': []}
    done = {(r['rank'], r['length']) for r in result['records']}
    try:
        for rank in args.ranks:
            weight = (3, 3, 2 * rank - 1, 0, 0)
            for length in range(2 * rank + 2, 2 * rank + 5):
                if (rank, length) in done:
                    continue
                stable = sum(dimensions(weight, length, 2 * (length + 1), 'Sp').values())
                finite = sum(dimensions(weight, length, 2 * rank, 'Sp').values())
                record = dict(rank=rank, native_weight=weight, length=length,
                              stable=stable, finite=finite, identities=stable-finite,
                              elapsed_seconds=time.monotonic()-start)
                result['records'].append(record)
                path.write_text(json.dumps(result, indent=2)+'\n')
                print(json.dumps(record), flush=True)
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
