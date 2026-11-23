"""Classical images of twenty middle alternating triples for one pair of exchanges."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import gzip
from itertools import combinations_with_replacement
import json
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import integer_section as integer
from exchange_scores import PAIRS


def point(i, pair):
    with np.errstate(all="ignore"):
        x = (integer.sample(i) @ integer.native26().reshape(52, 676)).reshape(4, 26, 26)
        low, high, exchange, scales = integer.tensors()
        assert scales == (28, 252, 84)
        right = integer.halves(x, high, [integer.CROSSED[2:]]).reshape((26,) * 4)
        for j in pair:
            right = integer.exact_apply(integer.action(exchange), right, PAIRS[j])
        middles = list(combinations_with_replacement(range(4), 3))
        # Work in batches to keep the integer limb arrays small.
        vals = []
        for k in range(0, len(middles), 5):
            arms = [(integer.CROSSED[0], m) for m in middles[k : k + 5]]
            left = -integer.halves(x.transpose(0, 2, 1), low, arms)
            vals.extend(map(int, integer.gram(left, right.reshape(1, -1)).ravel()))
    return [i, vals]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--indices", type=int, nargs=2, required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        data = list(
            pool.map(partial(point, pair=args.indices), range(720), chunksize=2)
        )
    args.output.write_bytes(
        gzip.compress(json.dumps(data, separators=(",", ":")).encode(), mtime=0)
    )
