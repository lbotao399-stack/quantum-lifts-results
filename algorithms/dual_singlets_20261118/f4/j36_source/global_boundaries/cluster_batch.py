"""Extend the same exact evaluation matrix at each sampling checkpoint."""

import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
import numpy as np
from global_source_boundaries import colour, minor, point, write

if __name__ == "__main__":
    cases = (0, 3, 4, 1, 2, 5, 7)
    size, width, prime = 2496, 1008, 1009
    workers = min(int(os.environ["SLURM_CPUS_PER_TASK"]), 48)
    colour(prime)
    data = np.empty((size, len(cases) * width), dtype=np.uint16)
    limits = {320, 768, 1344, 1920, size}
    first = None
    with ProcessPoolExecutor(workers, mp_context=get_context("fork")) as pool:
        tasks = ((i, prime, cases, width) for i in range(size))
        for i, row in enumerate(pool.map(point, tasks, chunksize=2), 1):
            data[i - 1] = row
            if i in limits:
                rows, cols, det = minor(data[:i], prime)
                record = dict(
                    prime=prime,
                    points=i,
                    cases=cases,
                    width=width,
                    rows=rows,
                    columns=cols,
                    determinant=det,
                    rank_lower_bound=len(rows),
                )
                write(Path("global-%d.json" % i), record)
                print(i, len(rows), det, flush=True)
