"""Adjacent rational rank bounds for the possible USp(4) level-32 scalar dressing."""
import argparse
from pathlib import Path

from so7_hair_bound import projected_image


def run(length, out, workers):
    # Either pair of lower bounds adding to 13,639 proves H at length 11 vanishes.
    # The two stopping ranks are forecasts obtained by inserting the eleven
    # known even products into the Euler recursion. They are never used as data.
    if length == 10:
        widths, stop, seed = {1: 1304, 3: 4200, 5: 3000}, 4896, 893437
    elif length == 11:
        widths, stop, seed = {0: 352, 2: 5100, 4: 7200}, 8743, 914897
    else:
        raise ValueError("expected source length 10 or 11")
    return projected_image((0, 0, 8, 4, 4), length, 4, widths, stop, seed,
                           out, workers, "image.json.gz", checkpoint=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, choices=(10, 11), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run(args.length, args.out, args.workers)
