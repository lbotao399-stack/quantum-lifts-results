"""Bound odd G2 BMN cohomology using rational minors of closed Cartan products."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from su3_bmn_parity_products import products_at_length


def calculate(weight, trials):
    table = json.loads((HERE / "5690901/exceptional.json").read_bytes())
    entry = next(row for row in table["rows"] if tuple(row["weight"]) == weight)
    index, = entry["euler"]
    records = []
    for length in range(max(weight), sum(weight)+1):
        if length % 2 == 0:
            row = products_at_length(weight, length, trials, group="G2")
            records.append(row)
            print(weight, length, row["rational_product_rank"], flush=True)
    even = sum(row["rational_product_rank"] for row in records)
    return dict(group="G2", weight=weight, level=2*sum(weight), trials_per_length=trials,
                ordinary_index=index, even_product_lower_bound=even,
                odd_cohomology_lower_bound=max(0, even-index), lengths=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, nargs=3, required=True)
    parser.add_argument("--trials", type=int, default=256)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.trials < 1:
        parser.error("Enable the certificate checks and specify positive trials")
    result = calculate(tuple(args.weight), args.trials)
    args.out.write_text(json.dumps(result, separators=(",", ":"))+"\n")
    print("Odd cohomology lower bound", result["odd_cohomology_lower_bound"])
