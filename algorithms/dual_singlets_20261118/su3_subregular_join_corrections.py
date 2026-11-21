"""Combine two rational first corrections and replay their common scalar-chart equation."""
import argparse
import json
from pathlib import Path

from su3_subregular_chart_replay import verify
from su3_subregular_image import write
from su3_subregular_root_coupled import d_variables, read_coefficients
from su3_subregular_root_linear import decode, root_degree


def join(high, linear, through_d):
    higher = read_coefficients(decode(high))
    lower = read_coefficients(decode(linear))
    if any(len(d_variables(w)) < through_d for w in higher):
        raise ValueError("The high-d correction contains a lower component")
    if any(root_degree(w) > 1 for w in lower):
        raise ValueError("The linear correction contains a higher off-diagonal ghost component")
    higher.update((w, c) for w, c in lower.items() if len(d_variables(w)) < through_d)
    return [[w.hex(), str(c)] for w, c in sorted(higher.items())]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image", type=Path, required=True)
    p.add_argument("--high", type=Path, required=True)
    p.add_argument("--linear", type=Path, required=True)
    p.add_argument("--through-d", type=int, choices=(1, 2, 3), required=True)
    p.add_argument("--primitive", type=Path)
    p.add_argument("--residual", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("Enable assertions for the source and coefficient checks")
    recipe = join(args.high, args.linear, args.through_d)
    reference = Path(__file__).resolve().parent / "taylor_hooks/su3_scalar_levi/subregular"
    result = verify(args.image, recipe, 4, args.through_d, reference, also_through_roots=1,
                    residual_path=args.residual)
    if args.primitive:
        write(args.primitive, recipe)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
