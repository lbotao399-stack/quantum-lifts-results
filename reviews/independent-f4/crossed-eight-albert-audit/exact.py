"""Evaluate each pole with the graph's prescribed exterior basis, over Q."""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
from pathlib import Path

from contract import Colour, contract

HERE = Path(__file__).resolve().parent
POLES = [((19, 29, 31, 45, 46, 47), (42, 18, 17, 12, 11, 9)),
         ((19, 20, 29, 45, 46, 47), (42, 39, 18, 12, 11, 9))]


def value(job):
    side, i, graph, reverse, legacy = job
    tree = graph["tree"]
    if reverse:
        tree = [tree[0], list(reversed(tree[1]))]
    colour = Colour(json.loads((HERE / "colour.json").read_bytes()))
    pair = [contract(colour, graph, graph["bits"], tree, pole, canonical=not legacy) for pole in POLES]
    return {"side": side, "graph": i, "poles": list(map(str, pair))}


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--graph-zero", action="store_true", help="Run the order-reversal control only")
    cli.add_argument("--workers", type=int, default=1)
    args = cli.parse_args()
    data = json.loads((HERE / "audit-data.json").read_bytes())
    left, right = (data["graphs"][s] for s in ("left", "right"))
    if args.graph_zero:
        for legacy in (True, False):
            results = [value(("left", 0, left[0], swap, legacy))["poles"] for swap in (False, True)]
            expected = [["-5750", "24200"], ["5750", "-24200"]] if legacy else [["5750", "-24200"]]*2
            assert results == expected
            print(json.dumps({"exterior_order": "tree" if legacy else "graph", "pole_pairs": results}), flush=True)
        return
    assert all(left[i] == right[i] for i in range(12))
    jobs = [("left", i, g, False, False) for i, g in enumerate(left)]
    jobs.extend(("right", i, g, False, False) for i, g in enumerate(right) if i >= 12)
    rows = {"left": {}, "right": {}}
    with ProcessPoolExecutor(args.workers) as executor:
        for result in executor.map(value, jobs):
            side, i = result["side"], result["graph"]
            rows[side][i] = list(map(Fraction, result["poles"]))
            print(json.dumps(result), flush=True)
    rows["right"].update({i: rows["left"][i] for i in range(12)})
    for side in ("left", "right"):
        scores = [g["w"] * (-4*2392*7113106*u + 2*2299*1850212*v)
                  for g, (u, v) in zip(data["graphs"][side], [rows[side][i] for i in range(len(rows[side]))])]
        normalized = sum(scores) / (2*Fraction(28, 3)**4)
        print(json.dumps({"side": side, "raw": str(sum(scores)), "unsubtracted_score": str(normalized)}), flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    main()
