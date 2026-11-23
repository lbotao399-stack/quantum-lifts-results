"""Evaluate both first marks of the three six-cubic pivots in the order-703 minor."""
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
import json
from pathlib import Path

from contract import Colour, contract
from exact import POLES

HERE = Path(__file__).resolve().parent


@cache
def colour():
    return Colour(json.loads((HERE / 'colour.json').read_bytes()))


def evaluate(job):
    i, graph = job
    poles = [contract(colour(), graph, graph['bits'], graph['tree'], p) for p in POLES]
    return {'graph': i, 'poles': list(map(str, poles))}


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('--workers', type=int, default=1)
    args = cli.parse_args()
    data = json.loads((HERE / 'section-graphs.json').read_bytes())
    scores = []
    with ProcessPoolExecutor(args.workers) as pool:
        for row in pool.map(evaluate, enumerate(data['graphs'])):
            assert row['graph'] == len(scores)
            u, v = map(Fraction, row['poles'])
            scores.append(-4*2392*7113106*u + 2*2299*1850212*v)
            print(json.dumps(row), flush=True)
    for row, banks in data['rows'].items():
        normalized = [sum((w*scores[i] for i,w in bank), Fraction()) / (2*Fraction(28,3)**3)
                      for bank in banks]
        assert normalized[0] == normalized[1]
        print(json.dumps({'pivot': int(row), 'left': str(normalized[0]), 'right': str(normalized[1])}), flush=True)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    main()
