"""Decompose the complete SO(7) incoming source by field type and jet order."""
import argparse
import json
import sys
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path

from check_incoming import RUNTIME, cochain_dimension

sys.path.insert(0, str(RUNTIME))
from g2_full_classical_screen_v1 import contents
from exact_dims import _multiplicity
from fortuity.full.letters import Q_letter, is_ferm
from sp_full_profiles import profile


MU = (3, 3, 2, 2, 2)


def dimensions(mu, length):
    answer = Counter()
    rho = (2, 1, 0)
    for offset, sign in ((0, 1), (1, -1)):
        for order in permutations(range(3)):
            inversions = sum(order[i] > order[j] for i in range(3) for j in range(i+1, 3))
            flavours = tuple(mu[k+2] + rho[k] - rho[order[k]] for k in range(3))
            shifted = (mu[0]+offset, mu[1]-offset) + flavours
            for content in contents(shifted, length):
                counts = Counter(content)
                bosons = tuple(sorted((n for l, n in counts.items() if not is_ferm(l)), reverse=True))
                fermions = tuple(sorted((n for l, n in counts.items() if is_ferm(l)), reverse=True))
                answer[profile(content)] += sign * (-1)**inversions * _multiplicity(bosons, fermions, 0, 7)
    assert min(answer.values()) >= 0 and sum(answer.values()) == cochain_dimension(mu, length)
    return {p: n for p, n in sorted(answer.items()) if n}


def records(dims):
    return [dict(index=i, profile=p, dimension=n, fermions=sum(v for (s, d), v in p if s % 2 == 0))
            for i, (p, n) in enumerate(dims.items())]


def build():
    sources = dimensions(MU, 10)
    targets = dimensions(MU, 11)
    source_index = {p: i for i, p in enumerate(sources)}
    target_index = {p: i for i, p in enumerate(targets)}
    edges = defaultdict(set)
    for content in contents(MU, 10):
        p = profile(content)
        if p in sources:
            for pos, letter in enumerate(content):
                for _, left, right in Q_letter(letter):
                    out = profile(content[:pos] + content[pos+1:] + (left, right))
                    if out in targets:
                        edges[source_index[p]].add(target_index[out])
    src = records(sources)
    for row in src:
        row["targets"] = sorted(edges[row["index"]])
    assert sum(sources.values()) == 42991 and sum(targets.values()) == 14832
    return dict(mu=MU, length=10, sources=src, targets=records(targets))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: [[r["index"], r["fermions"], r["dimension"], r.get("targets")]
                            for r in result[key]] for key in ("sources", "targets")}))
