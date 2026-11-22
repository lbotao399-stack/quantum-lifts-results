"""Evaluate the polarized Pfaffian trace logarithm with a fixed physical slot."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from verify import PHI, GHOSTS, so_canon
from orthogonal_glueball_second_sources import MATRICES


@lru_cache(None)
def orientations(word):
    alphabet = PHI + GHOSTS

    def rotations(seq):
        return [seq[i:] + seq[:i] for i in range(len(seq))]

    forward = rotations(word)
    orbits = [word] if word[::-1] in forward else [word, word[::-1]]
    factor = Fraction(len(set(forward)), 2 * len(word))
    refs = []
    for original in orbits:
        cuts = []
        for j in range(len(word)):
            w = original[-j:] + original[:-j] if j else original
            sign, canonical = so_canon((tuple(alphabet[k] for k in w),))
            assert sign and canonical == (tuple(alphabet[k] for k in word),)
            cuts.append((w, sign))
        refs.append(cuts)
    return refs, factor


def attach_rooted(rel):
    longest = max(map(len, rel.catalog))
    # One outer minus gives (-1)^(number of tiles+1). Each tile supplies one minus.
    table = defaultdict(list)
    for block, choices in rel.catalog.items():
        for b, bit, degree, used, c in choices:
            for i, j, v in MATRICES[b]:
                table[block, i].append((j, bit, degree, used, -c * v))
    good = rel.project.divisors

    @lru_cache(None)
    def trace(word):
        ors, factor = orientations(word)
        total = Counter()
        for cuts in ors:
            for overlap, (w, orientation) in enumerate(cuts[:longest]):
                n = len(w)
                states = [defaultdict(int) for _ in range(n + 1)]
                states[0][0, 0, 0, 0, 0] = states[0][1, 1, 0, 0, 0] = 1
                for i in range(n):
                    if not states[i]:
                        continue
                    for length in range(
                        1 + overlap if i == 0 else 1, min(longest, n - i) + 1
                    ):
                        block = w[i : i + length]
                        target = states[i + length]
                        adj = (table.get((block, 0), ()), table.get((block, 1), ()))
                        if not any(adj):
                            continue
                        for (start, end, m, x, z), v in states[i].items():
                            if not v:
                                continue
                            for nxt, b, y, t, c in adj[end]:
                                if m & b or z & t or x + y not in good:
                                    continue
                                s = -1 if b and (m // (b * 2)).bit_count() % 2 else 1
                                target[start, nxt, m | b, x + y, z | t] += v * c * s
                for (a, b, m, x, z), v in states[-1].items():
                    if a == b and v:
                        total[m, x, z] -= orientation * v
        return {k: v * factor for k, v in total.items() if v}

    rel.trace.cache_clear()
    rel.trace = trace
