"""Add the third-jet differential and two second jets in one colour word."""

from collections import Counter
from functools import cache
from itertools import combinations
import numpy as np
from global_source_boundaries import (
    Evaluation,
    HEIGHTS,
    arm_bank,
    arm,
    subsets,
)

# One invariant word has both second jets or the third jet. Its partner is unmarked.
CUTS = (
    (4, 0, 3, 0),
    (4, 4, 0, 0),
    (4, 1, 3, 0),
    (4, 0, 3, 1),
    (4, 0, 2, 1),
    (3, 1, 3, 0),
    (4, 1, 2, 0),
    (4, 2, 0, 1),
    (4, 2, 3, 0),
    (3, 0, 3, 1),
)


def encode(pattern, tokens):
    order = [i for t in tokens if not isinstance(t, str) for i in t]
    assert sorted(order) == list(range(len(pattern)))
    inversions = sum(
        i > j and pattern[i] == pattern[j]
        for k, i in enumerate(order)
        for j in order[k + 1 :]
    )
    return (-1) ** inversions, tuple(
        t if isinstance(t, str) else tuple(pattern[i] for i in t) for t in tokens
    )


@cache
def split(spec):
    pattern, tokens = spec
    terms = Counter()
    passed = 0
    for k, t in enumerate(tokens):
        if isinstance(t, str) or len(t) == 1:
            continue
        if len(t) == 2:
            i, j = t
            pieces = [(1, ((i,), (j,))), (-1, ((j,), (i,)))]
        else:
            i, j, h = t
            pieces = []
            for a, b, c in ((i, j, h), (j, h, i), (h, i, j)):
                pieces.extend(((1, ((a,), (b, c))), (-1, ((b, c), (a,)))))
        for sg, words in pieces:
            phase, term = encode(pattern, tokens[:k] + words + tokens[k + 1 :])
            terms[term] += sg * phase * (-1) ** passed
        passed += int(len(t) == 2)
    return tuple((v, k) for k, v in terms.items() if v)


def banks(case):
    rng = np.random.default_rng(713224 + case)
    cuts = CUTS[case]
    other = tuple(h - k for h, k in zip(HEIGHTS, cuts))

    def pattern(weight):
        return tuple(
            int(x)
            for x in rng.permutation(
                [i for i, n in enumerate(weight) for _ in range(n)]
            )
        )

    parents = []
    scalar = []
    for n in range(42):
        pat = pattern(cuts)
        count = len(pat)
        kind = n % 3
        if kind == 0:
            k = int(rng.integers(count - 2))
            blocks = [tuple(range(k, k + 3))]
        else:
            options = [
                (i, j) for i, j in combinations(range(count - 1), 2) if j > i + 1
            ]
            i, j = options[int(rng.integers(len(options)))]
            blocks = [(i, i + 1), (j, j + 1)]
        for i in range(count):
            if all(i not in t for t in blocks):
                blocks.append((i,))
        blocks.sort()
        if rng.integers(2):
            stop = int(rng.integers(2, len(blocks)))
            blocks.insert(stop, ("tr", "jo", "ad")[int(rng.integers(3))])
        parents.append((pat, tuple(blocks)))
    for n in range(65):
        pat = pattern(other)
        words = [(i,) for i in pat]
        if n % 3:
            stop = int(rng.integers(2, len(pat)))
            words.insert(stop, ("tr", "jo", "ad")[int(rng.integers(3))])
        scalar.append(tuple(words))
    return cuts, tuple(parents), tuple(scalar)


banks = cache(banks)


def marked_bank(case):
    rng = np.random.default_rng(60331974 + case)
    cuts = CUTS[case]
    comp = tuple(h - k for h, k in zip(HEIGHTS, cuts))
    # In the shortest unadorned halves use a smaller basis to avoid repeated words.
    left = arm_bank(cuts, 64, rng, True)
    right = arm_bank(comp, 48, rng, True)
    return cuts, left, right


marked_bank = cache(marked_bank)


class JacobiEvaluation(Evaluation):
    def channels(self, a, b):
        n, m = len(a), len(b)
        af = a.reshape(n, 676)
        bf = b.reshape(m, 676)
        out = [af @ b.transpose(0, 2, 1).reshape(m, 676).T]
        for name in ("ad", "jo"):
            gram = self.frames[name][2]
            x = self.coordinates(af, name) @ gram % self.p
            out.append(x @ self.coordinates(bf, name).T)
        out.append(a.trace(axis1=1, axis2=2)[:, None] * b.trace(axis1=1, axis2=2))
        return np.stack(out) % self.p

    def internal(self, case):
        cuts, parents, scalar = banks(case)
        terms = [split(s) for s in parents]
        out = 0
        value = self.evaluator()
        for a, b, phase in subsets(cuts):
            derivatives = np.stack(
                [
                    sum((c * value(w, a) for c, w in row), self.zero) % self.p
                    for row in terms
                ]
            )
            t = np.stack([value(w, b) for w in scalar])
            out += phase * self.channels(derivatives, t)
            value.cache_clear()
        return np.asarray(out % self.p, dtype=np.uint16).ravel()

    def external(self, case):
        cuts, left, right = marked_bank(case)
        out = 0
        value = self.evaluator()

        def frame(bank, mask):
            marked = []
            diff = []
            for spec in bank:
                f, d, r, sg = arm(spec)
                marked.append(value(f, mask))
                diff.append((value(d, mask) - sg * value(r, mask)) % self.p)
            return np.stack(marked), np.stack(diff)

        for a, b, phase in subsets(cuts):
            f, d = frame(left, a)
            g, e = frame(right, b)
            out += phase * (self.channels(d, g) - self.channels(f, e))
            value.cache_clear()
        return np.asarray(out % self.p, dtype=np.uint16).ravel()
