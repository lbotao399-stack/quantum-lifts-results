"""Evaluate alternating third-jet and two-second-jet boundaries directly on colour matrices."""

from functools import cache
from fractions import Fraction
import gzip
import json
from itertools import combinations, product
import numpy as np

from global_source_boundaries import (
    Evaluation,
    HERE,
    HEIGHTS,
    OFFSETS,
    PAIRS,
    arm,
    inverse,
    layout,
    subsets,
)
from global_source_jacobi import banks, marked_bank
from global_source_scan import FAMILIES as JACOBI_FAMILIES, previous, sample_columns

TRIPLES = tuple(combinations(range(4), 3))


@cache
def marked_extensions(token, mask):
    role, *groups = token
    available = [
        [i for i in range(HEIGHTS[g]) if (mask >> (OFFSETS[g] + i)) & 1] for g in groups
    ]
    terms = []
    for labels in product(*available):
        if len(set(labels)) != len(labels):
            continue
        previous_mask = mask
        inversions = sum(x > y for j, x in enumerate(labels) for y in labels[j + 1 :])
        for g, i in reversed(tuple(zip(groups, labels))):
            previous_mask ^= 1 << (OFFSETS[g] + i)
            within = (previous_mask >> OFFSETS[g]) & ((1 << HEIGHTS[g]) - 1)
            inversions += (within >> (i + 1)).bit_count()
        if role == "x":
            index = labels[0]
        elif role == "j":
            index = 16 + TRIPLES.index(tuple(sorted(labels)))
        else:
            index = (4 if role == "y" else 10) + PAIRS.index(tuple(sorted(labels)))
        terms.append((previous_mask, (-1) ** inversions, index))
    return terms


def typed(tokens):
    return tuple(
        t if isinstance(t, str) else (("x" if len(t) == 1 else "y"), *t) for t in tokens
    )


@cache
def derivative(word):
    """Apply Q0 on the ordered predecessor before doing any Taylor alternation."""
    terms = []
    phase = 1
    for k, token in enumerate(word):
        if isinstance(token, str):
            continue
        if token[0] == "part":
            _, channel, inside = token
            for sign, value in derivative(inside):
                terms.append(
                    (
                        phase * sign,
                        word[:k] + (("part", channel, value),) + word[k + 1 :],
                    )
                )
            phase *= (-1) ** jet_parity(inside)
        elif token[0] == "z":
            terms.append((phase, word[:k] + (("j", *token[1:]),) + word[k + 1 :]))
        elif token[0] == "y":
            terms.append((phase, word[:k] + (("k", *token[1:]),) + word[k + 1 :]))
            phase = -phase
    return tuple(terms)


@cache
def jet_parity(word):
    return (
        sum(
            jet_parity(t[2]) if t[0] == "part" else int(t[0] == "y")
            for t in word
            if not isinstance(t, str)
        )
        % 2
    )


@cache
def group_degrees(word):
    groups = []
    for t in word:
        if isinstance(t, str):
            continue
        if t[0] == "part":
            weights = group_degrees(t[2])
            groups.extend(j for j, k in enumerate(weights) for _ in range(k))
        else:
            groups.extend(t[1:])
    return tuple(groups.count(g) for g in range(4))


@cache
def branch_splits(mask, degrees):
    available = [
        [i for i in range(HEIGHTS[g]) if mask >> (OFFSETS[g] + i) & 1] for g in range(4)
    ]
    terms = []
    for groups in product(
        *(combinations(labels, k) for labels, k in zip(available, degrees))
    ):
        right = sum(
            1 << (OFFSETS[g] + i) for g, group in enumerate(groups) for i in group
        )
        inversions = sum(
            i > j
            for whole, chosen in zip(available, groups)
            for i in whole
            if i not in chosen
            for j in chosen
        )
        terms.append((mask ^ right, right, (-1) ** inversions))
    return terms


@cache
def invariant_metric(p):
    h = np.zeros((26, 26), dtype=np.int64)
    representation = json.loads(
        gzip.decompress((HERE / "representation26.json.gz").read_bytes())
    )
    for i, j, x in representation["form"]:
        q = Fraction(x)
        h[i, j] = q.numerator * pow(q.denominator, -1, p) % p
    return h.astype(float), inverse(h, p).astype(float)


@cache
def indexed_banks(family, case):
    if family == "base":
        cuts, a, b = layout(case)
    elif family == "out":
        cuts, a, b = marked_bank(case)
    else:
        cuts, a, b = banks(case)
        parents = []
        for pattern, tokens in a:
            words = tuple(
                t
                if isinstance(t, str)
                else (("x", "y", "z")[len(t) - 1], *(pattern[i] for i in t))
                for t in tokens
            )
            parents.append(words)
        return cuts, tuple(parents), tuple(typed(w) for w in b)
    return cuts, tuple(typed(arm(x)[0]) for x in a), tuple(typed(arm(x)[0]) for x in b)


class MarkEvaluation(Evaluation):
    def __init__(self, seed, prime, brackets=False):
        super().__init__(seed, prime, brackets)
        jets = list(self.g)

        def comm(a, b):
            return (a @ b - b @ a) % prime

        jets.extend(comm(jets[i], jets[j]) for i, j in PAIRS)
        for i, j, k in TRIPLES:
            first = comm(jets[i], jets[4 + PAIRS.index((j, k))])
            second = comm(jets[j], jets[4 + PAIRS.index((i, k))])
            third = comm(jets[k], jets[4 + PAIRS.index((i, j))])
            jets.append((first - second + third) % prime)
        self.jets = np.stack(jets)
        # One product includes no more than 24 signed terms of inner dimension 26.
        assert 24 * 676 * prime**2 < 2**53

    def transpose(self, x):
        h, inv = invariant_metric(self.p)
        return (inv @ x.swapaxes(-2, -1) % self.p) @ h % self.p

    def decorate(self, x, name):
        if name in ("sym", "skew"):
            partner = self.transpose(x)
            if name == "skew":
                return (
                    (x - partner) * pow(2, -1, self.p) - super().decorate(x, "ad")
                ) % self.p
            return (
                (x + partner) * pow(2, -1, self.p)
                - super().decorate(x, "jo")
                - self.eye * (x.trace() * pow(26, -1, self.p) % self.p)
            ) % self.p
        return super().decorate(x, name)

    def words(self):
        @cache
        def evaluate(word, mask):
            if not word:
                assert not mask
                return self.eye
            prefix, token = word[:-1], word[-1]
            if isinstance(token, str):
                return self.decorate(evaluate(prefix, mask), token)
            if token[0] == "part":
                _, channel, inner = token
                terms = branch_splits(mask, group_degrees(inner))
                if not terms:
                    return self.zero
                left = np.concatenate(
                    [evaluate(prefix, l) * s for l, r, s in terms], axis=1
                )
                right = np.concatenate(
                    [projected(inner, r, channel) for l, r, s in terms]
                )
                return left @ right % self.p
            terms = marked_extensions(token, mask)
            if not terms:
                return self.zero
            if len(terms) == 1:
                m, sign, index = terms[0]
                return sign * (evaluate(prefix, m) @ self.jets[index]) % self.p
            left = np.concatenate([evaluate(prefix, m) for m, _, _ in terms], axis=1)
            right = np.concatenate(
                [self.jets[index] * sign for _, sign, index in terms]
            )
            return (left @ right) % self.p

        @cache
        def projected(word, mask, channel):
            return self.decorate(evaluate(word, mask), channel)

        clear = evaluate.cache_clear

        def reset():
            clear()
            projected.cache_clear()

        evaluate.cache_clear = reset
        return evaluate

    def selected(self, family, case, indices):
        cuts, a, b = indexed_banks(family, case)
        return self.evaluate_bank(cuts, a, b, family != "in", indices)

    def evaluate_bank(self, cuts, a, b, external, indices):
        nch = int(max(indices) // (len(a) * len(b))) + 1
        chan, ai, bi = np.unravel_index(indices, (nch, len(a), len(b)))
        use_a, ia = np.unique(ai, return_inverse=True)
        use_b, ib = np.unique(bi, return_inverse=True)
        aw = [a[i] for i in use_a]
        bw = [b[i] for i in use_b]

        def derivatives(ws, mask):
            return np.stack(
                [
                    sum(
                        (sign * evaluate(term, mask) for sign, term in derivative(w)),
                        self.zero,
                    )
                    % self.p
                    for w in ws
                ]
            )

        out = np.zeros(len(indices))
        evaluate = self.words()
        for l, r, sg in subsets(cuts):
            da = derivatives(aw, l)
            vb = np.stack([evaluate(w, r) for w in bw])
            if external:
                va = np.stack([evaluate(w, l) for w in aw])
                db = derivatives(bw, r)
            for c in np.unique(chan):
                take = np.flatnonzero(chan == c)

                def pairing(x, y):
                    if c == 3:
                        return (
                            x.trace(axis1=1, axis2=2)[ia[take]]
                            * y.trace(axis1=1, axis2=2)[ib[take]]
                        )
                    xf = x.reshape(len(x), -1)
                    if c in (0, 4):
                        first = xf if c == 0 else self.transpose(x).reshape(len(x), -1)
                        second = y.transpose(0, 2, 1).reshape(len(y), -1)
                    else:
                        name = ("ad", "jo")[c - 1]
                        gram = self.frames[name][2]
                        first = self.coordinates(xf, name) @ gram % self.p
                        second = self.coordinates(y.reshape(len(y), -1), name)
                    return (first @ second.T)[ia[take], ib[take]]

                v = pairing(da, vb)
                if external:
                    v -= pairing(va, db)
                out[take] += sg * (v % self.p)
            evaluate.cache_clear()
        return np.asarray(out % self.p, dtype=np.uint16)


@cache
def retained(certificate):
    import json
    from pathlib import Path

    old = np.asarray(json.loads(Path(certificate).read_bytes())["columns"])
    answer = []
    for case, pos, indices in previous():
        keep, src = np.nonzero(old[old < 2398, None] == pos[None, :])
        # Old pivots were stored in the order chosen by the elimination.
        answer.append(("base", case, np.flatnonzero(old < 2398)[keep], indices[src]))
    for k, (kind, case) in enumerate(JACOBI_FAMILIES):
        pos = np.flatnonzero((old >= 2398 + 850 * k) & (old < 2398 + 850 * (k + 1)))
        if len(pos):
            answer.append(
                (kind, case, pos, sample_columns(kind, case)[(old[pos] - 2398) % 850])
            )
    return answer


def old_row(evaluation, certificate):
    groups = retained(certificate)
    row = np.zeros(sum(len(pos) for _, _, pos, _ in groups), dtype=np.uint16)
    for kind, case, pos, indices in groups:
        row[pos] = evaluation.selected(kind, case, indices)
    return row
