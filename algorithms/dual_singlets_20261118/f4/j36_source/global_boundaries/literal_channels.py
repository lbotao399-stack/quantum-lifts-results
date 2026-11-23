"""Check projected boundary words by explicit column permutations and the graded Leibniz rule."""

from itertools import permutations, product
import json
import numpy as np
from colour_channel_search import recipe, schedules
from marked_jet_blocks import MarkEvaluation, PAIRS, group_degrees


def evaluate_word(e, word, assignment):
    cursors = [iter(v) for v in assignment]

    def commutator(a, b):
        return (a @ b - b @ a) % e.p

    def letter(label):
        role, *groups = label
        indices = [next(cursors[g]) for g in groups]
        if role == "x":
            return e.jets[indices[0]], e.zero, 0
        if len(set(indices)) < len(indices):
            return e.zero, e.zero, int(role == "y")

        def y(i, j):
            s = 1 if i < j else -1
            return s * e.jets[4 + PAIRS.index(tuple(sorted((i, j))))] % e.p

        if role == "y":
            i, j = indices
            return y(i, j), commutator(e.jets[i], e.jets[j]), 1
        i, j, k = indices
        jac = (
            commutator(e.jets[i], y(j, k))
            + commutator(e.jets[j], y(k, i))
            + commutator(e.jets[k], y(i, j))
        ) % e.p
        return e.zero, jac, 0

    def parent(tokens):
        value, variation, parity = e.eye, e.zero, 0
        for t in tokens:
            if t[0] == "part":
                v, d, odd = parent(t[2])
                v = e.decorate(v, t[1])
                d = e.decorate(d, t[1])
            else:
                v, d, odd = letter(t)
            variation = (variation @ v + (-1) ** parity * (value @ d)) % e.p
            value = (value @ v) % e.p
            parity = (parity + odd) % 2
        return value, variation, parity

    return parent(word)


def alternate(e, w):
    value = e.zero.copy()
    diff = e.zero.copy()
    perms = [list(permutations(range(k))) for k in group_degrees(w)]
    for choices in product(*perms):
        sign = (-1) ** sum(
            i > j
            for sequence in choices
            for k, i in enumerate(sequence)
            for j in sequence[k + 1 :]
        )
        v, d, parity = evaluate_word(e, w, choices)
        value = (value + sign * v) % e.p
        diff = (diff + sign * d) % e.p
    return value, diff, parity


def check(p):
    e = MarkEvaluation(15903, p)
    answer = []
    for case in (0, 1):
        cuts, left, right, external = recipe(case)
        # There is no split colour column in these two cases. Each half can be alternated independently.
        choices = schedules()[case]
        for ch in range(5):
            options = choices[choices // (len(left) * len(right)) == ch]
            if not len(options):
                continue
            quick = e.evaluate_bank(cuts, left, right, external, options)
            pos = int(np.flatnonzero(quick)[0])
            idx = int(options[pos])
            _, i, j = np.unravel_index(idx, (5, len(left), len(right)))
            a, da, parity = alternate(e, left[i])
            b, db, _ = alternate(e, right[j])

            def pair(x, y):
                if ch == 0:
                    return (x @ y).trace()
                if ch == 3:
                    return x.trace() * y.trace()
                if ch == 4:
                    return (e.transpose(x) @ y).trace()
                return (e.decorate(x, ("ad", "jo")[ch - 1]) @ y).trace()

            direct = int((pair(da, b) + (-1) ** parity * pair(a, db)) % p)
            assert direct == int(quick[pos]) and direct
            answer.append([case, ch, idx, direct])
    assert {row[1] for row in answer} == set(range(5))
    return answer


if __name__ == "__main__":
    assert __debug__
    with np.errstate(all="ignore"):
        print(json.dumps({str(p): check(p) for p in (1009, 1013)}, indent=2))
