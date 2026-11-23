"""Direct Pauli-free trace test of the two virtual edges and their classical variation."""

from itertools import permutations, product
import json
import numpy as np
from global_source_boundaries import colour, OFFSETS
from marked_jet_blocks import PAIRS
from mixed_colour_edges import EdgeEvaluation, differentiated
from mixed_source_search import model


def test(p):
    e = EdgeEvaluation(15903, p)
    native, frames = colour(p)
    jordan = frames["jo"][0].reshape(26, 26, 26)

    def comm(x, y):
        return (x @ y - y @ x) % p

    def y(i, j):
        if i == j:
            return e.zero
        return (1 if i < j else -1) * e.jets[4 + PAIRS.index(tuple(sorted((i, j))))] % p

    def literal(w, assigned, i, j):
        pos = w.index("cross")
        its = [iter(v) for v in assigned]
        values = []
        for k, token in enumerate(t for t in w if t != "cross"):
            role, *groups = token if k < pos else token[1:]
            ids = [next(its[g]) for g in groups]
            if role == "x":
                values.append((e.jets[ids[0]], None, 0))
                continue
            if role == "y":
                a, b = ids
                values.append((y(a, b), comm(e.jets[a], e.jets[b]), 1))
                continue
            a, b, c = ids
            values.append(
                (
                    e.zero,
                    (
                        comm(e.jets[a], y(b, c))
                        + comm(e.jets[b], y(c, a))
                        + comm(e.jets[c], y(a, b))
                    )
                    % p,
                    0,
                )
            )
        total = 0
        parity = 0
        for k, (v, d, odd) in enumerate(values):
            if d is not None:
                mats = [z[0] for z in values]
                mats[k] = d
                prefix = e.eye
                for mat in mats[:pos]:
                    prefix = prefix @ mat % p
                jcol = np.eye(26)[i].copy()
                a_col = native[j].astype(float) % p
                for tok, mat in zip(reversed(w[pos + 1 :]), reversed(mats[pos:])):
                    if tok[0] == "adj":
                        a_col = comm(mat, a_col)
                    else:
                        jcol = mat @ jcol % p
                ji = (jcol @ jordan.reshape(26, 676) % p).reshape(26, 26)
                total += (-1) ** parity * int(((ji @ prefix % p) @ a_col % p).trace())
            parity = (parity + odd) % 2
        return total % p

    result = []
    for case in (0, 1):
        cuts, left, _, _ = model(case)
        labels = [tuple(range(n)) for n in cuts]
        mask = sum(1 << (OFFSETS[g] + i) for g, c in enumerate(cuts) for i in range(c))
        get, clear = e.linear_word()
        for n in range(len(left)):
            word = left[n]
            raw = (
                sum(
                    (s * get(term, mask) for s, term in differentiated(word)),
                    e.rectzero,
                )
                % p
            )
            hits = np.argwhere(raw)
            if not len(hits):
                continue
            i, j = map(int, hits[len(hits) // 2])
            total = 0
            for assn in product(*(permutations(v) for v in labels)):
                inversions = sum(
                    a > b
                    for seq in assn
                    for m, a in enumerate(seq)
                    for b in seq[m + 1 :]
                )
                total += (-1) ** inversions * literal(word, assn, i, j)
            assert total % p == int(raw[i, j])
            result.append([case, n, i, j, total % p])
            if sum(r[0] == case for r in result) == 2:
                break
        clear()
    assert len(result) == 4
    return result


if __name__ == "__main__":
    assert __debug__
    with np.errstate(all="ignore"):
        print(json.dumps({str(p): test(p) for p in (1009, 1013)}, indent=2))
