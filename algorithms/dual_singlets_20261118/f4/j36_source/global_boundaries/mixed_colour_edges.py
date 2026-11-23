"""Trace pairings with a 26-dimensional and an adjoint colour edge across a cut."""

from functools import cache
import numpy as np

from global_source_boundaries import colour, subsets
from marked_jet_blocks import MarkEvaluation, derivative, marked_extensions


@cache
def transition_colour(p):
    gens, frames = colour(p)
    jordan = frames["jo"][0].reshape(26, 26, 26).astype(np.int64)
    # For a word A, store Tr(L_i A T_j) = Tr(T_j L_i A).
    products = (gens[None] @ jordan[:, None]) % p
    coefficients = products.transpose(0, 1, 3, 2).reshape(26 * 52, 676).T
    # Columns of ad_T are defined by [T,T_j] in the fixed native basis.
    trace_dual = frames["ad"][1]
    inverse_gram = frames["ad"][2]
    brackets = (gens[:, None] @ gens[None, :] - gens[None, :] @ gens[:, None]) % p
    maps = (
        (brackets.reshape(52 * 52, 676).astype(float) @ trace_dual % p)
        @ inverse_gram
        % p
    )
    maps = maps.reshape(52, 52, 52).transpose(0, 2, 1)
    assert np.array_equal(
        (maps.transpose(0, 2, 1).reshape(52 * 52, 52) @ frames["ad"][0] % p),
        brackets.reshape(52 * 52, 676) % p,
    )
    # The native basis of the defining representation also labels the Jordan matrices.
    # Check [T,L_i] = L_(T i) before using action on the 26-dimensional edge.
    lhs = (gens[:, None] @ jordan[None, :] - jordan[None, :] @ gens[:, None]).reshape(
        52 * 26, 676
    ) % p
    rhs = gens.transpose(0, 2, 1).reshape(52 * 26, 26) @ frames["jo"][0] % p
    assert np.array_equal(lhs, rhs)
    return coefficients.astype(float), maps.reshape(52, 52 * 52)


class EdgeEvaluation(MarkEvaluation):
    def __init__(self, seed, p, brackets=False):
        super().__init__(seed, p, brackets)
        self.bridge, colour_maps = transition_colour(p)
        flat = self.jets.reshape(len(self.jets), 676)
        coords = self.coordinates(flat, "ad") @ self.frames["ad"][2] % p
        assert np.array_equal(coords @ self.frames["ad"][0] % p, flat)
        self.maps = (coords @ colour_maps % p).reshape(len(self.jets), 52, 52)
        self.rectzero = np.zeros((26, 52))

    def linear_word(self):
        base = self.words()

        @cache
        def value(word, mask):
            token = word[-1]
            if token == "cross":
                m = base(word[:-1], mask)
                return (m.reshape(676) @ self.bridge % self.p).reshape(26, 52)
            side, role, *groups = token
            terms = marked_extensions((role, *groups), mask)
            if not terms:
                return self.rectzero
            if side == "adj":
                left = np.concatenate(
                    [value(word[:-1], old) * sg for old, sg, n in terms], axis=1
                )
                right = np.concatenate([self.maps[n] for old, sg, n in terms])
            else:
                left = np.concatenate(
                    [self.jets[n].T * sg for old, sg, n in terms], axis=1
                )
                right = np.concatenate([value(word[:-1], old) for old, sg, n in terms])
            return (left @ right) % self.p

        def clear():
            base.cache_clear()
            value.cache_clear()

        return value, clear

    def block_edges(self, configuration, indices):
        cuts, left, right, external = configuration
        i, j = np.unravel_index(indices, (len(left), len(right)))
        part_a, ia = np.unique(i, return_inverse=True)
        part_b, ib = np.unique(j, return_inverse=True)
        a = [left[k] for k in part_a]
        b = [right[k] for k in part_b]
        value, clear = self.linear_word()

        def diff(ws, mask):
            return np.stack(
                [
                    sum(
                        (s * value(w, mask) for s, w in differentiated(v)),
                        self.rectzero,
                    )
                    % self.p
                    for v in ws
                ]
            )

        def pair(u, v):
            transformed = (
                (self.frames["jo"][2] @ u % self.p) @ self.frames["ad"][2] % self.p
            )
            return (transformed.reshape(len(u), -1) @ v.reshape(len(v), -1).T)[
                ia, ib
            ] % self.p

        total = np.zeros(len(indices))
        for m, n, sg in subsets(cuts):
            dq = diff(a, m)
            r = np.stack([value(w, n) for w in b])
            out = pair(dq, r)
            if external:
                q = np.stack([value(w, m) for w in a])
                dr = diff(b, n)
                out -= pair(q, dr)
            total += sg * out
            clear()
        return np.asarray(total % self.p, dtype=np.uint16)


@cache
def differentiated(word):
    # Work in the untagged jet algebra; the side of each linear colour action is fixed.
    at = word.index("cross")
    stripped = word[:at] + tuple(token[1:] for token in word[at + 1 :])
    terms = []
    for sign, dw in derivative(stripped):
        after = tuple((word[k + 1][0], *dw[k]) for k in range(at, len(dw)))
        terms.append((sign, dw[:at] + ("cross",) + after))
    return tuple(terms)
