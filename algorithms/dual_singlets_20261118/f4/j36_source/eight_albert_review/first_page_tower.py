"""Reconstruct the Ward equations for the seventh and eighth F4 glueball powers."""

from functools import cache, lru_cache
from itertools import combinations, product
from math import factorial, gcd, lcm, prod
from collections import defaultdict
from pathlib import Path
import json
import sys

from flint import fmpz_mat, fmpq_mat, nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "first_page"))
from colour import Colour, Q, items, wedge, prepend, add, compact, apply  # noqa: E402
from ward import Descent  # noqa: E402
from intermediate_rank import tensor_multiplicity, exchange_trace, pure_current  # noqa: E402


class Harmonics:
    """Enumerate Borel abelian ideals, then solve the finite extremal tensor equations."""

    def __init__(self, c):
        self.c = c
        self.high = {}
        self.labels = {}
        self.low = {}
        self.norm = {}
        self.low_sign = {}
        self.phase = {}
        self.mirror = [(i, Q(v)) for i, v in c.data["longest_weyl"]]
        above = {
            i: {k for e in c.data["raising"] for k in c.br[e][i]} for i in c.positive
        }
        layer = {frozenset()}
        for n in range(1, 10):
            layer = {
                s | {i}
                for s in layer
                for i in c.positive
                if i not in s and above[i] <= s and not any(c.br[i][j] for j in s)
            }
            self.high[n] = []
            ls = []
            for s in layer:
                h = wedge(sorted(s))[0]
                labels = []
                for e, f in zip(c.data["raising"], c.data["lowering"]):
                    a = c.br[e][f]
                    den = sum(x * c.weights[e][k] for k, x in a.items())
                    label = Q(2 * sum(x * c.weight(h)[k] for k, x in a.items()), den)
                    assert label.denominator == 1 and label >= 0
                    labels.append(int(label))
                ls.append((labels, h))
            for j, (labels, h) in enumerate(sorted(ls, reverse=(n == 5))):
                self.high[n].append(h)
                self.labels[n, j] = labels
                t = items(h)
                neg = [c.opposite[i] for i in t]
                self.low[n, j], self.low_sign[n, j] = wedge(neg)
                self.norm[n, j] = prod(c.metric[i, k] for i, k in zip(t, neg))
                self.phase[n, j] = prod(self.mirror[i][1] for i in t)
            if n in c.high:
                assert c.high[n] == self.high[n]

    @cache
    def weight_triples(self):
        dic = defaultdict(list)
        for x in combinations(range(52), 3):
            w = wedge(x)[0]
            dic[self.c.weight(w)].append(w)
        return dic

    def seed_basis(self, n, i, j):
        c = self.c
        weight = tuple(
            a - b for a, b in zip(c.weight(self.high[n][i]), c.weight(self.high[n][j]))
        )
        slots = self.weight_triples().get(weight, [])
        if not slots:
            return []
        rows = defaultdict(dict)
        for k, m in enumerate(slots):
            x, y, z = items(m)
            for a, b, d in [(x, y, z), (y, z, x), (z, x, y)]:
                for bval, coeff in c.br[a][b].items():
                    add(rows["br", bval, d], k, coeff)
            for stage, rs, ls in [
                ("e", c.data["raising"], self.labels[n, j]),
                ("f", c.data["lowering"], self.labels[n, i]),
            ]:
                for t, (r, power) in enumerate(zip(rs, ls)):
                    row = {m: 1}
                    for _ in range(power + 1):
                        row = apply(lambda s: c.act(r, s), row)
                    for a, b in row.items():
                        add(rows[stage, t, a], k, b)
            if i == j:
                moved = [self.mirror[x] for x in items(m)]
                dest, s = wedge([x for x, _ in moved])
                add(rows["mirror", dest], k, int(s * prod(x for _, x in moved)))
                add(rows["mirror", m], k, -1)
        A = fmpz_mat(
            [[int(row.get(k, 0)) for k in range(len(slots))] for row in rows.values()]
        )
        ker, dim = A.nullspace()
        assert (A * ker).is_zero()

        out = []
        for t in range(dim):
            z = [int(ker[k, t]) for k in range(len(slots))]
            d = gcd(*z)
            out.append({a: b // d for a, b in zip(slots, z) if b})
        return out


class Equations(Descent):
    """Use the native metric and Casimir recurrence on each exterior-factor degree."""

    def __init__(self, m):
        self.m = m
        self.n = m - 1
        self.c = c = Colour()
        self.H = h = Harmonics(c)
        self.high_weight = {n: [c.weight(x) for x in v] for n, v in h.high.items()}
        self.low_weight = {
            n: [c.weight(h.low[n, k]) for k in range(len(v))] for n, v in h.high.items()
        }
        self.mirror = h.mirror
        self.colour_duals = [
            [(i, c.metric[i, j]) for i in range(52) if (i, j) in c.metric]
            for j in range(52)
        ]
        self.lam = self.target_kernel()
        self.blocks = {}
        self.cols = {}
        seeds = []
        for i in range(len(h.high[self.n])):
            for j in range(i + 1):
                block = h.seed_basis(self.n, i, j)
                self.blocks[i, j] = len(block)
                for vec in block:
                    seeds.append((i, j, vec))
        self.all_seeds = seeds
        self.d = len(seeds)
        self.zero = fmpq_mat(1, self.d)
        for k, (i, j, vec) in enumerate(seeds):
            for mask, v in vec.items():
                key = i, j, mask
                if key not in self.cols:
                    self.cols[key] = fmpq_mat(1, self.d)
                self.cols[key][0, k] = v

    def target_kernel(self):
        c, h, m = self.c, self.H, self.m
        n = m - 1
        cols = len(h.high[m])
        rows = []
        for i, u in enumerate(h.high[n]):
            for j in range(len(h.high[n])):
                v = h.low[n, j]
                for r, w in enumerate(c.weights):
                    if any(x + y + z for x, y, z in zip(c.weight(u), c.weight(v), w)):
                        continue
                    a = [Q(0) for _ in range(cols)]
                    for pair, z in c.delta[r].items():
                        f, t = items(pair)
                        for e, k, s in [(f, t, z), (t, f, -z)]:
                            U, b = prepend(e, u)
                            V, d = prepend(k, v)
                            if b and d:
                                for slot in range(cols):
                                    a[slot] += b * d * s * self.pure(slot, U, V)
                    den = lcm(*(int(x.denominator) for x in a))
                    rows.append([int(den * x) for x in a])
        ker, size = fmpz_mat(rows).nullspace()
        assert size == 1
        sol = [int(ker[i, 0]) for i in range(cols)]
        div = gcd(*sol) * (1 if next(x for x in sol if x) > 0 else -1)
        sol = [x // div for x in sol]
        self.target_rows = rows
        self.irreps = []
        for u in h.high[m]:
            w = tuple(2 * a + b for a, b in zip(c.weight(u), c.rho2))
            self.irreps.append(
                int(
                    prod(
                        c.scalar_product(w, c.weights[r])
                        / c.scalar_product(c.rho2, c.weights[r])
                        for r in c.positive
                    )
                )
            )
        self.value = (
            (-1) ** (m * (m - 1) // 2)
            * factorial(m)
            * sum(a * b for a, b in zip(self.irreps, sol))
        )
        assert self.value
        if m == 6:
            assert self.value == -15313094596800 and sol == [2392, 2299]
        return sol

    @lru_cache(160000)
    def pure(self, i, u, v):
        c, h, m = self.c, self.H, self.m
        w = c.weight(u)
        if w in self.high_weight[m]:
            return (
                h.norm[m, i] * h.low_sign[m, i]
                if (u, v) == (h.high[m][i], h.low[m, i])
                else 0
            )
        fac = self.denominator(m, w, True)
        if not fac:
            return 0
        total = 0
        for e, f, g in c.roots:
            a = c.act(e, u)
            if a:
                for vv, y in c.act(f, v).items():
                    for uu, x in a.items():
                        total += g * x * y * self.pure(i, uu, vv)
        return fac * total

    @lru_cache(120000)
    def target(self, u, v):
        c = self.c
        if c.scalar_product(c.weight(u), c.rho2) < 0:
            u, v = v, u
        return sum(x * self.pure(i, u, v) for i, x in enumerate(self.lam))

    @lru_cache(25000)
    def seed_row(self, i, j, t):
        h = self.H
        n = self.n
        if i < j:
            u, s = self.transformed(t)
            return (s * h.phase[n, i] / h.phase[n, j]) * self.seed_row(j, i, u)
        total = self.zero
        for mask, v in self.dual(t).items():
            if (i, j, mask) in self.cols:
                total = total + v * self.cols[i, j, mask]
        return total

    @lru_cache(48000)
    def lower_row(self, i, j, v, t):
        c, h, n = self.c, self.H, self.n
        w = c.weight(v)
        if any(a + b + d for a, b, d in zip(w, c.weight(t), self.high_weight[n][i])):
            return self.zero
        if w in self.low_weight[n]:
            return (
                h.low_sign[n, j] * self.seed_row(i, j, t)
                if v == h.low[n, j]
                else self.zero
            )
        fac = self.denominator(n, w, False)
        if not fac:
            return self.zero
        total = self.zero
        for e, f, g in c.roots:
            for vv, x in c.act(f, v).items():
                for tt, y in c.act(e, t).items():
                    total = total + (g * x * y) * self.lower_row(i, j, vv, tt)
        return fac * total

    def probes(self):
        c, h, n = self.c, self.H, self.n
        pairs = {}
        for a, b in combinations(range(52), 2):
            pairs.setdefault(c.weight(wedge((a, b))[0]), []).append((a, b))
        prob = []
        for i, u in enumerate(h.high[n]):
            for j in range(len(h.high[n - 1])):
                low = h.low[n - 1, j]
                for r, w in enumerate(c.weights):
                    need = tuple(
                        -a - b - z for a, b, z in zip(c.weight(u), c.weight(low), w)
                    )
                    prob.extend((i, j, r, a, b) for a, b in pairs.get(need, []))
        return prob

    def row(self, probe):
        i, j, r, a, b = probe
        c, h, n = self.c, self.H, self.n
        u = h.high[n][i]
        v = h.low[n - 1, j]
        active = wedge((a, b))[0]
        part = 0
        hom = self.zero
        for pair, z in c.delta[r].items():
            x, y = items(pair)
            for first, second, coeff in [(x, y, z), (y, x, -z)]:
                vv, s = prepend(first, v)
                t, ss = prepend(second, active)
                if not s or not ss:
                    continue
                coeff *= s * ss
                for k in range(len(h.high[n])):
                    hom = hom + coeff * self.lower_row(i, k, vv, t)
                for (f, d), q in c.mu(t).items():
                    U, us = prepend(f, u)
                    V, vs = prepend(d, vv)
                    if us and vs:
                        part += coeff * q * us * vs * self.target(U, V)
        return (list(hom.entries()) if hasattr(hom, "entries") else hom), Q(part)


def residue(x, p):
    return int(x.numerator) * pow(int(x.denominator), -1, p) % p


class Scalar(Equations):
    def __init__(self, power, sol):
        super().__init__(power)
        self.zero = Q(0)
        self.cols = {}
        for (i, j, data), val in zip(self.all_seeds, sol, strict=True):
            for t, x in data.items():
                add(self.cols, (i, j, t), val * x)


class Mixed(Scalar):
    """Move the adjoint action off the scalar test argument using invariance of the metric."""

    def __init__(self, n, coordinates):
        super().__init__(n, coordinates)
        self.vectors = {}
        for (i, j, seed), coeff in zip(self.all_seeds, coordinates, strict=True):
            obj = self.vectors.setdefault((i, j), {})
            for m, x in seed.items():
                add(obj, m, x * coeff)
        for i, (j, d) in enumerate(self.mirror):
            assert self.mirror[j][0] == i and self.mirror[j][1] * d == 1
        for i in range(len(self.H.high[self.n])):
            for j in range(i + 1, len(self.H.high[self.n])):
                v = {}
                t = self.H.phase[self.n, i] / self.H.phase[self.n, j]
                for m, x in self.vectors[j, i].items():
                    a, s = self.transformed(m)
                    add(v, a, t * s * x)
                self.vectors[i, j] = compact(v)

    @staticmethod
    def add_scaled(out, vec, factor):
        for k, x in vec.items():
            add(out, k, x * factor)

    def add_action(self, out, root, vec, factor):
        for k, x in vec.items():
            for a, y in self.c.act(root, k).items():
                add(out, a, x * y * factor)

    @lru_cache(30000)
    def lower_vector(self, i, j, v):
        c, h, n = self.c, self.H, self.n
        w = c.weight(v)
        if w in self.low_weight[n]:
            return (
                {m: a * h.low_sign[n, j] for m, a in self.vectors[i, j].items()}
                if v == h.low[n, j]
                else {}
            )
        fac = self.denominator(n, w, False)
        if not fac:
            return {}
        total = {}
        for e, f, g in c.roots:
            for vv, x in c.act(f, v).items():
                self.add_action(total, e, self.lower_vector(i, j, vv), -g * x * fac)
        return compact(total)

    @lru_cache(130000)
    def upper_vector(self, i, j, u, v):
        c, h, n = self.c, self.H, self.n
        w = c.weight(u)
        if w in self.high_weight[n]:
            return self.lower_vector(i, j, v) if u == h.high[n][i] else {}
        fac = self.denominator(n, w, True)
        if not fac:
            return {}
        total = {}
        for e, f, g in c.roots:
            for uu, x in c.act(e, u).items():
                for vv, y in c.act(f, v).items():
                    self.add_scaled(
                        total, self.upper_vector(i, j, uu, vv), g * x * y * fac
                    )
                self.add_action(total, f, self.upper_vector(i, j, uu, v), -g * x * fac)
        return compact(total)

    def test(self, i, j, r):
        c, h, n = self.c, self.H, self.n
        u0 = h.high[n - 1][i]
        v0 = h.low[n - 1, j]
        res = {}
        for pair, z in c.delta[r].items():
            a, b = items(pair)
            for x, y, q in ((a, b, z), (b, a, -z)):
                u, s = prepend(x, u0)
                v, t = prepend(y, v0)
                if not s or not t:
                    continue
                if c.scalar_product(c.weight(v), c.rho2) > c.scalar_product(
                    c.weight(u), c.rho2
                ):
                    u, v = v, u
                for i1 in range(len(h.high[n])):
                    for j1 in range(len(h.high[n])):
                        self.add_scaled(res, self.upper_vector(i1, j1, u, v), s * t * q)
        return compact(res)


def independent_multiplicities(equations):
    h = equations.H
    n = equations.n
    (scalar,) = pure_current(3, 0)
    sizes = []
    for i, j in equations.blocks:
        u, v = tuple(h.labels[n, i]), tuple(h.labels[n, j])
        multiplicity = tensor_multiplicity(u, v, scalar)
        dimension = (
            (multiplicity + exchange_trace(u, scalar)) // 2 if i == j else multiplicity
        )
        assert dimension == equations.blocks[i, j]
        sizes.append(dimension)
    return sizes


def ward_rows(system, probes):
    rows = []
    for probe in probes:
        entries, particular = system.row(probe)
        rows.append(entries + [-particular])
    return rows


def det_mod(rows, p):
    return int(nmod_mat([[residue(x, p) for x in row] for row in rows], p).det())


def check_power(power, packet):
    model = Equations(power)
    dimensions = independent_multiplicities(model)
    rows = ward_rows(model, packet["probes"])
    result = dict(
        power=power,
        target_weights=model.lam,
        branch_dimensions=model.irreps,
        target_pairing=model.value,
        homogeneous_dimensions=dimensions,
    )
    before = [r[:-1] for r in rows]
    if power == 8:
        matrix = rows + ward_rows(model, [packet["extra_probe"]])
        minors = {
            str(p): [det_mod(before, p), det_mod(matrix, p)] for p in (1009, 1013)
        }
        assert minors == {"1009": [156, 375], "1013": [382, 29]}
        assert len(matrix) == 92 and len(dimensions) == 6 and sum(dimensions) == 91
        result["homogeneous_and_augmented_ward_minors"] = minors
        return result
    assert sum(dimensions) == 68
    coeffs = fmpq_mat(before).solve(fmpq_mat([[r[-1]] for r in rows])).entries()
    assert list(map(str, coeffs)) == packet["coefficients"]
    mix = Mixed(power, coeffs)
    probes = mix.probes()
    nontrivial = 0
    for probe in probes:
        actual, particular = mix.row(probe)
        assert actual + particular == 0
        nontrivial += bool(particular)
    counts = []
    h = mix.H
    triples = h.weight_triples()
    for i, j in product(range(len(h.high[power - 2])), repeat=2):
        count = 0
        for r, w in enumerate(mix.c.weights):
            weight = tuple(
                -a - b - c
                for a, b, c in zip(
                    w,
                    mix.c.weight(h.high[power - 2][i]),
                    mix.c.weight(h.low[power - 2, j]),
                )
            )
            assert not mix.test(i, j, r)
            count += len(triples.get(weight, []))
        counts.append(count)
    assert (len(probes), nontrivial, sum(counts)) == (1117, 1074, 18948)
    result.update(
        differentiated_scalar_tests=len(probes),
        nonzero_terms_cancelled=nontrivial,
        mixed_components_by_extremal_pair=counts,
    )
    return result


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable Python assertions")
    packet = json.loads((HERE / "tower_probes.json").read_bytes())
    which = [int(p) for p in sys.argv[1:]] or [7, 8]
    checks = Descent().check_psi()
    results = [check_power(p, packet[str(p)]) for p in which]
    print(json.dumps(dict(universal_scalar_check=checks, powers=results), indent=2))
