"""Abelian gauge + polynomial W: Q0, Q1, Q2, Q3 over Q (exact) or F_p.

Port of work/intrinsic_higher_lift_2026-09-19/engine.py (Theory, wick,
selected, qn, graphmoment) with three changes only:
  1. coefficient ring is a parameter (P=0: Fraction; P>0: integers mod P);
  2. per-chiral weights (W must be homogeneous of degree D in them);
  3. Q_n, n >= 3, uses higher_Q_article_2026-09-20/higher_moments.fanmoment(n,...)
     with the same ordered-edge prefactor (-1)^n 2^(n-1); for deg W<=3
     every graph of Q3 is a fan (article Prop. 2.1(i) / 3.10).  For
     deg W = 4 the three-leaf book (quartic hub, three cubic pages) is added
     as in higher_Q_article_2026-09-20/quartic_q3.py (bookwick, bookmoment,
     prefactor -4/3); gate: gate_q3_book.py.  deg W >= 5: Q3 refuses.
  4. Q_n, n >= 4, for deg W <= 3 only: every vertex has valence <= 3, so the
     internal tree is a path and every graph is a fan; gate: gate_q4.py.
Physical Q_n = (64 pi^2)^(-n) times stored Q_n (as in the source engine).
Pure standard library (no numpy).
"""
from fractions import Fraction as F
from collections import defaultdict
from functools import lru_cache
from itertools import product, permutations, combinations
from math import factorial

P_DEFAULT = 2147483647


class Ring:
    def __init__(self, P=0):
        self.P = P

    def __call__(self, x):
        P = self.P
        if not P:
            return F(x)
        x = F(x)
        d = x.denominator % P
        if d == 0:
            raise ZeroDivisionError('denominator divisible by p')
        return x.numerator % P * pow(d, P - 2, P) % P

    def inv(self, a):
        return 1 / a if not self.P else pow(a, self.P - 2, self.P)


class Theory:
    def __init__(self, charges, W, P=0, weights=None, D=None):
        self.R = Ring(P)
        self.P = P
        self.charges = charges
        self.rank = len(next(iter(charges.values())))
        self.gs = tuple(charges)
        self.bs = {g: 'B' + g for g in self.gs}
        self.odd = set(self.bs.values()) | {'c' + str(a) for a in range(self.rank)}
        self.fields = tuple(sorted(self.gs + tuple(self.bs.values()) + tuple('b' + str(a) for a in range(self.rank)) + tuple('c' + str(a) for a in range(self.rank))))
        self.partners = {**self.bs, **{b: g for g, b in self.bs.items()}, **{'b' + str(a): 'c' + str(a) for a in range(self.rank)}, **{'c' + str(a): 'b' + str(a) for a in range(self.rank)}}
        self.wt = dict(weights) if weights else {g: 1 for g in self.gs}
        degs = {sum(self.wt[g] for g in w) for w in W}
        if D is None:
            assert len(degs) == 1, ('W not homogeneous in the given weights', degs)
            D = degs.pop()
        else:
            assert degs <= {D}, ('W not homogeneous of degree D', degs, D)
        self.D = D
        assert all(0 < self.wt[g] < D for g in self.gs), 'need 0<wt(gamma)<D'
        self.degW = max((len(w) for w in W), default=0)
        R = self.R
        self.W = {}
        for w, c in W.items():
            k = tuple(sorted((g, 0, 0) for g in w))
            self.W[k] = self.W.get(k, 0) + R(c)
        self.W = self.fin(self.W)
        self.vertices = [(tuple(g for g, r, s in w), c) for w, c in self.W.items()]
        for g, qs in charges.items():
            for a, q in enumerate(qs):
                if q:
                    self.vertices.append(((self.bs[g], 'c' + str(a), g), R(q)))
        self.base = {}
        for g, qs in charges.items():
            self.base[g] = self.add(*(self.scale(self.mul(self.v('c' + str(a)), self.v(g)), q) for a, q in enumerate(qs)))
            self.base[self.bs[g]] = self.add(self.fd(self.W, g), *(self.scale(self.mul(self.v('c' + str(a)), self.v(self.bs[g])), -q) for a, q in enumerate(qs)))
        for a in range(self.rank):
            self.base['c' + str(a)] = {}
            self.base['b' + str(a)] = self.add(*(self.scale(self.mul(self.v(self.bs[g]), self.v(g)), -qs[a]) for g, qs in charges.items()))

    # ---------- ring helpers
    def fin(self, out):
        P = self.P
        if P:
            return {w: c % P for w, c in out.items() if c % P}
        return {w: c for w, c in out.items() if c}

    def parity(self, u):
        return int(u[0] in self.odd)

    def canon(self, w):
        os_ = [u for u in w if self.parity(u)]
        if len(os_) != len(set(os_)):
            return (), 0
        sg = (-1) ** sum(self.parity(a) * self.parity(b) for i, a in enumerate(w) for b in w[i + 1:] if a > b)
        return tuple(sorted(w)), sg

    def add(self, *ps):
        out = defaultdict(int)
        for p in ps:
            for w, c in p.items():
                out[w] += c
        return self.fin(out)

    def scale(self, p, c):
        c = self.R(c) if isinstance(c, F) else c
        return self.fin({w: x * c for w, x in p.items()})

    def mul(self, *ps):
        out = {(): self.R(1)}
        for p in ps:
            new = defaultdict(int)
            for a, c in out.items():
                for b, d in p.items():
                    w, s = self.canon(a + b)
                    if s:
                        new[w] += s * c * d
            out = self.fin(new)
        return out

    def v(self, g, r=0, s=0):
        return {((g, r, s),): self.R(1)}

    def deriv(self, p, axis):
        out = defaultdict(int)
        for w, c in p.items():
            for i, u in enumerate(w):
                v = list(u)
                v[axis + 1] += 1
                k, s = self.canon(w[:i] + (tuple(v),) + w[i + 1:])
                if s:
                    out[k] += s * c
        return self.fin(out)

    def jet(self, p, r, s):
        for axis, num in enumerate((r, s)):
            for _ in range(num):
                p = self.deriv(p, axis)
        return p

    def fd(self, p, g):
        out = defaultdict(int)
        for w, c in p.items():
            for i, u in enumerate(w):
                if u == (g, 0, 0):
                    out[w[:i] + w[i + 1:]] += c * (-1) ** (int(g in self.odd) * sum(self.parity(x) for x in w[:i]))
        return self.fin(out)

    @lru_cache(maxsize=1 << 15)
    def q0field(self, u):
        return self.jet(self.base[u[0]], *u[1:])

    def q0(self, p):
        out = defaultdict(int)
        for w, c in p.items():
            sgn = 1
            for i, u in enumerate(w):
                img = self.q0field(u)
                if img:
                    t = self.mul({w[:i]: 1}, img, {w[i + 1:]: 1})
                    for k, d in t.items():
                        out[k] += sgn * c * d
                if self.parity(u):
                    sgn = -sgn
        return self.fin(out)

    def widiff(self, p, var):
        out = defaultdict(int)
        for w, c in p.items():
            for i, u in enumerate(w):
                if u == var:
                    out[w[:i] + w[i + 1:]] += c * (-1) ** (int(var[1] in self.odd) * sum(x[1] in self.odd for x in w[:i]))
        return self.fin(out)

    def wiedge(self, p, i, j, slot=None):
        out = {}
        for a, b in self.partners.items():
            right = (j, b) if slot is None else (j, b, slot)
            out = self.add(out, self.widiff(self.widiff(p, right), (i, a)))
        return out

    @lru_cache(maxsize=1 << 14)
    def wick(self, n, ns, sig):
        root = n + 2
        choices = [[t for t in self.vertices if self.partners[ns[sig[v]]] in t[0]] for v in range(n + 1)]
        result = defaultdict(int)
        for verts in product(*choices):
            initial = tuple((i, g) for i, (w, c) in enumerate(verts, 1) for g in w) + tuple((root, g, k) for k, g in enumerate(ns))
            coefficient = self.R(1)
            for w, c in verts:
                coefficient = coefficient * c
            p = {initial: coefficient}
            edges = []
            for i in range(1, n + 1):
                edges.extend([(i, i + 1, None), (i, root, sig[i - 1])])
            edges.append((n + 1, root, sig[n]))
            for i, j, slot in reversed(edges):
                p = self.wiedge(p, i, j, slot)
                if not p:
                    break
            for w, c in p.items():
                assert all(len(u) == 2 and u[0] < root for u in w)
                residual = tuple(tuple(g for vertex, g in w if vertex == v) for v in range(1, n + 2))
                result[residual] += c
        return tuple(self.fin(result).items())

    @lru_cache(maxsize=1 << 13)
    def selected(self, n, w):
        if n == 3 and self.degW > 4:
            raise NotImplementedError('Q3 for deg W>=5: spider/book with an active hub is not implemented')
        if n > 3 and self.degW > 3:
            raise NotImplementedError('Q_n, n >= 4, deg W >= 4: non-fan trees not implemented')
        out = {}
        ns = tuple(u[0] for u in w)
        for sig in permutations(range(n + 1)):
            coeffs = self.wick(n, ns, sig)
            jets_in = tuple(w[i][1:] for i in sig)
            moments = (graphmoment(n, jets_in) if n <= 2 else fanmoment(n, jets_in)) if coeffs else {}
            for res, c in coeffs:
                for jets, k in moments.items():
                    if any(not ww and any(jj) for ww, jj in zip(res, jets)):
                        continue
                    polys = [self.jet(self.mul(*(self.v(g) for g in ww)), *jj) for ww, jj in zip(res, jets)]
                    out = self.add(out, self.scale(self.mul(*polys), self.R((-1) ** n * 2 ** (n - 1) * k) * c))
            if n == 3 and self.degW == 4:
                coeffs = self.bookwick(ns, sig)
                if coeffs:
                    moments = bookmoment(jets_in)
                    for res, c in coeffs:
                        for jets, k in moments.items():
                            if any(not ww and any(jj) for ww, jj in zip(res, jets)):
                                continue
                            polys = [self.jet(self.mul(*(self.v(g) for g in ww)), *jj) for ww, jj in zip(res, jets)]
                            out = self.add(out, self.scale(self.mul(*polys), self.R(F(-4, 3) * self.book_sign * k) * c))
        return out

    book_sign = 1   # gate_q3_book.py flips it to test sensitivity

    @lru_cache(maxsize=1 << 13)
    def bookwick(self, ns, sig):
        """quartic_q3.bookwick for the abelian vertex list: pages 1,2,3 (cubic vertices),
        hub 4 (a W vertex, all four legs contracted), root 5 = the insertion."""
        root = 5
        edges = ((1, 4, None), (1, 5, sig[0]), (2, 4, None), (2, 5, sig[1]),
                 (3, 4, None), (3, 5, sig[2]), (4, 5, sig[3]))
        choices = []
        for v in range(4):
            choices.append([(w, c) for w, c in self.vertices
                            if self.partners[ns[sig[v]]] in w
                            and (all(g in self.gs for g in w) if v == 3 else len(w) == 3)])
        result = defaultdict(int)
        for vs in product(*choices):
            seq = tuple((v, g) for v, (w, c) in enumerate(vs, 1) for g in w) + tuple((root, g, k) for k, g in enumerate(ns))
            cc = self.R(1)
            for w, c in vs:
                cc = cc * c
            pp = {seq: cc}
            for i, j, slot in reversed(edges):
                pp = self.wiedge(pp, i, j, slot)
                if not pp:
                    break
            for w, c in pp.items():
                assert all(len(u) == 2 and u[0] < root for u in w)
                res = tuple(tuple(g for v, g in w if v == i) for i in range(1, 5))
                assert not res[3]
                result[res] += c
        return tuple(self.fin(result).items())

    def qn(self, n, p):
        if n == 0:
            return self.q0(p)
        out = defaultdict(int)
        for w, c in p.items():
            for ids in combinations(range(len(w)), n + 1):
                selected = tuple(w[i] for i in ids)
                result = self.selected(n, selected)
                if not result:
                    continue
                sg = (-1) ** sum(self.parity(w[i]) * sum(self.parity(w[j]) for j in range(i) if j not in ids[:k]) for k, i in enumerate(ids))
                rest = tuple(w[i] for i in range(len(w)) if i not in ids)
                for k_, d in self.mul(result, {rest: 1}).items():
                    out[k_] += c * sg * d
        return self.fin(out)

    @staticmethod
    def fmt(p):
        return [{'coefficient': str(c), 'monomial': list(w)} for w, c in sorted(p.items())]


# ---------------- moment polynomials (verbatim logic of the sources)
def padd(*ps):
    out = defaultdict(F)
    for p in ps:
        for a, c in p.items():
            out[a] += c
    return {a: c for a, c in out.items() if c}


def pscale(p, c):
    return {a: b * c for a, b in p.items() if b * c}


def pmul(p, q):
    out = defaultdict(F)
    for a, c in p.items():
        for b, d in q.items():
            out[tuple(x + y for x, y in zip(a, b))] += c * d
    return {a: c for a, c in out.items() if c}


@lru_cache(maxsize=1 << 12)
def graphmoment(n, jets):
    """engine.py graphmoment, n in (1,2)."""
    assert n in (1, 2)
    N = 4 + 2 * (n + 1)
    one = {(0,) * N: F(1)}

    def var(i):
        a = [0] * N
        a[i] = 1
        return {tuple(a): F(1)}
    a, b, u, v = [var(i) for i in range(4)]
    lam = [[var(4 + 2 * i + j) for j in range(2)] for i in range(n + 1)]

    def minus(p):
        return padd(one, pscale(p, -1))

    def wedge(i, j):
        return padd(pmul(lam[i][0], lam[j][1]), pscale(pmul(lam[i][1], lam[j][0]), -1))
    if n == 1:
        density = wedge(0, 1)
        currents = [[padd(pmul(minus(u), lam[0][j]), pmul(v, lam[1][j])) for j in range(2)],
                    [padd(pmul(u, lam[0][j]), pmul(minus(v), lam[1][j])) for j in range(2)]]
    else:
        density = pmul(b, pmul(padd(pmul(minus(u), wedge(0, 1)), pmul(v, wedge(0, 2))), padd(pmul(a, wedge(0, 2)), wedge(1, 2))))
        y23 = [padd(pmul(u, padd(pmul(a, lam[0][j]), lam[1][j])), pscale(pmul(v, lam[2][j]), -1)) for j in range(2)]
        y12 = [padd(pmul(pmul(a, b), lam[0][j]), pscale(pmul(minus(b), lam[1][j]), -1), pmul(minus(b), y23[j])) for j in range(2)]
        currents = [[padd(lam[0][j], pscale(y12[j], -1)) for j in range(2)],
                    [padd(lam[1][j], y12[j], pscale(y23[j], -1)) for j in range(2)],
                    [padd(lam[2][j], y23[j]) for j in range(2)]]
    p = density
    for i, jet in enumerate(jets):
        for j, num in enumerate(jet):
            for _ in range(num):
                p = pmul(p, currents[i][j])
    out = defaultdict(F)
    for ex, c in p.items():
        A, B, U, V = ex[:4]
        mom = F(factorial(U) * factorial(V), (A + 1) * (B + 1) * factorial(U + V + 2))
        jj = tuple(tuple(ex[4 + 2 * i:6 + 2 * i]) for i in range(n + 1))
        out[jj] += c * mom
    return {j: c for j, c in out.items() if c}


@lru_cache(maxsize=1 << 12)
def fanmoment(n, jets, active=None):
    """higher_Q_article_2026-09-20/higher_moments.py fanmoment (all n)."""
    nv = n + 1
    np_ = 2 * n
    N = np_ + 2 * nv
    one = {(0,) * N: F(1)}

    def var(i):
        e = [0] * N
        e[i] = 1
        return {tuple(e): F(1)}

    def minus(p):
        return padd(one, pscale(p, -1))

    def wedge(p, q):
        return padd(pmul(p[0], q[1]), pscale(pmul(p[1], q[0]), -1))

    def vsum(*vs):
        return [padd(*(v[a] for v in vs)) for a in range(2)]

    def vmul(p, v):
        return [pmul(p, x) for x in v]
    aa = [var(2 * i) for i in range(n - 1)]
    bb = [var(2 * i + 1) for i in range(n - 1)]
    u, v = var(np_ - 2), var(np_ - 1)
    lam = [[var(np_ + 2 * i + a) if active is None or i in active else {} for a in range(2)] for i in range(nv)]
    effective = [lam[0]]
    for i in range(n - 1):
        effective.append(vsum(lam[i + 1], vmul(aa[i], effective[-1])))
    currents = [None] * nv
    currents[n - 1] = vsum(vmul(minus(u), effective[n - 1]), vmul(v, lam[n]))
    currents[n] = vsum(vmul(u, effective[n - 1]), vmul(minus(v), lam[n]))
    density = wedge(effective[n - 1], lam[n])
    for i in reversed(range(n - 1)):
        density = pmul(density, pmul(bb[i], wedge(effective[i], currents[i + 1])))
        currents[i] = vsum(vmul(minus(aa[i]), effective[i]), vmul(minus(bb[i]), currents[i + 1]))
        currents[i + 1] = vmul(bb[i], currents[i + 1])
    for i, jet in enumerate(jets):
        for axis, count in enumerate(jet):
            for _ in range(count):
                density = pmul(density, currents[i][axis])
    out = defaultdict(F)
    for ex, c in density.items():
        value = c
        for power in ex[:np_ - 2]:
            value /= power + 1
        U, V = ex[np_ - 2:np_]
        value *= F(factorial(U) * factorial(V), factorial(U + V + 2))
        jj = tuple(tuple(ex[np_ + 2 * i:np_ + 2 * i + 2]) for i in range(nv))
        out[jj] += value
    return {j: c for j, c in out.items() if c}


@lru_cache(maxsize=1 << 10)
def bookmoment(jets):
    """higher_Q_article_2026-09-20/higher_moments.py bookmoment (verbatim logic)."""
    n = 3; nv = 4; np_ = 6; N = 12
    one = {(0,) * N: F(1)}

    def var(i):
        e = [0] * N
        e[i] = 1
        return {tuple(e): F(1)}

    def minus(p):
        return padd(one, pscale(p, -1))
    aa = [var(2 * i) for i in range(3)]
    bb = [var(2 * i + 1) for i in range(3)]
    lam = [[var(np_ + 2 * i + a) for a in range(2)] for i in range(3)]
    common = [padd(*(pmul(aa[i], lam[i][axis]) for i in range(3))) for axis in range(2)]
    complement = padd(one, *(pscale(p, -1) for p in bb))
    currents = [[padd(pmul(minus(aa[i]), lam[i][a]), pmul(bb[i], common[a])) for a in range(2)] for i in range(3)]
    currents.append([pmul(complement, p) for p in common])
    density = dict(one)
    for i in range(3):
        density = pmul(density, padd(pmul(lam[i][0], common[1]), pscale(pmul(lam[i][1], common[0]), -1)))
    for i, jet in enumerate(jets):
        for a, count in enumerate(jet):
            for _ in range(count):
                density = pmul(density, currents[i][a])
    out = defaultdict(F)
    for ex, c in density.items():
        value = c
        for i in range(3):
            value *= F(factorial(ex[2 * i + 1]), ex[2 * i] + 1)
        value /= factorial(4 + ex[1] + ex[3] + ex[5])
        jj = tuple(tuple(ex[np_ + 2 * i:np_ + 2 * i + 2]) for i in range(3)) + ((0, 0),)
        out[jj] += value
    return {j: c for j, c in out.items() if c}
