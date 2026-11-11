"""Exact rational (PRIME=0) or modular (PRIME=p) verifier for abelian d_r lift certificates.

Certificate format (identical to cert_dr.py):
  chain x = (x_0, ..., x_{r-1}) and dual functional phi = (phi_1, ..., phi_r) with rational entries [num, den]:
  (V1) Every monomial of x_k lies in the stage-k block B(k), and sum_{a+b=k} Q_a x_b = 0 for k = 0..r-1.
  (V2) For every monomial y in B(m) (m = 1..r): sum_{k=m}^r phi_k(Q_{k-m} y) = 0.
  (V3) t = sum_{k=1}^r Q_k x_{r-k} satisfies Q_0 t = 0 and phi_r(t) != 0.

Usage:
  PRIME=0 python3 verify_abelian_cert.py FAMILY 'KEY' R CERT.json.gz
"""
import gc, gzip, json, os, resource, sys, time
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
for cand in (HERE, HERE.parent):
    if (cand / 'qengine.py').exists():
        sys.path.insert(0, str(cand))
        break

from qengine import Theory, graphmoment, fanmoment, bookmoment

BUILTIN_FAMILIES = {
    'singlet': {
        'charges': {'x': (1,), 'y': (-1,), 's': (0,)},
        'W': {('s', 'x', 'y'): F(1)},
        'weights': {'x': 1, 'y': 1, 's': 1},
    },
    'singlet_cubic_self': {
        'charges': {'x': (1,), 'y': (-1,), 's': (0,)},
        'W': {('s', 'x', 'y'): F(1), ('s', 's', 's'): F(1, 3)},
        'weights': {'x': 1, 'y': 1, 's': 1},
    },
}


def _nullspace_int(M, n):
    from math import gcd
    lcm = lambda a, b: a * b // gcd(a, b)
    rows = [list(map(F, r)) for r in M]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        a = rows[r][c]
        rows[r] = [x / a for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                b = rows[i][c]
                rows[i] = [x - b * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    basis = []
    for fc in [c for c in range(n) if c not in piv]:
        v = [F(0)] * n
        v[fc] = F(1)
        for i, pc in enumerate(piv):
            v[pc] = -rows[i][fc]
        L = 1
        for x in v:
            L = lcm(L, x.denominator)
        basis.append(tuple(int(x * L) for x in v))
    return basis


class Graded:
    """Self-contained block grading (identical to hunt.Graded, without external imports)."""
    def __init__(self, name, P):
        if name in BUILTIN_FAMILIES:
            f = BUILTIN_FAMILIES[name]
        else:
            import families as FAM
            f = FAM.get(name)
        self.name = name
        self.T = T = Theory(f['charges'], f['W'], P=P, weights=f['weights'])
        gs = list(T.gs)
        M = []
        for m in f['W']:
            row = [0] * len(gs)
            for g in m:
                row[gs.index(g)] += 1
            M.append(row)
        flav = _nullspace_int(M, len(gs)) if M else [tuple(int(i == j) for j in range(len(gs))) for i in range(len(gs))]
        self.flav = flav
        D = T.D
        self.wt, self.gh, self.q, self.f = {}, {}, {}, {}
        zq = (0,) * T.rank
        zf = (0,) * len(flav)
        for i, g in enumerate(gs):
            self.wt[g], self.gh[g], self.q[g] = T.wt[g], 0, tuple(T.charges[g])
            self.f[g] = tuple(v[i] for v in flav)
            B = T.bs[g]
            self.wt[B], self.gh[B], self.q[B] = D - T.wt[g], -1, tuple(-x for x in T.charges[g])
            self.f[B] = tuple(-v[i] for v in flav)
        for a in range(T.rank):
            self.wt[f'b{a}'], self.gh[f'b{a}'], self.q[f'b{a}'], self.f[f'b{a}'] = D, -2, zq, zf
            self.wt[f'c{a}'], self.gh[f'c{a}'], self.q[f'c{a}'], self.f[f'c{a}'] = 0, 1, zq, zf
        self.lin = {g: (self.gh[g],) + self.q[g] + self.f[g] for g in T.fields}
        nonc = [g for g in T.fields if not g.startswith('c')]
        k = len(self.lin[T.fields[0]])
        self.amax = [max(abs(self.lin[g][j]) / self.wt[g] for g in nonc) for j in range(k)]
        self.cmax = [abs(self.lin['c0'][j]) for j in range(k)]

    def basis(self, W, R, S, lin=None):
        T = self.T
        if min(W, R, S) < 0:
            return []
        vs = tuple(sorted((g, r, s) for g in T.fields for r in range(R + 1) for s in range(S + 1) if not (g.startswith('c') and r + s == 0)))
        L = [self.lin[u[0]] for u in vs]
        WT = [self.wt[u[0]] for u in vs]
        PAR = [T.parity(u) for u in vs]
        k = len(self.amax)
        out = []
        amax, cmax = self.amax, self.cmax

        def rec(i, w, r_, s_, cur, word):
            if w == 0 and r_ == 0 and s_ == 0:
                if lin is None or tuple(cur) == lin:
                    out.append(word)
                return
            if lin is not None:
                jr = r_ + s_
                for j in range(k):
                    if abs(lin[j] - cur[j]) > amax[j] * w + cmax[j] * jr + 1e-9:
                        return
            for j in range(i, len(vs)):
                u = vs[j]
                if WT[j] > w or u[1] > r_ or u[2] > s_:
                    continue
                nl = L[j]
                for t in range(k):
                    cur[t] += nl[t]
                rec(j + PAR[j], w - WT[j], r_ - u[1], s_ - u[2], cur, word + (u,))
                for t in range(k):
                    cur[t] -= nl[t]

        rec(0, W, R, S, [0] * k, ())
        return out


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


def unkey(s):
    return tuple(tuple(u) for u in json.loads(s))


def make_fast_ops(T):
    odd = frozenset(T.odd)
    one = T.R(1)

    def canon(w):
        os_ = [u for u in w if u[0] in odd]
        n = len(os_)
        sg = 1
        if n > 1:
            if len(set(os_)) != n:
                return (), 0
            inv = 0
            for i in range(n - 1):
                a = os_[i]
                for j in range(i + 1, n):
                    if a > os_[j]:
                        inv += 1
            if inv & 1:
                sg = -1
        return tuple(sorted(w)), sg

    def mul2(p, q):
        if not p or not q:
            return {}
        out = defaultdict(int)
        for w1, c1 in p.items():
            for w2, c2 in q.items():
                k, s = canon(w1 + w2)
                if s:
                    out[k] += (c1 * c2 if s == 1 else -c1 * c2)
        return T.fin(out)

    @lru_cache(maxsize=2048)
    def vjet(ww, jj):
        if not ww:
            return {(): one} if jj == (0, 0) else {}
        base = {(): one}
        for g in ww:
            base = mul2(base, T.v(g))
        return T.jet(base, *jj)

    @lru_cache(maxsize=2048)
    def fanmoment_cached(n, jets_in, active):
        return fanmoment(n, jets_in, active)

    @lru_cache(maxsize=4096)
    def selected_fast(n, w):
        if n == 3 and T.degW > 4:
            raise NotImplementedError
        if n > 3 and T.degW > 3:
            raise NotImplementedError
        acc = defaultdict(int)
        ns = tuple(u[0] for u in w)
        pref = T.R((-1) ** n * 2 ** (n - 1))
        for sig in permutations(range(n + 1)):
            coeffs = T.wick(n, ns, sig)
            jets_in = tuple(w[i][1:] for i in sig)
            if coeffs:
                if n <= 2:
                    moments = graphmoment(n, jets_in)
                    for res, c in coeffs:
                        act_idx = [i for i, ww in enumerate(res) if ww]
                        for jets, k in moments.items():
                            if any(not ww and (jj[0] or jj[1]) for ww, jj in zip(res, jets)):
                                continue
                            poly = {(): one}
                            for i in act_idx:
                                poly = mul2(poly, vjet(res[i], jets[i]))
                                if not poly:
                                    break
                            if poly:
                                sc = pref * T.R(k) * c
                                for m, d in poly.items():
                                    acc[m] += d * sc
                else:
                    for res, c in coeffs:
                        active = tuple(i for i, ww in enumerate(res) if ww)
                        moments = fanmoment_cached(n, jets_in, active)
                        for jets, k in moments.items():
                            poly = {(): one}
                            for i in active:
                                poly = mul2(poly, vjet(res[i], jets[i]))
                                if not poly:
                                    break
                            if poly:
                                sc = pref * T.R(k) * c
                                for m, d in poly.items():
                                    acc[m] += d * sc
            if n == 3 and T.degW == 4:
                bcoeffs = T.bookwick(ns, sig)
                if bcoeffs:
                    moments = bookmoment(jets_in)
                    for res, c in bcoeffs:
                        act_idx = [i for i, ww in enumerate(res) if ww]
                        for jets, k in moments.items():
                            if any(not ww and (jj[0] or jj[1]) for ww, jj in zip(res, jets)):
                                continue
                            poly = {(): one}
                            for i in act_idx:
                                poly = mul2(poly, vjet(res[i], jets[i]))
                                if not poly:
                                    break
                            if poly:
                                sc = T.R(F(-4, 3) * T.book_sign * k) * c
                                for m, d in poly.items():
                                    acc[m] += d * sc
        return T.fin(acc)

    def q0_fast(p):
        out = defaultdict(int)
        for w, c in p.items():
            sc = c
            for i, u in enumerate(w):
                img = T.q0field(u)
                if img:
                    a, b = w[:i], w[i + 1:]
                    for m, d in img.items():
                        k, s = canon(a + m + b)
                        if s:
                            out[k] += (sc * d if s == 1 else -sc * d)
                if u[0] in odd:
                    sc = -sc
        return T.fin(out)

    def qn_fast(n, p):
        if n == 0:
            return q0_fast(p)
        out = defaultdict(int)
        for w, c in p.items():
            lw = len(w)
            pars = [1 if u[0] in odd else 0 for u in w]
            for ids in combinations(range(lw), n + 1):
                sel = tuple(w[i] for i in ids)
                result = selected_fast(n, sel)
                if not result:
                    continue
                idset = set(ids)
                sg_exp = sum(pars[i] * sum(pars[j] for j in range(i) if j not in ids[:k]) for k, i in enumerate(ids))
                sg = -c if (sg_exp & 1) else c
                rest = tuple(w[i] for i in range(lw) if i not in idset)
                for m, d in result.items():
                    k_, s = canon(m + rest)
                    if s:
                        out[k_] += (sg if s == 1 else -sg) * d
        return T.fin(out)

    def clear_caches():
        selected_fast.cache_clear()
        fanmoment_cached.cache_clear()
        vjet.cache_clear()
        gc.collect()

    return q0_fast, qn_fast, clear_caches


def _rss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 2) if sys.platform == 'darwin' else r / 1024


def verify(fam, key, r, certf):
    P = int(os.environ.get('PRIME', '0'))
    max_rss_mb = float(os.environ.get('MAX_RSS_MB', '2048'))
    t0 = time.time()
    G = Graded(fam, P)
    T = G.T
    q0f, qnf, clear_caches = make_fast_ops(T)
    D = T.D
    W, R, S, gh, q, f = key
    q, f = tuple(q), tuple(f)
    one = T.R(1)

    with gzip.open(certf, 'rt') as fh:
        cert = json.load(fh)
    assert cert['r'] == r and cert['key'] == key

    def cvt(a, b):
        return T.R(F(a, b))

    xs = [{unkey(w): cvt(a, b) for w, (a, b) in d.items() if cvt(a, b)} for d in cert['x']]
    phis = {k: {unkey(w): cvt(a, b) for w, (a, b) in cert['phi'][str(k)].items() if cvt(a, b)} for k in range(1, r + 1)}

    def B(m):
        return G.basis(W - m * D, R + m, S + m, (gh + 2 * m,) + q + f)

    # (V1) Grading check and chain equations sum_{a+b=k} Q_a x_b == 0 for k = 0..r-1
    for k in range(r):
        bset = set(B(k))
        assert xs[k], f'empty chain element x_{k}'
        for w in xs[k]:
            assert w in bset, f'monomial {w} outside stage-{k} basis'
        eq = T.add(*(qnf(a, xs[k - a]) if a else q0f(xs[k]) for a in range(k + 1)))
        assert not eq, f'V1 failed at k={k}'
    log(f'V1 pass (chain equations k=0..{r-1})')

    # (V3) Target t = sum_{k=1}^r Q_k x_{r-k} is Q0-closed and phi_r(t) != 0
    t_vec = T.add(*(qnf(k, xs[r - k]) for k in range(1, r + 1)))
    assert not q0f(t_vec), 'target t not Q0-closed'
    val = sum(phis[r].get(w, 0) * c for w, c in t_vec.items())
    if P:
        val %= P
    assert val != 0, 'V3 failed: phi_r(t) == 0'
    log(f'V3 pass; phi_{r}(t) = {val}')

    # (V2) Orthogonality on every monomial y of stage m = 1..r
    rep = {'family': fam, 'key': key, 'r': r, 'P': P, 'primes': len(cert.get('primes', [])), f'phi_{r}(t)': str(val)}
    for m in range(1, r + 1):
        mons = B(m)
        tm = time.time()
        for idx, w in enumerate(mons):
            ym = {w: one}
            tot = 0
            for k in range(m, r + 1):
                pk = phis[k]
                if not pk:
                    continue
                im = qnf(k - m, ym) if k > m else q0f(ym)
                for u, c in im.items():
                    if u in pk:
                        tot += pk[u] * c
            if P:
                tot %= P
            assert tot == 0, f'V2 failed at stage m={m}, monomial {idx}'
            if idx and idx % 5000 == 0:
                if _rss_mb() > max_rss_mb:
                    clear_caches()
                    assert _rss_mb() <= max_rss_mb, f'RSS exceeded {max_rss_mb} MB guard'
                if idx % 10000 == 0:
                    log(f'V2 stage {m}: {idx}/{len(mons)} ({time.time()-tm:.1f}s, rss={_rss_mb():.0f}MB)')
        rep[f'V2_stage{m}_monomials'] = len(mons)
        log(f'V2 stage {m} pass ({len(mons)} monomials, {time.time()-tm:.1f}s, rss={_rss_mb():.0f}MB)')
        clear_caches()

    rep['verified'] = True
    rep['seconds'] = round(time.time() - t0, 2)
    print(json.dumps(rep))


if __name__ == '__main__':
    if len(sys.argv) == 6 and sys.argv[1] == 'verify':
        _, _, fam, keys, rs, certf = sys.argv
    else:
        _, fam, keys, rs, certf = sys.argv
    verify(fam, json.loads(keys), int(rs), certf)
