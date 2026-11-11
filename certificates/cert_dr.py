"""Exact certificate for a d_r lift on the SU(2)-invariant subcomplex.

Certificate = chain (x_0..x_{r-1}) and a functional phi = (phi_1..phi_r) on the H = 0 rows of the stage blocks:
  (V1) x_k invariant (H = 0, L x_k = 0) in the stage-k block,   sum_{a+b=k} Q_a x_b = 0  (k = 0..r-1)
  (V2) for every H = 0 monomial m of the stage-m block (m = 1..r):  sum_{k=m}^{r} phi_k(Q_{k-m} m) = 0
  (V3) phi_r(t) != 0,   t = sum_{k=1}^{r} Q_k x_{r-k}
(V2) is checked on all monomials, so t is outside the span of the stacked images of the full (hence also of the
invariant) stage blocks: d_r[x_0] != 0 in the spectral sequence of the invariant subcomplex and of the full complex.

  PRIME=p python3 cert_dr.py prime THEORY 'KEY' R OUT.json.gz   chain + phi mod p (phi = phi0 o Reynolds)
  python3 cert_dr.py crt OUT.json.gz IN1.json.gz IN2.json.gz ...   CRT + rational reconstruction
  PRIME=0 python3 cert_dr.py verify THEORY 'KEY' R CERT.json.gz    exact checks (V1)-(V3)
"""
import gzip, json, os, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nonab import lin_, rep_matrices, SU2
from elim import Elim, EXE
from hunt2 import combo, tag
from hunt3 import stackcol


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


def wkey(w):
    return json.dumps([list(u) for u in w], separators=(',', ':'))


def unkey(s):
    return tuple(tuple(u) for u in json.loads(s))


def maps(T, name):
    """lowering L (= the raise_map of nonab, H-charge -2) and raising U (from rho(F)) on letters"""
    from nonab_hunt import THEORIES
    th = THEORIES[name]
    Lm, Um = {}, {}
    for nm, (kind, w) in th['matter'].items():
        comps, rho, ch = rep_matrices(kind)
        for X, M in (('E', Lm), ('F', Um)):
            for (i, j), v in rho[X].items():
                M.setdefault(nm + i, []).append((nm + j, v))
                M.setdefault('B' + nm + j, []).append(('B' + nm + i, -v))
    _, rad, _ = rep_matrices('adj')
    for X, M in (('E', Lm), ('F', Um)):
        for (C, B_), v in rad[X].items():
            M.setdefault('c' + C, []).append(('c' + B_, v))
            M.setdefault('b' + B_, []).append(('b' + C, -v))
    return Lm, Um


def act(T, rm, p):
    out = {}
    for w, c in p.items():
        for i, u in enumerate(w):
            for g2, a in rm.get(u[0], ()):
                k, s = T.canon(w[:i] + ((g2, u[1], u[2]),) + w[i + 1:])
                if s:
                    out[k] = out.get(k, 0) + s * c * T.R(a)
    return T.fin(out)


def blocks(T, key, r):
    D = T.D
    W, R, S, gh = key[:4]
    f = tuple(key[4:])
    st = lambda m, dgh=0: (W - m * D, R + m, S + m, (gh + 2 * m + dgh, 0) + f)
    return st


# ---------------------------------------------------------------- prime mode
def prime_mode(name, key, r, out):
    from nonab_hunt import build, P
    T, G = build(name)
    el = Elim(P, 'c' if os.path.exists(EXE) else 'py')
    st = blocks(T, key, r)

    def IB(m, dgh=0):
        W_, R_, S_, lin = st(m, dgh)
        return G.inv_basis(W_, R_, S_, lin, el)
    tic = time.time()
    src, prev = IB(0), IB(0, -1)
    _, res = el.run([], [T.q0(o) for o in src], combos=True)
    cyc = [lin_(T, src, combo(None, co, 't')) for co in res if co is not None]
    _, res = el.run([T.q0(o) for o in prev], cyc, combos=False)
    classes = [o for o, r_ in zip(cyc, res) if r_ is None]
    log('E1', len(classes))
    chains = [[o] for o in classes]
    prims, imgs = {}, {}
    for s in range(1, r + 1):
        prims[s] = IB(s)
        log('stage', s, 'inv dim', len(prims[s]))
        for m in range(1, s + 1):
            if (m, s - m) not in imgs:
                imgs[(m, s - m)] = [T.qn(s - m, x) if s - m else T.q0(x) for x in prims[m]]
        tests = [T.add(*(T.qn(k, ch[s - k]) for k in range(1, s + 1))) for ch in chains]
        base, owner = [], []
        for m in range(1, s + 1):
            for j in range(len(prims[m])):
                base.append(stackcol([(k, imgs[(m, k - m)][j]) for k in range(m, s + 1)]))
                owner.append((m, j))
        rk, res = el.run(base, [tag(s, t) for t in tests], combos=True)
        log('d%d rank' % s, sum(co is None for co in res))
        if s == r:
            break
        new = []
        for co in res:
            if co is None:
                continue
            tc, bc = combo(None, co, 't'), combo(None, co, 'b')
            z = {m: [] for m in range(1, s + 1)}
            for b, c in bc.items():
                m, j = owner[b]
                z[m].append(T.scale(prims[m][j], c))
            ch = [T.add(*(T.scale(chains[i][0], c) for i, c in tc.items()))]
            for m in range(1, s):
                ch.append(T.add(*(T.scale(chains[i][m], c) for i, c in tc.items()), *z[m]))
            ch.append(T.add(*z[s]))
            new.append(ch)
        chains = new
    hit = [i for i, co in enumerate(res) if co is None]
    assert hit, 'no d_r lift here'
    i0 = hit[0]
    chain, t = chains[i0], tag(r, tests[i0])
    # phi0: phi0 . base = 0, phi0 . t = 1  (transposed system, one column per row key)
    rows = {}
    for j, col in enumerate(base):
        for k, c in col.items():
            rows.setdefault(k, {})[j] = c
    for k, c in t.items():
        rows.setdefault(k, {})['t'] = c
    rkeys = list(rows)
    tcols = [{('t' if j == 't' else j): c for j, c in rows[k].items()} for k in rkeys]
    log('transposed system', len(tcols), 'columns')
    _, res = el.run(tcols, [{'t': 1}], combos=True)
    assert res[0] is not None, 'phi0 does not exist'
    phi0 = {}
    for (kind, i), c in res[0].items():
        if kind == 'b' and c % P:
            phi0[rkeys[i]] = (-c) % P
    # check phi0 on the invariant base and t
    for col in base:
        assert sum(phi0.get(k, 0) * c for k, c in col.items()) % P == 0
    assert sum(phi0.get(k, 0) * c for k, c in t.items()) % P == 1
    log('phi0 terms', len(phi0))
    # Reynolds: phi = phi0 o prod_j (1 - U L / (j(j+1))) on each stage row block (H = 0)
    Lm, Um = maps(T, name)
    phi = {}
    for k in range(1, r + 1):
        W_, R_, S_, lin = st(k, 1)
        V0 = G.basis(W_, R_, S_, lin)
        Vm = G.basis(W_, R_, S_, (lin[0], -2) + lin[2:])
        Lv = [act(T, Lm, {w: 1}) for w in V0]
        Uu = [act(T, Um, {u: 1}) for u in Vm]
        cur = {w: c for (kk, w), c in phi0.items() if kk == k}
        j = 0
        while True:
            psi = {}
            for u, im in zip(Vm, Uu):
                s_ = sum(cur.get(w, 0) * c for w, c in im.items()) % P
                if s_:
                    psi[u] = s_
            if not psi:
                break
            j += 1
            inv = pow(j * (j + 1), P - 2, P)
            new = {}
            for w, im in zip(V0, Lv):
                s_ = (cur.get(w, 0) - inv * sum(psi.get(u, 0) * c for u, c in im.items())) % P
                if s_:
                    new[w] = s_
            cur = new
            assert j < 40, 'Reynolds did not terminate'
        log('stage', k, 'rows', len(V0), 'spins', j, 'phi terms', len(cur))
        for w, c in cur.items():
            phi[(k, w)] = c
    assert sum(phi.get(k, 0) * c for k, c in t.items()) % P == 1
    for col in base:
        assert sum(phi.get(k, 0) * c for k, c in col.items()) % P == 0
    rec = {'P': P, 'theory': name, 'key': key, 'r': r, 'E1': len(classes), 'seconds': round(time.time() - tic, 1),
           'x': [{wkey(w): c for w, c in x.items()} for x in chain],
           'phi': {str(k): {wkey(w): c for (kk, w), c in phi.items() if kk == k} for k in range(1, r + 1)}}
    with gzip.open(out, 'wt') as fh:
        json.dump(rec, fh)
    log('wrote', out, rec['seconds'], 's')


# ---------------------------------------------------------------- CRT
def isqrt(n):
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def ratrec(a, m):
    """Wang rational reconstruction of a mod m, |num|, den <= sqrt(m/2)"""
    bound = isqrt(m // 2)
    r0, r1, s0, s1 = m, a % m, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return F(r1, s1)


def crt_mode(out, files):
    recs = []
    for fn in files:
        with gzip.open(fn, 'rt') as fh:
            recs.append(json.load(fh))
    Ps = [x['P'] for x in recs]
    assert len(set(Ps)) == len(Ps)
    M = 1
    for p in Ps:
        M *= p

    def rec_entry(vals):
        a, m = 0, 1
        for p, v in zip(Ps, vals):
            # a = a mod m, v mod p
            t = ((v - a) * pow(m % p, p - 2, p)) % p
            a, m = a + m * t, m * p
        return ratrec(a, m)

    def rec_dict(ds):
        keys = set()
        for d in ds:
            keys |= set(d)
        out, bad, hmax = {}, 0, 0
        for k in keys:
            q = rec_entry([d.get(k, 0) for d in ds])
            if q is None:
                bad += 1
                continue
            if q:
                out[k] = [q.numerator, q.denominator]
                hmax = max(hmax, abs(q.numerator).bit_length() + q.denominator.bit_length())
        return out, bad, hmax
    r = recs[0]['r']
    cert = {'theory': recs[0]['theory'], 'key': recs[0]['key'], 'r': r, 'primes': Ps, 'x': [], 'phi': {}}
    stats = []
    for i in range(r):
        d, bad, h = rec_dict([x['x'][i] for x in recs])
        cert['x'].append(d)
        stats.append(('x%d' % i, len(d), bad, h))
    for k in map(str, range(1, r + 1)):
        d, bad, h = rec_dict([x['phi'][k] for x in recs])
        cert['phi'][k] = d
        stats.append(('phi%s' % k, len(d), bad, h))
    print(json.dumps({'primes': len(Ps), 'log2M': M.bit_length(), 'entries': stats}))
    with gzip.open(out, 'wt') as fh:
        json.dump(cert, fh)


# ---------------------------------------------------------------- exact verification
def verify_mode(name, key, r, certf):
    from nonab_hunt import build, P
    assert P == 0, 'run with PRIME=0'
    T, G = build(name)
    st = blocks(T, key, r)
    with gzip.open(certf, 'rt') as fh:
        cert = json.load(fh)
    x = [{unkey(w): F(a, b) for w, (a, b) in d.items()} for d in cert['x']]
    phi = {k: {unkey(w): F(a, b) for w, (a, b) in cert['phi'][str(k)].items()} for k in range(1, r + 1)}
    Lm, Um = maps(T, name)
    rep = {'theory': name, 'key': key, 'r': r, 'primes': len(cert['primes'])}
    # (V1)
    for k in range(r):
        W_, R_, S_, lin = st(k)
        assert x[k], 'empty chain element'
        for w in x[k]:
            assert G.key_of(w) == (W_, R_, S_, lin), ('block', k, w)
        assert not act(T, Lm, x[k]), ('not invariant', k)
        eq = T.add(T.q0(x[k]), *(T.qn(a, x[k - a]) for a in range(1, k + 1)))
        assert not eq, ('chain identity', k)
    nq = T.add(*(T.qn(a, x[r - a]) for a in range(1, r + 1)))
    val = sum(phi[r].get(w, 0) * c for w, c in nq.items())
    rep['phi_r(t)'] = str(val)
    assert val != 0
    log('V1, V3 pass; phi_r(t) =', val)
    # (V2) on all H = 0 monomials of the stage blocks
    for m in range(1, r + 1):
        W_, R_, S_, lin = st(m)
        mons = G.basis(W_, R_, S_, lin)
        tic = time.time()
        for idx, w in enumerate(mons):
            p = {w: F(1)}
            tot = 0
            for k in range(m, r + 1):
                im = T.qn(k - m, p) if k > m else T.q0(p)
                ph = phi[k]
                tot += sum(ph.get(u, 0) * c for u, c in im.items())
            assert tot == 0, ('V2 fails', m, w)
            if idx % 5000 == 0:
                log('V2 stage', m, idx, '/', len(mons), round(time.time() - tic), 's')
        rep['V2_stage%d_monomials' % m] = len(mons)
    rep['verified'] = True
    print(json.dumps(rep))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'prime':
        prime_mode(sys.argv[2], json.loads(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    elif mode == 'crt':
        crt_mode(sys.argv[2], sys.argv[3:])
    elif mode == 'verify':
        verify_mode(sys.argv[2], json.loads(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
