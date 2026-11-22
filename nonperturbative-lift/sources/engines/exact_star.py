"""Quantum product A * B of local operators in components, from the product-forest graphs.

    [hbar^n] (A * B) = KAPPA^n * sum over forest topologies (1 / |Aut|) * sum over assignments of letters of A, B to
                       the external edges of [ raw-cut kernel moment ] x [ ordered fibre contractions ]

Graphs: n cubic interaction vertices forming a forest of full binary trees, leaves on the inputs A, B
(census.py); one surviving letter at the root of every tree; kernels from product_kernels.py
(tree charts, exact, arbitrary jets).  Sources: lambda_root -> derivative on the surviving letter, lambda_A -> total
derivative on the uncontracted letters of A; B is the base point.  A jet d^m on a contracted input letter
= (-1)^|m| m! x Taylor coefficient of the kernel in the displacement of its edge.
Contractions: engine convention (qengine.Theory.wiedge), edges in lexicographic order, applied right to left to
the product I_1 ... I_n A B.  Normalisation: engine (stored) units, Q_n(physical) = (64 pi^2)^(-n) Q_n(stored).
"""
import sys, itertools, time, random
from math import factorial
from functools import lru_cache
import exact_product as E
from sunlib import *

import census as CEN
import product_kernels as PK

KAPPA = [T.R(-2)]         # fixed at one loop by the Leibniz identity (t1.py); the same constant at every order


def use(theory):
    """switch to another ITheory (cubic interaction, same prime)"""
    global T
    T = theory
    wick.cache_clear()


def fin(d): return T.fin(d)
def add(*ps): return T.add(*ps)
def sc(p, a): a = T.R(a); return fin({w: c * a for w, c in p.items()})
def sub(p, q): return fin({**p, **{w: p.get(w, 0) - c for w, c in q.items()}})
def mul(*ps): return T.mul(*ps)


@lru_cache(maxsize=None)
def topologies(n):
    out = []
    for f in CEN.forests(n):
        edges, roots = CEN.edges_of(f)
        out.append((tuple(map(tuple, edges)), tuple(roots), CEN.aut_forest(f), tuple(map(CEN.key, f))))
    return tuple(out)


@lru_cache(maxsize=None)
def _kern(n, ti):
    edges, roots, aut, name = topologies(n)[ti]
    syms = [S.symbols('l%d_0 l%d_1' % (r, r)) for r in roots] + [S.symbols('lA_0 lA_1')]
    mom = {r: syms[i] for i, r in enumerate(roots)}
    mom[n + 1] = syms[-1]
    return SourceFactoredKernel(n, list(edges), mom), [x for s in syms for x in s]


@lru_cache(maxsize=None)
def moments_direct(n, ti, jets):
    """the same moments from the chart of the whole forest (slow; cross-check of product_kernels)"""
    edges, roots, aut, name = topologies(n)[ti]
    k, syms = _kern(n, ti)
    ext = [e for e in edges[:-1] if e[1] > n]
    powers, fac = {}, 1
    for e, (r, s) in zip(ext, jets):
        if r: powers[(e, 0)] = r
        if s: powers[(e, 1)] = s
        fac *= (-1) ** (r + s) * factorial(r) * factorial(s)
    c = S.expand(k.coefficient(powers) * fac)
    out = {}
    if c != 0:
        for exps, coef in S.Poly(c, *syms).terms():
            jj = tuple((exps[2 * i], exps[2 * i + 1]) for i in range(len(roots) + 1))
            out[jj] = F(int(coef.p), int(coef.q))
    return out


def moments_exact(n, ti, jets):
    """jets: (r, s) for every external edge, in edge order.  Returns {(jj_root_1, ..., jj_A): Fraction}."""
    assert PK.forests(n) == CEN.forests(n)
    return PK.moments(n, ti, jets)


@lru_cache(maxsize=None)
def moments(n, ti, jets):
    return tuple((jj, T.R(c)) for jj, c in moments_exact(n, ti, jets).items())


@lru_cache(maxsize=1 << 16)
def wick(n, ti, namesA, namesB, aA, aB):
    edges, roots, aut, name = topologies(n)[ti]
    A, B = n + 1, n + 2
    phys = edges[:-1]
    extA = [e for e in phys if e[1] == A]; extB = [e for e in phys if e[1] == B]
    slot = {}
    for e, k in zip(extA, aA): slot[e] = k
    for e, k in zip(extB, aB): slot[e] = k
    need = {v: [] for v in range(1, n + 1)}
    for e in extA: need[e[0]].append(T.partners[namesA[slot[e]]])
    for e in extB: need[e[0]].append(T.partners[namesB[slot[e]]])
    choices = []
    for v in range(1, n + 1):
        ch = []
        for w, c in T.vertices:
            l = list(w); ok = True
            for g in need[v]:
                if g in l: l.remove(g)
                else: ok = False; break
            if ok: ch.append((w, c))
        if not ch: return ()
        choices.append(ch)
    result = {}
    for verts in itertools.product(*choices):
        initial = tuple((i, g) for i, (w, c) in enumerate(verts, 1) for g in w) \
            + tuple((A, g, k) for k, g in enumerate(namesA)) + tuple((B, g, k) for k, g in enumerate(namesB))
        coef = 1
        for w, c in verts: coef = coef * c
        p = {initial: coef}
        for (i, j) in reversed(phys):
            p = T.wiedge(p, i, j, slot.get((i, j)))
            if not p: break
        for w, c in p.items():
            result[w] = (result.get(w, 0) + c)
    return tuple((w, c) for w, c in result.items() if c)


def dtot_pos(word, pos, k):
    """total derivative d_1^{k1} d_2^{k2} acting on the letters at the positions pos"""
    cur = {word: 1}
    for axis in (0, 1):
        for _ in range(k[axis]):
            new = {}
            for w, c in cur.items():
                for i in pos:
                    u = list(w[i]); u[axis + 1] += 1
                    w2 = w[:i] + (tuple(u),) + w[i + 1:]
                    new[w2] = new.get(w2, 0) + c
            cur = new
    return cur.items()


def star_n(n, Ap, Bp, only=None):
    """[hbar^n] of A * B (stored normalisation).  only: restrict to one topology index."""
    if n == 0: return mul(Ap, Bp)
    if not Ap or not Bp: return {}
    if n > sum(u[0][0] == 'b' for u in next(iter(Ap))) + sum(u[0][0] == 'b' for u in next(iter(Bp))): return {}
    out = {}
    tops = topologies(n)
    kn = KAPPA[0] ** n
    for wa, ca in Ap.items():
        namesA = tuple(u[0] for u in wa)
        for wb, cb in Bp.items():
            namesB = tuple(u[0] for u in wb)
            for ti, (edges, roots, aut, name) in enumerate(tops):
                if only is not None and ti != only: continue
                A, B = n + 1, n + 2
                phys = edges[:-1]
                ext = [e for e in phys if e[1] > n]
                nA = sum(e[1] == A for e in ext); nB = len(ext) - nA
                if nA > len(wa) or nB > len(wb): continue
                pref = ca * cb * kn * inv(aut)
                for aA in itertools.permutations(range(len(wa)), nA):
                    for aB in itertools.permutations(range(len(wb)), nB):
                        res = wick(n, ti, namesA, namesB, aA, aB)
                        if not res: continue
                        itA, itB = iter(aA), iter(aB)
                        jets = tuple((lambda u: (u[1], u[2]))(wa[next(itA)] if e[1] == A else wb[next(itB)]) for e in ext)
                        mom = moments(n, ti, jets)
                        if not mom: continue
                        for w, c in res:
                            for jj, k in mom:
                                letters, apos, ok = [], [], True
                                for tag in w:
                                    if len(tag) == 2:
                                        ri = roots.index(tag[0])
                                        letters.append((tag[1], jj[ri][0], jj[ri][1]))
                                    elif tag[0] == A:
                                        apos.append(len(letters)); letters.append(wa[tag[2]])
                                    else:
                                        letters.append(wb[tag[2]])
                                # a root whose surviving letter is absent cannot carry derivatives (does not occur for cubic vertices)
                                if not apos and any(jj[-1]): continue
                                for word, cc in dtot_pos(tuple(letters), tuple(apos), jj[-1]):
                                    kk, s = T.canon(word)
                                    if s: out[kk] = (out.get(kk, 0) + s * cc * c * k * pref)
    return fin(out)


def Qn(n, p): return T.q0(p) if n == 0 else T.qn(n, p)


def par_of(Ap): return sum(T.parity(u) for u in next(iter(Ap))) % 2 if Ap else 0


def defect(n, Ap, Bp):
    """order hbar^n of  Q(A * B) - (QA) * B - (-1)^|A| A * (QB)  on the full letter algebra"""
    sg = -1 if par_of(Ap) else 1
    tot = {}
    for k in range(n + 1):
        m = n - k
        t = sub(sub(Qn(k, star_n(m, Ap, Bp)), star_n(m, Qn(k, Ap), Bp)), sc(star_n(m, Ap, Qn(k, Bp)), sg))
        tot = add(tot, t)
    return tot
