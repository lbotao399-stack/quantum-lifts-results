"""Exact-Q adapter for the preserved Claude product engine.

No original source is modified. SUN chooses SU(N), as in sunlib.py.
Only the coefficient ring and helpers are replaced; graph formulas are unchanged.
"""
from pathlib import Path
from fractions import Fraction as F
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'product'))
import sunlib as S

T = S.ITheory(S.fields, S.partners, S.I, 0, 1,
              ghosts=['c' + A for A in S.labels])
S.T = T
S._img.clear()
S.fin = T.fin
S.sc = T.scale
S.sub = lambda p, q: T.add(p, T.scale(q, F(-1)))
S.inv = lambda x: F(1, 1) / F(x)

# Rebuild the matrix entries in Q, rather than reconstructing their residues.
S.ENT = {True: {}, False: {}}
for i in range(S.N):
    for j in range(S.N):
        S.ENT[True][i, j] = [(S.labels[a], F(S.mats[a][i][j]))
                              for a in range(S.nl) if S.mats[a][i][j]]
        S.ENT[False][i, j] = []
        for a in range(S.nl):
            value = sum(S.kinv[a][b] * S.mats[b][i][j]
                        for b in range(S.nl))
            if value:
                S.ENT[False][i, j].append((S.labels[a], value))


def cyc_poly(letters):
    m = len(letters)
    out = {}
    entries = [S.ENT[S.UP[letter[0]]] for letter in letters]
    for indices in S.itertools.product(range(S.N), repeat=m):
        choices = [entries[j][indices[j], indices[(j + 1) % m]]
                   for j in range(m)]
        for colors in S.itertools.product(*choices):
            coeff = F(1)
            for _, value in colors:
                coeff *= value
            word, sign = T.canon(tuple((letters[j][0] + colors[j][0],
                                      letters[j][1], letters[j][2])
                                     for j in range(m)))
            if sign:
                out[word] = out.get(word, 0) + sign * coeff
    return T.fin(out)


S.cyc_poly = cyc_poly
import ops as O


def comps(polynomial):
    out = {}
    for word, value in polynomial.items():
        probes = [letter for letter in word if letter[0][0] == 'p']
        assert len(probes) == 1
        rest = tuple(letter for letter in word if letter[0][0] != 'p')
        color = probes[0][0][1:]
        out.setdefault(color, {})
        out[color][rest] = out[color].get(rest, 0) + value
    return {color: T.fin(poly) for color, poly in out.items()}


O.comps = comps
import mu1 as M1
M1.KAPPA = {'bc': 2, 'cb': -2, 'bb': -2}
import mu2 as M2

Q0 = S.Q0
Q1 = S.Q1
Q2 = S.Q2
Q3 = lambda p: T.qn(3, p)
mu1 = M1.mu1
mu2 = M2.mu2
mul = T.mul
add = T.add
sc = T.scale
sub = S.sub
v = T.v
cartan = O.cartan
Jp = S.Jp
deriv = O.deriv


if __name__ == '__main__':
    C0 = cyc_poly([('c', 1, 0), ('c', 0, 1)])
    for color in S.labels:
        b = v('b' + color)
        assert not M1.defect(C0, b)
        defect = sub(mu1(mul(C0, C0), b), sc(mul(C0, mu1(C0, b)), 2))
        print(color, 'multiplication-composition defect terms', len(defect))
    print('Exact-Q product adapter: SU(%s), C0/b Leibniz checks PASS' % S.N)
