"""Rational certification of the SU(2) hidden-product calculation.

Uses modular witnesses to select a small exact linear system. Every final
polynomial identity is recomputed and checked over fractions.Fraction.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, gcd
import pickle, json, time, resource, hashlib
import numpy as np
from flint import fmpq, fmpq_mat, ctx
ctx.threads = 1
import exact_product as E

ROOT = Path(__file__).resolve().parent
P = 2147483647
BOUND = isqrt(P // 2)
START = time.time()
REPORT = []
CHECK_PATH = ROOT / 'rational_checks.json'


def say(message):
    print(message, flush=True)


def rat(x):
    x = int(x) % P
    if not x:
        return F(0)
    r0, r1, s0, s1 = P, x, 0, 1
    while abs(r1) > BOUND:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q*r1, s1, s0 - q*s1
    if s1 < 0:
        r1, s1 = -r1, -s1
    if not 0 < s1 <= BOUND or gcd(r1, s1) != 1 or (r1-x*s1) % P:
        raise ValueError(('rational reconstruction failed', x))
    return F(r1, s1)


def poly(p):
    return {w: rat(c) for w, c in p.items()}


def load(name):
    with (ROOT / name).open('rb') as f:
        return pickle.load(f)


def check(name, residual):
    result = {'name': name, 'residual_terms': len(residual),
              'passed': not residual, 'seconds': time.time()-START}
    REPORT.append(result)
    say('%s: %s (%d terms, %.1fs)' %
        (name, 'PASS' if not residual else 'FAIL', len(residual), time.time()-START))
    CHECK_PATH.write_text(json.dumps(REPORT, indent=2))
    if residual:
        with (ROOT / 'failed_residual.pkl').open('wb') as f:
            pickle.dump(residual, f)
        raise AssertionError(name)


def modq(x):
    x = F(x)
    return x.numerator % P * pow(x.denominator % P, P-2, P) % P


def independent_coordinate_rows(columns):
    """Select a square nonsingular minor modulo P; exact verification follows."""
    monomials = sorted(set().union(*(p.keys() for p in columns)))
    basis = []
    selected = []
    for monomial in monomials:
        row = np.array([modq(p.get(monomial, 0)) for p in columns], dtype=np.int64)
        for pivot, old in basis:
            if row[pivot]:
                row = (row - row[pivot] * old % P) % P
        nz = np.flatnonzero(row)
        if len(nz):
            pivot = int(nz[0])
            row = row * pow(int(row[pivot]), P-2, P) % P
            basis.append((pivot, row))
            selected.append(monomial)
            if len(selected) == len(columns):
                return selected
    raise ValueError(('dependent selected columns', len(selected), len(columns)))


def rational_solve(columns, target):
    coordinates = independent_coordinate_rows(columns)
    fq = lambda x: fmpq(F(x).numerator, F(x).denominator)
    matrix = fmpq_mat([[fq(p.get(w, 0)) for p in columns] for w in coordinates])
    rhs = fmpq_mat([[fq(target.get(w, 0))] for w in coordinates])
    solution = matrix.solve(rhs)
    coefficients = [F(int(solution[i, 0].numerator), int(solution[i, 0].denominator))
                    for i in range(len(columns))]
    check('exact linear solve on every monomial',
          E.sub(target, E.add(*(E.sc(p, c) for p, c in zip(columns, coefficients)))))
    return coefficients, coordinates


def main():
    A = load('su2c4/c4A.pkl')
    B = load('su2c4/c4B.pkl')
    W = load('product/c2sq_witnesses.pkl')
    C0 = poly(A['C0'])
    a4, a3, a2_C4 = (poly(A[key]) for key in ('a4', 'a3', 'a2'))
    a2, a1, a0 = (poly(W[key]) for key in ('a2', 'a1', 'a0'))
    check('C0 agrees with trace formula', E.sub(C0, E.cyc_poly([('c',1,0),('c',0,1)])))
    check('C2: Q0 a2', E.Q0(a2))
    check('C2: Q0 a1 + Q1 a2', E.add(E.Q0(a1), E.Q1(a2)))
    check('C2: Q0 a0 + Q1 a1 + Q2 a2', E.add(E.Q0(a0), E.Q1(a1), E.Q2(a2)))
    check('C4: Q0 a4', E.Q0(a4))
    check('C4: Q0 a3 + Q1 a4', E.add(E.Q0(a3), E.Q1(a4)))
    check('C4: Q0 a2 + Q1 a3 + Q2 a4', E.add(E.Q0(a2_C4), E.Q1(a3), E.Q2(a4)))
    check('source normalization a2^2 = -C0 a4', E.add(E.mul(a2,a2), E.mul(C0,a4)))

    P4 = E.mul(a2, a2)
    P3 = E.add(E.sc(E.mul(a1,a2),2), E.mu1(a2,a2))
    P2 = E.add(E.sc(E.mul(a0,a2),2), E.mul(a1,a1), E.mu1(a2,a1),
               E.mu1(a1,a2), E.mu2(a2,a2))
    check('C2 product closure Nb=4', E.Q0(P4))
    check('C2 product closure Nb=3', E.add(E.Q1(P4), E.Q0(P3)))
    check('C2 product closure Nb=2', E.add(E.Q2(P4), E.Q1(P3), E.Q0(P2)))
    La4 = E.mul(C0, a4)
    La3 = E.add(E.mul(C0,a3), E.mu1(C0,a4))
    La2 = E.add(E.mul(C0,a2_C4), E.mu1(C0,a3), E.mu2(C0,a4))
    for name, actual in [('La4',La4),('La3',La3),('La2',La2)]:
        check('direct graph product equals stored '+name, E.sub(actual,poly(A[name])))

    k4p, k3p = poly(W['k4']), poly(W['k3'])
    k4, k3 = poly(B['k4']), poly(B['k3'])
    check('C2 product primitive Nb=4', E.sub(E.Q0(k4p),P4))
    check('C2 product primitive Nb=3', E.sub(E.Q0(k3p),E.sub(P3,E.Q1(k4p))))
    check('C0C4 product primitive Nb=4', E.sub(E.Q0(k4),La4))
    check('C0C4 product primitive Nb=3', E.sub(E.Q0(k3),E.sub(La3,E.Q1(k4))))
    r2p = E.sub(E.sub(P2,E.Q2(k4p)),E.Q1(k3p))
    r2 = E.sub(E.sub(La2,E.Q2(k4)),E.Q1(k3))
    check('C2 product reduced representative Q0 closed', E.Q0(r2p))
    check('C0C4 reduced representative Q0 closed', E.Q0(r2))
    say('Reduced representatives recomputed over Q: %d, %d terms' % (len(r2p),len(r2)))

    n0 = len(W['boundary_Q0_sources'])
    indices0 = [i for i,c in enumerate(W['boundary_coeffs'][:n0]) if c]
    indices1 = [i for i,c in enumerate(W['boundary_coeffs'][n0:-1]) if c]
    sources0 = [poly(W['boundary_Q0_sources'][i]) for i in indices0]
    sources1 = [poly(W['boundary_Q1_sources'][i]) for i in indices1]
    for i,p in enumerate(sources1):
        check('Q1-boundary source %d Q0 closed'%i,E.Q0(p))
    columns = [E.Q0(p) for p in sources0]+[E.Q1(p) for p in sources1]
    say('Solving rational boundary system with %d columns'%len(columns))
    coeffs,coordinates = rational_solve(columns,E.sub(r2p,E.sc(r2,3)))
    primitive0 = E.add(*(E.sc(p,c) for p,c in zip(sources0,coeffs[:len(sources0)])))
    primitive1 = E.add(*(E.sc(p,c) for p,c in zip(sources1,coeffs[len(sources0):])))
    check('hidden product boundary Q0 primitive + Q1 closed primitive',
          E.sub(E.sub(r2p,E.sc(r2,3)),E.add(E.Q0(primitive0),E.Q1(primitive1))))
    check('combined Q1 primitive is Q0 closed',E.Q0(primitive1))

    certificate = dict(C0=C0,C2=[a2,a1,a0],C4_partial=[a4,a3,a2_C4],
        C2_squared=[P4,P3,P2],C0_C4=[La4,La3,La2],
        C2_squared_primitives=[k4p,k3p],C0_C4_primitives=[k4,k3],
        r2p=r2p,r2=r2,boundary_Q0=primitive0,boundary_Q1=primitive1,
        coefficient=F(3),source_normalization_lambda=F(-1),
        selected_Q0_indices=indices0,selected_Q1_indices=indices1,
        rational_coefficients=coeffs,minor_coordinates=coordinates)
    path=ROOT/'c4_hidden_product_exact.pkl'
    with path.open('wb') as f:pickle.dump(certificate,f)
    summary=dict(checks=len(REPORT),passed=all(x['passed'] for x in REPORT),
        seconds=time.time()-START,maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        relation='[r2(C2 star C2)] = 3 [r2(C0 star C4_source)]',
        normalized_relation='C2^2 = -3 C0 C4',
        scope='Exact Q polynomial boundary identity at Nb=2. C4 full perturbative completion and non-exactness are separate checks.')
    (ROOT/'certificate_summary.json').write_text(json.dumps(summary,indent=2))
    say(json.dumps(summary))


if __name__=='__main__':
    main()
