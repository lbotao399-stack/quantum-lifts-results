"""Exact tests of the effective determinant and finite-circle relations in chapters 16 and 19."""
import json
import sympy as s


def check():
    mu, mass, scale, eta = s.symbols('mu mass scale eta', nonzero=True)
    M = s.Matrix(3, 3, s.symbols('M0:9'))
    q, t = s.symbols('q0:3'), s.symbols('t0:3')
    cubic = (s.Matrix(q).T * M * s.Matrix(t))[0] / mu
    def euler(f):
        return sum(x*s.diff(f,x) for x in M)-sum(x*s.diff(f,x) for x in q)
    W = cubic - scale*M.det()/(mu**4*mass)
    assert s.expand(euler(cubic)) == 0
    assert s.expand(euler(W) + 3*scale*M.det()/(mu**4*mass)) == 0
    flat = dict.fromkeys(q+t, 0)
    assert all(s.diff(cubic, x).subs(flat)==0 for x in [*M,*q,*t])
    assert M.det().subs(dict(zip(M,s.eye(3))))==1
    Y, B, Bt = s.symbols('Y B Bt')
    for n in (2, 3):
        meson = s.Matrix(n, n, s.symbols(f'x:{n*n}'))
        V = Y*(B*Bt-meson.det()+eta)
        residual = Y*s.diff(V,Y) - B*s.diff(V,B)-sum(x*s.diff(V,x) for x in meson)/n - eta*Y
        assert s.expand(residual)==0
        origin = {**dict.fromkeys(meson,0), B:0, Bt:0, eta:0}
        assert all(s.diff(V,x).subs(origin)==0 for x in [*meson, B, Bt, Y])
    return {'determinant_size':3, 'Euler_residual':0,'perturbatively_flat_evaluation_det':1,
            'circle_sizes':[2,3], 'circle_Euler_residuals':[0,0], 'circle_Y_free_at_eta_zero':True,
            'scope':'effective polynomial identities in these matrix sizes; the physical superpotential is input; no microscopic instanton insertion'}


if __name__ == '__main__':
    if not __debug__: raise SystemExit('Checks require enabled assertions (remove -O).')
    print(json.dumps(check(),indent=2))
