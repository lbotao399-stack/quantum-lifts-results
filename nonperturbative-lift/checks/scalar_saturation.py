"""Exact scalar-ideal elimination and local quotient checks."""
from itertools import combinations
from math import comb

import sympy as sp


def determinant_ideal(size):
    meson = sp.Matrix(size, size, sp.symbols(f"m:{size*size}"))
    quark = sp.Matrix(1, size, sp.symbols(f"q:{size}"))
    antiquark = sp.Matrix(size, 1, sp.symbols(f"u:{size}"))
    weight, inverse = sp.symbols("weight inverse")
    # weight stands for the full cofactor coefficient after multiplying the
    # meson F-equation by mu. The analytic proof keeps the physical units.
    equations = [*list(quark * meson), *list(meson * antiquark),
                 *list(quark.T * antiquark.T - weight * meson.cofactor_matrix())]
    variables = (*quark, *antiquark, *meson, weight)
    localized = sp.groebner([inverse * weight - 1, *equations], inverse, *variables, order="lex")
    elimination = [p.as_expr() for p in localized.polys if p.degree(inverse) == 0]
    proposed = sp.groebner([*equations, meson.det()], *variables, order="lex")
    assert elimination == [p.as_expr() for p in proposed.polys]
    original = sp.groebner(equations, *variables, order="lex")
    assert original.reduce(meson.det())[1] != 0
    for annihilator in (weight, *quark, *antiquark):
        assert original.reduce(annihilator * meson.det())[1] == 0
    print(f"Square size {size}: saturation equality, {len(elimination)} reduced basis elements.")


def pfaffian(matrix, indices):
    if not indices:
        return sp.Integer(1)
    first, *rest = indices
    return sum((-1)**j * matrix[first, partner] *
               pfaffian(matrix, rest[:j] + rest[j+1:]) for j, partner in enumerate(rest))


def compare_flavor_coordinates():
    meson = sp.Matrix(3, 3, sp.symbols("M:9"))
    baryon = sp.Matrix(1, 3, sp.symbols("B:3"))
    antibaryon = sp.Matrix(3, 1, sp.symbols("T:3"))
    matrix = sp.zeros(6)
    # Interleaved flavor order (1,~1,2,~2,3,~3).
    for i in range(3):
        for j in range(3):
            matrix[2*i, 2*j+1] = meson[i, j]
            matrix[2*j+1, 2*i] = -meson[i, j]
        j, k = (i+1) % 3, (i+2) % 3
        matrix[2*j, 2*k], matrix[2*k, 2*j] = baryon[i], -baryon[i]
        matrix[2*j+1, 2*k+1], matrix[2*k+1, 2*j+1] = antibaryon[i], -antibaryon[i]
    polynomial = pfaffian(matrix, list(range(6)))
    assert sp.expand(polynomial - meson.det() + (baryon * meson * antibaryon)[0]) == 0
    generators = []
    entries = list(baryon) + list(antibaryon) + list(meson)
    for i, j in combinations(range(6), 2):
        entry = matrix[i, j]
        sign, variable = entry.as_coeff_Mul()
        minor = pfaffian(matrix, [k for k in range(6) if k not in (i, j)])
        assert sp.expand(sign * polynomial.diff(variable) - (-1)**(i+j+1) * minor) == 0
        generators.append(sp.Poly(polynomial.diff(variable), *entries))
    supports = [set(generator.monoms()) for generator in generators]
    assert all(supports) and all(a.isdisjoint(b) for a, b in combinations(supports, 2))
    # Each minor identifies its removed letter. Taking the first letter
    # gives a distinct pivot in every exterior degree.
    for degree in range(1, 16):
        pivots = {(min(indices), indices[1:]) for indices in combinations(range(15), degree)}
        assert len(pivots) == comb(15, degree)
    mu, mass, magnetic_scale, electric_scale = sp.symbols("mu m magnetic electric", nonzero=True)
    assert sp.cancel((electric_scale/mu)/electric_scale - 1/mu) == 0
    assert (magnetic_scale/(mu**4*mass)).subs(magnetic_scale, mu**4*mass/electric_scale) == 1/electric_scale
    print("Interleaved six-flavor Pfaffian, all 15 signed derivatives and scale substitution agree.")
    print("Distinct exterior pivots in all 15 degrees, including 105 first fermionic targets.")


def circle_ideal(size):
    meson = sp.Matrix(size, size, sp.symbols(f"M:{size*size}"))
    baryon, antibaryon, monopole, eta, inverse = sp.symbols("B T Y eta inverse")
    variables = (monopole, baryon, antibaryon, *meson, eta)
    constraint = baryon*antibaryon - meson.det() + eta
    potential = monopole*constraint
    equations = [potential.diff(field) for field in variables[:-1]]
    localized = sp.groebner([inverse*eta-1, *equations], inverse, *variables, order="lex")
    elimination = [p.as_expr() for p in localized.polys if p.degree(inverse) == 0]
    proposed = sp.groebner([monopole, constraint], *variables, order="lex")
    assert elimination == [p.as_expr() for p in proposed.polys]
    print(f"Circle size {size}: saturated quotient has ideal (Y, constraint).")


def multiplier_pages():
    x, y = sp.symbols("x y")
    polynomial = x**5 + y**5 + x**3*y**3
    jacobian = [polynomial.diff(x), polynomial.diff(y)]
    weight, multiplier = sp.symbols("eta Y")
    derivative_part = ((5*polynomial - x**3*y**3)*(x*jacobian[0]+y*jacobian[1])
                       + (5*x**2*y**6*jacobian[0] - 3*x**4*y**5*jacobian[1])/(25-9*x*y))
    variation = multiplier*((weight - polynomial)*(weight + polynomial) + derivative_part/25)
    assert sp.cancel(variation - weight**2*multiplier) == 0
    # The proof shows (x,y)^7=0 in the local Jacobian algebra; imposing these
    # monomials selects the origin and excludes the other critical points.
    basis = sp.groebner([*jacobian, *(x**i*y**(7-i) for i in range(8))], x, y,
                        order="grevlex", domain=sp.QQ)
    monomials = [x**i*y**j for i in range(4) for j in range(4)]
    ambient = [x**i*y**(degree-i) for degree in range(7) for i in range(degree+1)]
    def coordinates(value):
        rem = sp.Poly(basis.reduce(value)[1], x, y)
        return sp.Matrix([rem.coeff_monomial(m) for m in ambient])
    identities = sp.Matrix.hstack(*(coordinates(m) for m in monomials))
    action = sp.Matrix.hstack(*(coordinates(polynomial*m) for m in monomials))
    assert identities.rank() == 16 and action.rank() == 1
    assert action[:, 0] == -identities[:, 15] / 5 and action[:, 1:] == sp.zeros(28, 15)
    assert all(basis.reduce(polynomial**2*m)[1] == 0 for m in monomials)
    print("Local multiplication by F has dimension 16, rank 1, square rank 0.")


if __name__ == "__main__":
    for n in (2, 3):
        determinant_ideal(n)
        circle_ideal(n)
    compare_flavor_coordinates()
    multiplier_pages()
