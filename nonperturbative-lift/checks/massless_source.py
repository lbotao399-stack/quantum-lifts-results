"""Complete zero-jet source of the Cartan obstruction in four-flavor magnetic SU(2)."""
from collections import Counter
from itertools import product

import sympy as sp


flavors = range(4)
colors = range(2)
heavy = 3


def letter(species):
    return lambda *indices: (species, *indices)


v, t, meson = map(letter, ("v", "t", "M"))
beta_v, beta_t, beta_m = map(letter, ("pv", "pt", "bM"))


def add(polynomial, coefficient, *fields):
    # Every output contains exactly one odd letter.
    polynomial[tuple(sorted(fields))] += coefficient


def charged_columns(columns, triangle):
    for species in ("v", "t"):
        for i, k, ell in product(flavors, repeat=3):
            # bM_(i,4) v_ell beta_v^k, or its right-flavor counterpart.
            if sorted((i, ell)) != sorted((heavy, k)):
                continue
            key = (species, i, k, ell)
            column = columns[key] = Counter()
            triangle[key] = int(i == heavy and ell == k)
            row, col = (i, heavy) if species == "v" else (heavy, i)
            for alpha, sigma in product(colors, repeat=2):
                scalar = v(ell, alpha) if species == "v" else t(ell, alpha)
                odd = beta_v(k, alpha) if species == "v" else beta_t(k, alpha)
                add(column, -1, v(row, sigma), t(col, sigma), scalar, odd)
            for alpha, j in product(colors, flavors):
                if species == "v":
                    add(column, 1, beta_m(row, col), v(ell, alpha), meson(k, j), t(j, alpha))
                else:
                    add(column, 1, beta_m(row, col), t(ell, alpha), meson(j, k), v(j, alpha))


def meson_columns(columns):
    for i, j, k, ell, r, s in product(flavors, repeat=6):
        if (i, j) >= (k, ell):  # exterior order of the two odd meson letters
            continue
        if sorted((i, k)) != sorted((heavy, r)) or sorted((j, ell)) != sorted((heavy, s)):
            continue
        column = columns[("M", i, j, k, ell, r, s)] = Counter()
        for alpha in colors:
            add(column, -1, meson(r, s), v(i, alpha), t(j, alpha), beta_m(k, ell))
            add(column, 1, meson(r, s), v(k, alpha), t(ell, alpha), beta_m(i, j))


def gauge_column(columns):
    # Twice the traceless projection in dimension two. Nonzero common factors
    # on this column and the matter columns are units for this kernel test.
    column = columns[("gauge",)] = Counter()
    for i, alpha, sigma in product(flavors, colors, colors):
        add(column, 2, t(heavy, alpha), v(i, alpha), beta_v(i, sigma), v(heavy, sigma))
        add(column, -1, t(heavy, alpha), v(heavy, alpha), v(i, sigma), beta_v(i, sigma))
        add(column, -2, t(heavy, alpha), beta_t(i, alpha), t(i, sigma), v(heavy, sigma))
        add(column, 1, t(heavy, alpha), v(heavy, alpha), beta_t(i, sigma), t(i, sigma))


def projected_evaluation(rows, i, species):
    # Take the coefficient of one odd variable, with one occupied heavy color
    # and, for a light flavor, use the other color. All unlisted fields vanish.
    alpha = 0 if i == heavy else 1
    scalar = v if species == "v" else t
    odd = beta_v if species == "v" else beta_t
    nonzero = {v(heavy, 0), t(heavy, 0), scalar(i, alpha), odd(i, alpha)}
    return sp.Matrix([[int(set(row) <= nonzero and odd(i, alpha) in row) for row in rows]])


def massive_evaluation():
    # mu is invertible, while m is also tested at zero.
    mu, mass = sp.Symbol("mu", nonzero=True), sp.Symbol("m")
    quarks = sp.Matrix(2, 4, sp.symbols("q:8"))
    antiquarks = sp.Matrix(4, 2, sp.symbols("t:8"))
    mesons = sp.Matrix(4, 4, sp.symbols("M:16"))
    potential = (quarks * mesons * antiquarks).trace() / mu + mass * mesons[3, 3]
    evaluation = dict.fromkeys([*quarks, *antiquarks, *mesons[:, 3], *mesons[3, :]], 0)
    evaluation.update({quarks[0, 3]: 1, antiquarks[3, 0]: -mu*mass})
    assert all(potential.diff(field).subs(evaluation) == 0
               for field in [*quarks, *antiquarks, *mesons])
    target = -3*mesons[:3, :3].det()/mu**3
    residue = sp.cancel(mass*target/(-mu*mass))
    assert sp.cancel(residue - 3*mesons[:3, :3].det()/mu**4) == 0
    assert residue.subs(dict(zip(mesons[:3, :3], sp.eye(3)))) == 3/mu**4
    print("All 32 tree F-equations vanish on the mass-polynomial family; determinant residue checked.")

    loop, coupling = sp.Symbol("hbar"), sp.Symbol("g", nonzero=True)
    triangle = 3*coupling**2/(64*sp.pi**2)
    # Columns e_m e_h Z, e_h C, e_m A. Rows e_h Z, e_m Z, C, A.
    total = sp.Matrix([[mass, mu*mass, 0], [-loop, 0, triangle*loop],
                       [0, loop, 0], [0, 0, mass]])
    residual = total * sp.Matrix([1, -1/mu, 1/triangle])
    assert residual == sp.Matrix([0, 0, -loop/mu, mass/triangle])
    assert sp.cancel(-mu*mass*residual[2] - triangle*loop*residual[3]) == 0
    print("The two-parameter Koszul variation cancels both intermediate rows with the source signs.")


def main():
    columns, triangle = {}, {}
    charged_columns(columns, triangle)
    meson_columns(columns)
    gauge_column(columns)
    rows = sorted({term for col in columns.values() for term, value in col.items() if value})
    tree = sp.Matrix([[col[row] for col in columns.values()] for row in rows])
    trace = sp.Matrix([[triangle.get(key, 0) for key in columns]])
    counts = Counter(key[0] for key in columns)
    assert counts == {"v": 7, "t": 7, "M": 24, "gauge": 1}
    assert tree.shape == (142, 39)
    assert tree.rank() == tree.col_join(trace).rank() == 38
    for i in flavors:
        evaluation = projected_evaluation(rows, i, "v") + projected_evaluation(rows, i, "t")
        diagonal = sp.Matrix([[int(key == ("v", heavy, i, i) or key == ("t", heavy, i, i))
                               for key in columns]])
        assert evaluation * tree == -diagonal
    print("39 invariant sources, 142 output monomials, tree rank 38, augmented rank 38.")
    print("All four flavor-resolved row identities hold on every column.")
    massive_evaluation()


if __name__ == "__main__":
    main()
