"""Solve finite polynomial Bianchi equations retaining every homogeneous mark.

The unknown positive-order coefficients are all six native adjoint two-forms.
Column elimination is exact over Q. No colour or cohomology coordinates are
discarded from the original equations.
"""

from itertools import combinations

from flint import fmpq, fmpq_mat


def simultaneous_marks(colour, quartet, curve, leading_forms, degree):
    if not __debug__:
        raise RuntimeError("enable assertions for certificate verification")
    if type(degree) is not int or degree < 1:
        raise ValueError("mark degree must be a positive integer")
    from nilpotent_cubic import raw_bracket
    from nilpotent_quadratic import add_term

    xs = [{((i,), c): fmpq(1) for i, c in enumerate(quartet)}, *curve[1:]]
    curvature = {}
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            for key, value in raw_bracket(colour, x, y).items():
                add_term(curvature, (i + j, key), value)
    assert not curvature
    native = [(axes, c) for axes in combinations(range(4), 2) for c in range(52)]
    unknowns = [(n, k) for n in range(1, degree + 1) for k in native]
    columns = []
    for order, key in unknowns:
        col = {}
        for j, x in enumerate(xs):
            for target, value in raw_bracket(colour, x, {key: fmpq(1)}).items():
                add_term(col, (order + j, target), value)
        columns.append(col)
    for _, form in leading_forms:
        col = {}
        for j, x in enumerate(xs):
            for target, value in raw_bracket(colour, x, form).items():
                add_term(col, (j, target), value)
        columns.append(col)
    rows = sorted({row for col in columns for row in col})
    index = {key: i for i, key in enumerate(rows)}
    matrix = fmpq_mat(len(rows), len(columns))
    for j, col in enumerate(columns):
        for key, value in col.items():
            matrix[index[key], j] = value
    reduced, rank = matrix.rref()
    pivots = [next(j for j in range(reduced.ncols()) if reduced[i, j]) for i in range(rank)]
    free_leading = [j for j in range(len(unknowns), len(columns)) if j not in pivots]
    leading = fmpq_mat(len(leading_forms), len(free_leading))
    certificates = []
    for pos, free in enumerate(free_leading):
        solution = {free: fmpq(1)}
        for i, pivot in enumerate(pivots):
            if reduced[i, free]:
                solution[pivot] = -reduced[i, free]
        series = [{} for _ in range(degree + 1)]
        for j, value in solution.items():
            if j < len(unknowns):
                n, key = unknowns[j]
                add_term(series[n], key, value)
            else:
                m = j - len(unknowns)
                leading[m, pos] = value
                for key, coefficient in leading_forms[m][1].items():
                    add_term(series[0], key, value * coefficient)
        # Verify every native polynomial coefficient, not an evaluation subset.
        full = {}
        for i, x in enumerate(xs):
            for j, y in enumerate(series):
                for key, value in raw_bracket(colour, x, y).items():
                    add_term(full, (i + j, key), value)
        assert not full
        certificates.append(series)
    return leading, certificates, {
        "native_equation_rows": len(rows),
        "unknown_higher_mark_coefficients": len(unknowns),
        "leading_mark_coefficients": len(leading_forms),
        "rank": rank,
        "allowed_leading_marks": len(free_leading),
        "all_native_coefficients_verified": True,
        "all_scalar_commutator_coefficients_verified": True,
    }
