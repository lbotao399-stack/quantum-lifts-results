"""Exact regular-scalar homotopy tests in the four-colour F4 polarization."""

from itertools import combinations
import json
from random import Random

from flint import fmpq, fmpq_mat

from colour import Colour
from commuting_source import levi_strata, simple_coordinates


def pivots(matrix):
    echelon, rank = matrix.rref()
    return [next(j for j in range(matrix.ncols()) if echelon[i, j]) for i in range(rank)]


def adjoint(colour, element):
    action = fmpq_mat(52, 52)
    for i, coefficient in enumerate(element.entries()):
        if coefficient:
            for j, bracket in enumerate(colour.br[i]):
                for k, value in bracket.items():
                    action[k, j] += coefficient * value
    return action


def bracket(colour, left, right):
    return adjoint(colour, left) * right


def minor_inverse(action):
    rows = pivots(action.transpose())
    columns = pivots(action)
    assert len(rows) == len(columns) == 48
    minor = fmpq_mat([[action[i, j] for j in columns] for i in rows])
    inverse = minor.inv()
    lifted = fmpq_mat(52, 52)
    for i, column in enumerate(columns):
        for j, row in enumerate(rows):
            lifted[column, row] = inverse[i, j]
    return lifted, rows, columns, minor.det()


def discriminant_inverse(action, coefficients, identity):
    # det(t-ad X) = t^4 p(t).  Evaluate (p(0)-p(t))/(p(0)t) by Horner.
    assert not any(coefficients[:4]) and coefficients[4] and not any(coefficients[5::2])
    polynomial = fmpq_mat(52, 52)
    for coefficient in reversed(coefficients[5:]):
        polynomial = action * polynomial - coefficient * identity
    return polynomial / coefficients[4]


def representative(colour, nodes, nilpotent_nodes):
    simple = colour.data["raising"]
    cartan = fmpq_mat([colour.weights[i] for i in simple])
    element = fmpq_mat(52, 1)
    for i, value in enumerate(cartan.solve(fmpq_mat([[int(i not in nodes)] for i in range(4)])).entries()):
        element[i, 0] = value
    for i in nilpotent_nodes:
        element[simple[i], 0] = 1
    return element


def check_homotopy(colour, element, action, inverse, identity, seed, levi=None):
    assert inverse * action * inverse == inverse
    scalar_projection = identity - inverse * action
    mark_projection = identity - action * inverse
    independent = pivots(scalar_projection)
    if levi is None:
        assert action * inverse * action == action and len(independent) == 4
        central = [fmpq_mat([[scalar_projection[i, j]] for i in range(52)]) for j in independent]
        assert all(not bracket(colour, x, y) for x, y in combinations(central, 2))
    else:
        assert scalar_projection == mark_projection == levi
    random = Random(seed)

    def vector():
        return fmpq_mat([[random.randrange(-3, 4)] for _ in range(52)])

    scalars = [element, vector(), vector(), vector()]
    pairs = list(combinations(range(4), 2))
    marks = {pair: vector() for pair in pairs}
    splits = {(i, j): bracket(colour, scalars[i], scalars[j]) for i, j in pairs}
    centres = [scalar_projection * x for x in scalars]
    normals = [x - centre for x, centre in zip(scalars, centres)]
    assert centres[0] == element

    def path_on_marks(values, time):
        result = {}
        for i, j in pairs:
            if i == 0:
                result[i, j] = mark_projection * values[i, j] + time * (action * inverse * values[i, j])
                continue
            a = inverse * values[0, i]
            b = inverse * values[0, j]
            linear = bracket(colour, a, centres[j]) + bracket(colour, centres[i], b)
            quadratic = bracket(colour, a, normals[j]) + bracket(colour, normals[i], b)
            result[i, j] = values[i, j] - (1 - time) * linear - ((1 - time**2) / 2) * quadratic
        return result

    chain_tests = derivative_tests = 0
    for time in (fmpq(0), fmpq(2, 3), fmpq(1)):
        at = [x + time * y for x, y in zip(centres, normals)]
        on_splits = path_on_marks(splits, time)
        before = path_on_marks(marks, time - 1)
        after = path_on_marks(marks, time + 1)
        for i, j in pairs:
            assert on_splits[i, j] == bracket(colour, at[i], at[j])
            first = fmpq_mat(52, 1) if i == 0 else inverse * marks[0, i]
            derivative = bracket(colour, first, at[j]) + bracket(colour, at[i], inverse * marks[0, j])
            assert (after[i, j] - before[i, j]) / 2 == derivative
            chain_tests += 1
            derivative_tests += 1
    assert path_on_marks(marks, fmpq(1)) == marks
    assert any(splits.values())
    if levi is None:
        assert not any(path_on_marks(splits, fmpq(0)).values())
        return dict(chain_tests=chain_tests, derivative_tests=derivative_tests,
                    abelian_centralizer_pairs=6, different_scalar_and_mark_projections=scalar_projection != mark_projection)

    at_zero = path_on_marks(marks, fmpq(0))
    assert any(path_on_marks(splits, fmpq(0)).values())
    transgressions = [fmpq_mat(52, 1), *(inverse * marks[0, i] for i in (1, 2, 3))]
    thirds = list(combinations(range(4), 3))
    integrals = {triple: fmpq_mat(52, 1) for triple in thirds}
    # The integrand Q0 H_t Q0 Z is quadratic in t. Simpson's rule is exact.
    for time, weight in ((fmpq(0), fmpq(1, 6)), (fmpq(1, 2), fmpq(2, 3)), (fmpq(1), fmpq(1, 6))):
        on_marks = path_on_marks(marks, time)
        on_splits = path_on_marks(splits, time)
        for triple in thirds:
            i, j, k = triple
            for a, pair, sign in ((i, (j, k), 1), (j, (i, k), -1), (k, (i, j), 1)):
                integrals[triple] += weight * sign * (
                    bracket(colour, normals[a], on_marks[pair]) -
                    bracket(colour, transgressions[a], on_splits[pair]))
    for i, j, k in thirds:
        original, end, restricted = (fmpq_mat(52, 1) for _ in range(3))
        for a, pair, sign in ((i, (j, k), 1), (j, (i, k), -1), (k, (i, j), 1)):
            original += sign * bracket(colour, scalars[a], marks[pair])
            end += sign * bracket(colour, centres[a], at_zero[pair])
            restricted += sign * bracket(colour, centres[a], levi * at_zero[pair])
        assert original - integrals[i, j, k] == end and levi * end == restricted
    return dict(levi_dimension=len(independent), chain_tests=chain_tests, derivative_tests=derivative_tests,
                third_jet_integrals=len(thirds))


def calculate():
    colour = Colour()
    strata = levi_strata(colour)
    identity = fmpq_mat(52, 52)
    for i in range(52):
        identity[i, i] = 1
    records = []
    slices = []
    roots = simple_coordinates(colour)
    for number, stratum in enumerate(strata["levi_subsystems"]):
        nodes = stratum["simple_nodes"]
        element = representative(colour, nodes, nodes)
        action = adjoint(colour, element)
        inverse, rows, columns, determinant = minor_inverse(action)
        coefficients = action.charpoly().coeffs()
        record = dict(simple_nodes=nodes, rows=rows, columns=columns, minor=str(determinant),
                      discriminant=str(coefficients[4] if len(coefficients) > 4 else 0),
                      minor_homotopy=check_homotopy(colour, element, action, inverse, identity, 19436 + number))
        if coefficients[4]:
            assert not nodes
            canonical = discriminant_inverse(action, coefficients, identity)
            record["invariant_homotopy"] = check_homotopy(colour, element, action, canonical, identity, 30652)
        records.append(record)
        if 0 < len(nodes) < 4:
            complement = [i for i, root in enumerate(roots) if any(root[j] for j in range(4) if j not in nodes)]
            nilpotent = nodes[:1] if len(nodes) > 1 else ()
            anchor = representative(colour, nodes, nilpotent)
            move = adjoint(colour, anchor)
            block = fmpq_mat([[move[i, j] for j in complement] for i in complement])
            back = block.inv()
            lift = fmpq_mat(52, 52)
            projection = fmpq_mat(identity)
            for i, row in enumerate(complement):
                projection[row, row] = 0
                for j, column in enumerate(complement):
                    lift[row, column] = back[i, j]
            tests = check_homotopy(colour, anchor, move, lift, identity, 71204 + number, projection)
            slices.append(dict(nodes=nodes, nilpotent_nodes=nilpotent, adjoint_rank=move.rank(),
                               complementary_determinant=str(block.det()), tests=tests))
    assert len(records) == 16
    assert records[-1]["simple_nodes"] == (0, 1, 2, 3) and records[-1]["discriminant"] == "0"
    assert records[-1]["minor_homotopy"]["different_scalar_and_mark_projections"]
    assert all(min(t["required_nilpotent_slots"]) > t["max_height"]
               for t in strata["levi_subsystems"] if len(t["simple_nodes"]) >= 2)
    return dict(ordinary_diagram=[4, 4, 4, 2], scalar_diagrams=strata["scalar_coefficients"],
                adjoint_rank=48, scalar_normal_coordinates=144, marked_generators=312,
                redundant_marked_generators=168, homogeneous_minor_powers=[10, 10, 10, 12],
                records=records, levi_slices=slices)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(), indent=2))
