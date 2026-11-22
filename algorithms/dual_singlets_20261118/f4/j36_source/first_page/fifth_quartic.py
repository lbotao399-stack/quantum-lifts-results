"""Contract the F4 long-root quartic in the lifted fermion measure of a unit instanton."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import comb, factorial, prod

import numpy as np

from colour import Colour
from fifth_instanton import calculate as instanton


def invariant_dimensions():
    # The weights of primitive three-forms are those of wedge^3 of the defining
    # C3 module, with one copy of the defining module removed.
    vector = [tuple(sign * (axis == j) for j in range(3))
              for axis in range(3) for sign in (-1, 1)]
    three_forms = Counter(tuple(map(sum, zip(*c))) for c in combinations(vector, 3))
    weights = list((three_forms - Counter(vector)).elements())
    assert len(set(weights)) == len(weights) == 14

    rho = (3, 2, 1)
    def invariant_count(choose, degree):
        character = Counter(tuple(map(sum, zip(*indices))) for indices in choose(weights, degree))
        answer = 0
        for perm in permutations(range(3)):
            sign = (-1)**sum(i > j for k, i in enumerate(perm) for j in perm[k+1:])
            for signs in product((-1, 1), repeat=3):
                delta = tuple(rho[i] - signs[i] * rho[perm[i]] for i in range(3))
                answer += sign * prod(signs) * character[delta]
        return answer

    exterior = [invariant_count(combinations, k) for k in range(1, 8)]
    symmetric = [invariant_count(combinations_with_replacement, k) for k in range(1, 7)]
    assert exterior == [0, 1, 0, 1, 0, 1, 0] and symmetric == [0, 0, 0, 1, 0, 0]
    return dict(exterior_degrees_1_to_7=exterior, symmetric_degrees_1_to_6=symmetric)


def native_quartic(c, endpoint):
    colours = endpoint["doublets"]
    low = endpoint["opposite_colour"]
    f = {low: endpoint["lowering_coefficient"]}
    def act(i, vector):
        result = Counter()
        for j, a in vector.items():
            for k, b in c.br[i][j].items():
                result[k] += a*b
        return result

    def pairing(left, right):
        return sum(a*b*c.metric.get((i, j), 0) for i, a in left.items() for j, b in right.items())

    # Keep twice the second polarization to make the construction integral.
    second = {}
    for i, j in product(colours, repeat=2):
        value = act(i, act(j, f))
        value.update(act(j, act(i, f)))
        second[i, j] = value

    @cache
    def diagonal(indices):
        value = f
        for _ in range(4):
            out = Counter()
            for i in indices:
                out.update(act(i, value))
            value = out
        return pairing(f, value)

    tensor = np.zeros((14,)*4, dtype=np.int64)
    components = []
    for indices in combinations_with_replacement(range(14), 4):
        i, j, k, l = word = tuple(colours[s] for s in indices)
        value = Fraction(pairing(second[i, j], second[k, l]) + pairing(second[i, k], second[j, l])
                         + pairing(second[i, l], second[j, k]), 12)
        # Independently polarize the diagonal four nested brackets by finite differences.
        alternate = sum((-1)**(4-m)*diagonal(term) for m in range(5)
                        for term in combinations(word, m))
        assert alternate == factorial(4)*value
        if value:
            assert (2*value).denominator == 1
            components.append([*word, str(value)])
            for perm in set(permutations(indices)):
                tensor[perm] = int(2*value)
    assert len(components) == 44
    return tensor, components


def inverse_pairs(endpoint):
    colours = endpoint["doublets"]
    partner = [None]*len(colours)
    numerators = [0]*len(colours)
    for left, right, b in endpoint["alternating_form"]:
        i, j = colours.index(left), colours.index(right)
        assert 4 % (2*b) == 0
        partner[i], partner[j] = j, i
        numerators[i], numerators[j] = 4//(2*b), -4//(2*b)
    assert all(i != j and partner[j] == i for i, j in enumerate(partner))
    return partner, numerators


def matchings(vertices):
    if not vertices:
        yield (), 1
        return
    a = vertices[0]
    for j, b in enumerate(vertices[1:]):
        tail = tuple(k for k in vertices[1:] if k != b)
        for pairs, sign in matchings(tail):
            yield ((a, b), *pairs), (-1)**j * sign


def graph_score(tensor, dual, matches):
    # Each edge carries exactly one independent colour and its symplectic partner.
    # Signs from moving vertices into flavour order are positive for four flavours.
    vertices = len(matches[0])*2
    sign = dual[1]
    arguments = []
    edge_numbers = [[0]*4 for _ in range(vertices)]
    crossed = [[False]*4 for _ in range(vertices)]
    for axis, matching in enumerate(matches):
        for j, (u, v) in enumerate(matching):
            edge_numbers[u][axis] = edge_numbers[v][axis] = axis*len(matching)+j
            crossed[v][axis] = True
    for i in range(vertices):
        local = tensor
        for axis in range(4):
            if crossed[i][axis]:
                # Contraction sums P^{ab} T(...,a,...) T(...,b,...).
                shape = [1]*4
                shape[axis] = 14
                local = local.take(dual[0], axis=axis)*np.asarray(sign).reshape(shape)
        arguments.extend([local, edge_numbers[i]])
    # With four vertices the absolute bound, before antisymmetrizing, is
    # 14^8 * 4^4 * 2^8; the contraction and its intermediates fit signed int64.
    assert vertices <= 4 and 14**8*4**4*2**8 < 2**63
    return int(np.einsum(*arguments, [], optimize="greedy"))


def saturated(tensor, dual, pfaffian, vertices):
    choices = tuple(matchings(tuple(range(vertices))))
    classes = []
    total = 0
    for order in combinations_with_replacement(range(len(choices)), 4):
        copies = factorial(4)//prod(factorial(n) for n in Counter(order).values())
        entries = [choices[i] for i in order]
        sign = prod(s for _, s in entries)
        raw = graph_score(tensor, dual, [m for m, _ in entries])
        score = Fraction(raw*sign, 2**vertices*4**(2*vertices))
        total += copies*score
        classes.append(dict(pairings=order, multiplicity=copies, signed_score=str(score)))
    assert sum(x["multiplicity"] for x in classes) == len(choices)**4
    power = 7 - vertices//2
    physical_colour = total * (factorial(power)*pfaffian)**4 / factorial(vertices)
    return dict(vertices=vertices, doublet_form_power_in_each_flavour=power,
                Wick_contractions=len(choices)**4, grouped=classes,
                contracted_quartic=str(total), top_coefficient=str(physical_colour))


def restricted_top(tensor, dual, endpoint, slots, vertices):
    # Extract the coefficient after keeping a specified set of colours in each
    # of the last three flavours. Direct exterior multiplication of the original
    # four-flavour vertex gives an independent test of all signs and factorials.
    width = len(slots)
    assert width == vertices
    def exterior_sign(left, right):
        swaps = 0
        while left:
            bit = left & -left
            swaps += (right & (bit-1)).bit_count()
            left ^= bit
        return (-1)**swaps

    def multiply(x, y):
        out = Counter()
        for a, s in x.items():
            for b, t in y.items():
                if not (a & b):
                    out[a|b] += exterior_sign(a, b)*s*t
        return {m: v for m, v in out.items() if v}

    currents = []
    for i in range(14):
        row = {}
        for a, b, d in product(range(width), repeat=3):
            value = int(tensor[i, slots[a], slots[b], slots[d]])
            if value:
                row[(1<<a)|(1<<(width+b))|(1<<(2*width+d))] = value
        currents.append(row)
    top = (1 << (3*width))-1
    inverse = Counter()
    pairs = [(i, j) for i, j in enumerate(dual[0]) if i < j]
    for i, j in pairs:
        for mask, v in multiply(currents[i], currents[j]).items():
            inverse[mask] += dual[1][i]*v
    inverse = {m: v for m, v in inverse.items() if v}
    half = vertices//2
    # The inverse form and both copies of the quartic have stored denominators 4,2,2.
    raised = {0: 1}
    for _ in range(half):
        raised = multiply(raised, inverse)
    coefficient = Fraction(raised.get(top, 0), 16**half)
    colours = endpoint["doublets"]
    b = {(colours.index(i), colours.index(j)): 2*v for i,j,v in endpoint["alternating_form"]}
    pf = orientation(pairs)*prod(b[p] for p in pairs)
    power = 7-half
    predicted = Fraction((-1)**half * factorial(power)*pf, factorial(half))*coefficient
    direct = 0
    for chosen in combinations(pairs, half):
        left = [p for p in pairs if p not in chosen]
        sites = sorted(i for p in chosen for i in p)
        sign = orientation(left + [tuple(sites)])
        rows = {0: 1}
        for i in sites:
            rows = multiply(rows, currents[i])
        direct += sign*prod(b[p] for p in left)*rows.get(top, 0)
    exact = Fraction((-1)**half * factorial(power)*direct, 2**vertices)
    assert predicted == exact
    return dict(sites_in_each_remaining_flavour=[colours[i] for i in slots],
                vertices=vertices, inverse_form_power=str(coefficient), direct_coefficient=str(exact))


def orientation(groups):
    order = [i for pair in groups for i in pair]
    return (-1)**sum(i > j for a, i in enumerate(order) for j in order[a+1:])


def embedded_su3(c, endpoint, tensor):
    low = endpoint["opposite_colour"]
    e = endpoint["raising_colour"]
    f = {low: Fraction(endpoint["lowering_coefficient"])}
    u, v = endpoint["alternating_form"][0][:2]
    assert set(c.br[u][v]) == {e}
    scale = Fraction(1, c.br[u][v][e])

    def bracket(x, y):
        out = Counter()
        for i, a in x.items():
            for j, b in y.items():
                for k, z in c.br[i][j].items():
                    out[k] += a*b*z
        return {k: v for k, v in out.items() if v}

    unit = {}
    for i, j in permutations(range(3), 2):
        m = np.zeros((3,3), dtype=np.int64)
        m[i,j] = 1
        unit[i,j] = m
    native = {(0,1): {u: Fraction(1)}, (1,2): {v: scale}, (0,2): {e: Fraction(1)}, (2,0): f}
    native[1,0] = bracket(native[1,2], f)
    native[2,1] = bracket(f, native[0,1])
    native[0,0] = bracket(native[0,1], native[1,0])
    native[1,1] = bracket(native[1,2], native[2,1])
    unit[0,0] = np.diag([1,-1,0])
    unit[1,1] = np.diag([0,1,-1])
    for key, a in unit.items():
        for other, b in unit.items():
            mat = a@b-b@a
            coordinates = {p: int(mat[p]) for p in permutations(range(3), 2)}
            coordinates.update({(0,0): int(mat[0,0]), (1,1): -int(mat[2,2])})
            expanded = Counter()
            for p, n in coordinates.items():
                for i, val in native[p].items():
                    expanded[i] += n*val
            assert {i:v for i,v in expanded.items() if v} == bracket(native[key], native[other])
            inner = sum(x*y*c.metric.get((i,j),0) for i,x in native[key].items()
                        for j,y in native[other].items())
            assert inner == int(np.trace(a@b))

    sites = list(map(endpoint["doublets"].index, (u,v)))
    entries = []
    for axes in product((0,1), repeat=4):
        indices = tuple(sites[i] for i in axes)
        value = Fraction(int(tensor[indices]),2) * scale**sum(axes)
        assert value == int(sum(axes) == 2)
        if value:
            entries.append(axes)
    adhm = Counter()
    for p in permutations(range(4)):
        # The epsilon sign and the Grassmann sign from sorting by flavour coincide.
        labels = tuple(int(p.index(i) in (1,3)) for i in range(4))
        adhm[labels] += orientation([p])**2
    assert sorted(adhm) == entries and set(adhm.values()) == {4}
    return dict(matrix_units={"E12":u,"E23":v,"E23_native_factor":str(scale),"E13":e,
                              "E31":low,"E31_native_factor":str(f[low])},
                tested_brackets=64, tested_metric_entries=64, restricted_metric="Tr3",
                native_diagonal_quartic="6*x^2*y^2", ordered_tensor_entries=6,
                SU3_ADHM_epsilon_coefficient_per_entry=4)


def direct_check(endpoint, result):
    import os
    from pathlib import Path
    import subprocess
    import tempfile

    colours = endpoint["doublets"]
    lines = [f"{colours.index(i)} {colours.index(j)} {2*v}" for i,j,v in endpoint["alternating_form"]]
    components = result["nonzero_native_components"]
    lines.append(str(len(components)))
    for *indices, value in components:
        twice = 2*Fraction(value)
        assert twice.denominator == 1
        lines.append(" ".join(map(str,[*(colours.index(i) for i in indices), int(twice)])))
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp)/"quartic"
        source = Path(__file__).with_name("fifth_quartic_direct.cpp")
        subprocess.run([os.environ.get("CXX","c++"), "-O3", "-std=c++17", str(source), "-o", str(exe)],check=True)
        output = subprocess.check_output([str(exe)],input="\n".join(lines)+"\n",text=True)
    direct = []
    for line in output.splitlines():
        m, power, tested, nonzero, top = map(int, line.split())
        assert tested == comb(7,m//2)**4
        direct.append(dict(vertices=m, power=power, colour_masks=tested, nonzero=nonzero, top_coefficient=str(top)))
    assert [r["vertices"] for r in direct] == [2,4,6]
    for row, reference in zip(direct,result["saturated_four_flavour_moments"]):
        assert row["top_coefficient"] == reference["top_coefficient"]
    return direct


def calculate(direct=False):
    endpoint = instanton()
    c = Colour()
    tensor, components = native_quartic(c, endpoint)
    dual = inverse_pairs(endpoint)
    cols = endpoint["doublets"]
    pairs = [(cols.index(i), cols.index(j)) for i, j, _ in endpoint["alternating_form"]]
    pfaffian = orientation(pairs)*prod(2*b for _, _, b in endpoint["alternating_form"])
    assert pfaffian == -1024
    # Two independent rooted subspaces, for each residual Grassmann degree.
    probes = []
    for m in (4, 6):
        for chosen in (pairs[:m//2], pairs[-m//2:]):
            slots = sorted(i for pair in chosen for i in pair)
            probes.append(restricted_top(tensor, dual, endpoint, slots, m))
    full = [saturated(tensor, dual, pfaffian, k) for k in (2, 4)]
    result = dict(symmetry=invariant_dimensions(), tensor_convention="polarization of B(f,(ad v)^4 f), B=Tr26/6",
                symmetric_coefficients_tested=comb(17,4), nonzero_native_components=components,
                ordered_nonzero=int(np.count_nonzero(tensor)), pfaffian=pfaffian,
                restricted_coefficients=probes, saturated_four_flavour_moments=full,
                embedded_SU3=embedded_su3(c, endpoint, tensor),
                physical_four_fermion_coefficient="nonzero on embedded SU3; native collective-coordinate scale not matched",
                physical_instanton_supercharge_pairing="open")
    if direct:
        result["direct_four_flavour_exterior"] = direct_check(endpoint, result)
    return result


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct", action="store_true", help="also rebuild the independent four-flavour exterior contraction")
    print(json.dumps(calculate(parser.parse_args().direct), indent=2))
