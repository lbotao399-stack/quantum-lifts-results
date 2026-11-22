"""Contract the 22 source colour monomials with three one-instanton scalar propagators."""

from itertools import permutations
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import json
from cubic_colours import load
from source_fluctuations import wedge, add, scale, matrix_product as mm, trace


def mode(a, j):
    """One Grassmann generator in flavour a and its local six-mode basis."""
    return {1 << (6 * a + j): 1}


def scalar(a, b):
    def col(g):
        return [scale(mode(g, 0), 4), scale(mode(g, 1), 4), mode(g, 4)]

    def row(g):
        return [scale(mode(g, 1), -4), scale(mode(g, 0), 4), mode(g, 5)]

    x, y, u, v = col(a), col(b), row(a), row(b)
    mat = [
        [add(wedge(y[i], u[j]), scale(wedge(x[i], v[j]), -1)) for j in range(3)]
        for i in range(3)
    ]
    # Twice the common dyad-normalized scalar; this keeps all mode coefficients integral.
    corr = scale(trace(mat), -1)
    mat = [[scale(z, 2) for z in row] for row in mat]
    for i in range(2):
        mat[i][i] = add(mat[i][i], corr)
    assert not trace(mat)
    return mat


def differentiate_profile(poly, spin, dot):
    ret = Counter()
    for word, c in poly.items():
        for flavour in range(4):
            for a, u in enumerate(spin):
                for d, b in enumerate(dot):
                    if not (u * b):
                        continue
                    i, j = 6 * flavour + a, 6 * flavour + 2 + d
                    if word & (1 << i) and not word & (1 << j):
                        change = (-1) ** ((word & ((1 << j) - (2 << i))).bit_count())
                        ret[word ^ (1 << i) ^ (1 << j)] -= c * change * u * b
    return {k: v for k, v in ret.items() if v}


def trace_product(a, b):
    return add(*(wedge(a[i][j], b[j][i]) for i in range(3) for j in range(3)))


def exterior_product(*p):
    z = {0: 1}
    for v in p:
        z = wedge(z, v)
    return z


def linear_combination(*p):
    return add(*(scale(x, k) for k, x in p))


def source_modes(mask):
    v = {0: 1}
    for i in range(12):
        if mask >> i & 1:
            a, j = divmod(i, 4)
            factor = (
                linear_combination((4, mode(a, j)), (-3, mode(a, j + 2)))
                if j < 2
                else scale(mode(a, j + 2), 5)
            )
            v = wedge(v, factor)
    return v


def transported_dual_basis():
    # Common factor of 150 on each transported matrix for integer arithmetic.
    out = []
    for i in range(3):
        for j in range(3):
            if i == j == 2:
                continue
            x = [
                [
                    Fraction(int(a == j and b == i))
                    - (Fraction(1, 3) if i == j and a == b else 0)
                    for b in range(3)
                ]
                for a in range(3)
            ]
            m = (Fraction(4, 5), Fraction(4, 5), Fraction(1))
            y = [[m[a] * x[a][b] * m[b] for b in range(3)] for a in range(3)]
            ty = sum(y[a][a] for a in range(3))
            tx = Fraction(16, 25) * (x[0][0] + x[1][1])
            for a in range(3):
                y[a][a] += -ty / 3 + Fraction(9, 32) * tx * (
                    int(a < 2) - Fraction(2, 3)
                )
            # Independent decomposition into the embedded triplet, two doublets and the singlet.
            mean = (x[0][0] + x[1][1]) / 2
            for k in range(2):
                for z in range(2):
                    singlet = mean if k == z else 0
                    assert y[k][z] == m[0] ** 2 * (x[k][z] - singlet) + singlet
                assert y[k][2] == m[0] * x[k][2] and y[2][k] == m[0] * x[2][k]
            assert y[2][2] == x[2][2] and not sum(y[a][a] for a in range(3))
            out.append([[{0: int(v * 150)} if v else {} for v in row] for row in y])
            assert all((150 * v).denominator == 1 for row in y for v in row)
    return out


class RadialProjection:
    def __init__(self, spin=(1, 0), dots=((1, 0), (0, 1))):
        self.A = scalar(0, 3)
        s, t = scalar(1, 3), scalar(2, 3)
        self.B, self.C = (
            [[differentiate_profile(v, spin, dot) for v in row] for row in mat]
            for dot, mat in zip(dots, (s, t))
        )
        self.U = transported_dual_basis()
        self.fields = self.U + [self.A]
        b, c = self.B, self.C
        self.b2 = trace_product(b, b)
        self.c2 = trace_product(c, c)
        self.bc = trace_product(b, c)
        self.disc = linear_combination(
            (1, exterior_product(self.b2, self.c2)),
            (-1, exterior_product(self.bc, self.bc)),
        )

    @lru_cache(None)
    def gram(self, k, second):
        x, y = self.fields[k], self.fields[second]
        b, c = self.B, self.C
        # 15 times symmetric polarization of the quadratic bracket of the five words.
        return linear_combination(
            (24, exterior_product(trace_product(x, y), self.disc)),
            (32, exterior_product(trace_product(x, b), trace_product(y, b), self.c2)),
            (-32, exterior_product(trace_product(x, b), trace_product(y, c), self.bc)),
            (-32, exterior_product(trace_product(x, c), trace_product(y, b), self.bc)),
            (32, exterior_product(trace_product(x, c), trace_product(y, c), self.b2)),
        )

    @lru_cache(None)
    def third(self, k, second, n):
        x, y, z = (self.fields[q] for q in (k, second, n))
        return add(trace_product(mm(x, y), z), trace_product(mm(y, x), z))

    def radial(self, key):
        # Split three labelled derivatives between the cubic trace and the quadratic bracket.
        terms = []
        A = 8
        for i in range(3):
            j, second = [key[t] for t in range(3) if t != i]
            terms.append(
                (
                    6,
                    exterior_product(
                        trace_product(mm(self.A, self.A), self.U[key[i]]),
                        self.gram(j, second),
                    ),
                )
            )
            terms.append(
                (6, exterior_product(self.third(A, j, second), self.gram(A, key[i])))
            )
        terms.append((3, exterior_product(self.third(*key), self.gram(A, A))))
        return linear_combination(*terms)

    def literal(self, key):
        fields = dict(
            A=self.A, B=self.B, C=self.C, **{str(i): v for i, v in enumerate(self.U)}
        )

        @lru_cache(None)
        def traced(word):
            mat = fields[word[0]]
            for c in word[1:]:
                mat = mm(mat, fields[c])
            return trace(mat)

        total = Counter()
        terms = json.loads(Path(__file__).with_name("result.json").read_text())[
            "three_jet_traces"
        ]
        for expression, coefficient in terms:
            factors = expression.split()
            slots = [
                (i, j)
                for i, w in enumerate(factors)
                for j, c in enumerate(w)
                if c == "A"
            ]
            assert len(slots) == 5
            for triple in permutations(slots, 3):
                change = {slot: str(c) for slot, c in zip(triple, key)}
                polynomials = [
                    traced(tuple(change.get((i, j), v) for j, v in enumerate(w)))
                    for i, w in enumerate(factors)
                ]
                total.update(
                    scale(exterior_product(*polynomials), 15 * Fraction(coefficient))
                )
        return {k: v for k, v in total.items() if v}

    def run(self, data, literal):
        joined = Counter()
        independent = Counter()
        intermediate = Counter()
        for key, coeffs in data.items():
            targ = self.radial(key)
            if literal:
                assert targ == self.literal(key), key
            assert all(z.bit_count() == 12 and (z >> 18).bit_count() == 6 for z in targ)
            s = add(*(scale(source_modes(m), v) for m, v in coeffs.items()))
            sub = wedge(s, targ)
            intermediate[key] = len(targ)
            independent[key] = len(sub)
            joined.update(sub)
        result = {k: v for k, v in joined.items() if v}
        assert all(
            (m >> 18) == 63
            and all(((m >> (6 * a)) & 63).bit_count() == 5 for a in range(3))
            for m in result
        )
        return result, independent, intermediate


def verify(literal=True):
    data = load()
    result = {}
    for spin in ((1, 0), (0, 1), (1, 1)):
        probe = RadialProjection(spin)
        total, joined, radial = probe.run(data, literal)
        # Remove the integer-arithmetic factors tracked in scalar(), transported_dual_basis(), gram()
        # and source_modes(). Physical common propagator/field factors are not included.
        cleared = 15 * 2**6 * 150**3 * 5**9
        result[str(spin)] = dict(
            radial_colour_monomials=len(+radial),
            paired_colour_monomials=len(+joined),
            exterior_coefficients={
                hex(k): str(Fraction(v, cleared)) for k, v in sorted(total.items())
            },
        )
    assert [len(d["exterior_coefficients"]) for d in result.values()] == [14, 0, 18]
    first = result[str((1, 0))]["exterior_coefficients"]
    assert Fraction(first["0xfdfdef"]) / Fraction(first["0xfefddf"]) == -1
    assert Fraction(first["0xffd7ef"]) / Fraction(first["0xfdfdef"]) == -Fraction(3, 4)
    return result


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
