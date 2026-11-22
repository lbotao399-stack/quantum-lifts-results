"""Test the archived F4 first adjoint on a six-letter trace times S^3."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache, lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import factorial
from pathlib import Path
from zipfile import ZipFile

from verify import ARCHIVE, contents


class Boundary:
    def __init__(self, colour):
        self.metric = {(i, j): Fraction(x) for i, j, x in colour['inverseMetric']}
        self.generators = []
        for generator in colour['nativeDualTranspose26']:
            entries = {}
            for i, j, value in generator:
                scaled = 2 * Fraction(value)
                assert scaled.denominator == 1
                entries[i, j] = int(scaled)
            self.generators.append(entries)
        self.rows = []
        for generator in self.generators:
            rows = defaultdict(list)
            for (i, j), value in generator.items():
                rows[i].append((j, value))
            self.rows.append(rows)

    def times_generator(self, left, index):
        result = defaultdict(int)
        rows = self.rows[index]
        for (i, k), a in left.items():
            for j, b in rows[k]:
                result[i, j] += a * b
        return {k: v for k, v in result.items() if v}

    @cache
    def bracket(self, i, j):
        left = self.times_generator(self.generators[i], j)
        right = self.times_generator(self.generators[j], i)
        return {k: d for k in left.keys() | right.keys()
                if (d := left.get(k, 0) - right.get(k, 0))}

    @lru_cache(maxsize=24000)
    def left(self, scalar, ghost, other_scalar, other_ghost):
        return self.times_generator(self.times_generator(self.bracket(scalar, ghost), other_scalar), other_ghost)

    @lru_cache(maxsize=12000)
    def right(self, scalar, i, j):
        first = self.times_generator(self.times_generator(self.generators[scalar], i), j)
        second = self.times_generator(self.times_generator(self.generators[scalar], j), i)
        return {(b, a): x for a, b in first.keys() | second.keys()
                if (x := first.get((a, b), 0) + second.get((a, b), 0))}

    @cache
    def minor(self, plus, minus):
        if not plus:
            return 1
        return sum((-1)**j * self.metric[plus[0], b] * self.minor(plus[1:], minus[:j] + minus[j+1:])
                   for j, b in enumerate(minus) if (plus[0], b) in self.metric)

    @cache
    def core(self, plus, minus, scalars):
        # Both epsilon signs cancel the Grassmann sign on exchanging that pair.
        value = 0
        for i, j in product(range(2), repeat=2):
            orientation = (-1)**(1 + i + j)
            u, v = plus[i], minus[j]
            for p, q, r in permutations(scalars):
                sign = (-1)**((p > q) + (p > r) + (q > r))
                right = self.right(r, plus[1-i], minus[1-j])
                if right:
                    for a, b in ((u, v), (v, u)):
                        left = self.left(p, a, q, b)
                        value += orientation * sign * sum(x * right.get(k, 0) for k, x in left.items())
        return Fraction(value, 2**7)

    @cache
    def coefficient(self, key):
        plus, minus, scalars = key[:5], key[5:10], key[10:]
        assert all(values == tuple(sorted(set(values))) for values in (plus, minus, scalars))
        value = Fraction()
        for u, us, up in split(plus):
            for v, vs, vp in split(minus):
                m = self.minor(us, vs)
                if m:
                    value -= factorial(3) * up * vp * m * self.core(u, v, scalars)
        return value


@cache
def split(colours):
    return [(tuple(colours[i] for i in selected), tuple(c for i, c in enumerate(colours) if i not in selected),
             (-1)**(sum(selected)-1)) for selected in combinations(range(5), 2)]


def residue(x, prime):
    x = Fraction(x)
    return x.numerator * pow(x.denominator, -1, prime) % prime


def graph(archive, name, prime):
    rows = contents(archive, 'full-score/' + name).decode().splitlines()
    actual_prime, metric_terms = map(int, rows[0].split())
    assert actual_prime == prime and int(rows[1+metric_terms]) == len(rows)-2-metric_terms
    result = []
    for line in rows[2+metric_terms:]:
        weight, label, *key = map(int, line.split())
        if label < 0:
            result.append((tuple(key), weight))
        else:
            assert weight == 0
    return result


def evaluate(archive):
    def read(name):
        return json.loads(contents(archive, name))

    boundary = Boundary(read('collected-O0/color-basis-native.json'))
    folded = {tuple(k): Fraction(x) for k, x in read('exact-score/eta0-folded-exact.json')['weights']}
    raw_source = {tuple(k): Fraction(x)/6 for k, x in read('cochains/inputs/exact-raw-anchors.json')['anchors']}
    raw_source.update({tuple(k): Fraction(x)/6 for k, x in read('exact-score/eta0-complete-exact-anchors.json')})
    normalized = read('exact-score/exact_normalized_class.json')
    eta0 = sum(w * boundary.coefficient(k) for k, w in folded.items())
    assert sum(w * raw_source[k] for k, w in folded.items()) == -Fraction(normalized['eta0_corrected'])

    rows = read('cochains/inputs/exact-coordinate-rows.json')['rows']
    coordinates = [sum(Fraction(w) * boundary.coefficient(tuple(k)) for k, w in row) for row in rows]
    records = []
    for prime in (1000003, 1000033):
        tables = contents(archive, f'cochains/inputs/convert-{prime}.jsonl').splitlines()
        weights = [w for line in tables for w in json.loads(line)['etaK_weight']]
        assert len(weights) == len(coordinates) == 40
        k_value = sum(w * residue(x, prime) for w, x in zip(weights, coordinates)) % prime
        def contract(edges, coefficient):
            return sum(w * residue(coefficient(k), prime) for k, w in edges) % prime
        k_edges = graph(archive, f'graph-{prime}.in', prime)
        assert k_value == contract(k_edges, boundary.coefficient)
        assert contract(k_edges, raw_source.__getitem__) == residue(-Fraction(normalized['etaK_corrected']), prime)
        for name in ([f'eta0-{prime}.in', 'eta0-1000003-unfolded.in'] if prime == 1000003 else [f'eta0-{prime}.in']):
            edges = graph(archive, name, prime)
            assert contract(edges, boundary.coefficient) == residue(eta0, prime)
        total = (residue(eta0, prime) + k_value) % prime
        assert k_value == 10 * residue(eta0, prime) % prime and total != 0
        records.append(dict(prime=prime, eta0=residue(eta0, prime), etaK=k_value, sum=total))
    # A second evaluation with generic native tree splitting gives this digest.
    values = [[list(k), str(boundary.coefficient(k))] for k in sorted(raw_source)]
    digest = sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
    assert len(values) == 2981 and digest == 'bd35f8393b86358cc4f329b3ee9a9b3f51b3267b620be85a31000a0c3cc96fbd'
    assert eta0 == Fraction(-98441322408000, 7)
    return dict(archive_sha256=ARCHIVE, archived_source_score=normalized['d2_coefficient'],
                archived_component_count=len(values), nonzero_archived_components=sum(x[1] != '0' for x in values),
                archived_components_sha256=digest, eta0=str(eta0), distinct_queries=boundary.coefficient.cache_info().currsize,
                eta0_unfolded_rows_at_first_prime=7276, etaK_raw_rows=477, prime_tests=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    if not __debug__ or sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE:
        parser.error('use the pinned archive with Python assertions enabled')
    with ZipFile(args.archive) as source:
        print(json.dumps(evaluate(source), indent=2))
