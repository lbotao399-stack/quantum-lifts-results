"""Rebuild the two extremal contributions to the F4 scalar endpoint over Q."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
import json
from math import prod
from pathlib import Path
from zipfile import ZipFile

from first_action import FirstAction, inputs
from native import matrix
from verify import ARCHIVE, contents
from ward_boundary import Boundary, graph, residue


def sign(order):
    return (-1)**sum(a > b for i, a in enumerate(order) for b in order[i+1:])


class ScalarLift:
    def __init__(self, colour):
        self.bracket = {}
        for i, j, terms in colour['bracketPairs']:
            self.bracket[i, j] = {k: int(v) for k, v in terms}
            self.bracket[j, i] = {k: -int(v) for k, v in terms}
        self.inverse = matrix(colour['inverseMetric'])
        self.weight = [tuple(self.bracket.get((i,j), {}).get(j, 0) for i in range(4))
                       for j in range(52)]
        self.pairs = self.by_weight(combinations_with_replacement(range(52), 2))
        self.triples = self.by_weight(combinations(range(52), 3))

    def by_weight(self, entries):
        groups = defaultdict(list)
        for key in entries:
            weight = tuple(map(sum, zip(*(self.weight[j] for j in key))))
            groups[weight].append(key)
        return groups

    @cache
    def omega(self, a, b):
        result = Counter()
        for (i,j), c in self.inverse.items():
            for u, x in self.bracket.get((i,a), {}).items():
                for v, y in self.bracket.get((j,b), {}).items():
                    result[u,v] += c*x*y
        return result

    @cache
    def casimir(self, a, b):
        result = Counter({(a,b): 36})
        for pair, c in self.omega(a,b).items():
            result[tuple(sorted(pair))] += 2*c
        return result

    @cache
    def bracket_triple(self, triple):
        a,b,c = triple
        result = Counter()
        for u,v,w in ((a,b,c), (b,c,a), (c,a,b)):
            for t, value in self.bracket.get((u,v), {}).items():
                result[tuple(sorted((t,w)))] += value
        return result

    def scalar_values(self, u, v):
        weight = tuple(map(sum, zip(self.weight[u], self.weight[v])))
        pairs = self.pairs[weight]
        row = {ab: (self.omega(*ab)[u,v] + self.omega(*ab)[v,u])/Fraction(2)
               for ab in pairs}
        # The polynomial takes the inverse eigenvalue of k^t k at each Casimir.
        eigen = {0: -27, 26: -14, 40: -7}
        projected = Counter()
        for c, value in eigen.items():
            component = row.copy()
            for h in eigen.keys() - {c}:
                component = {ab: (sum(x*component[q] for q, x in self.casimir(*ab).items())
                                  - h*component[ab])/(c-h) for ab in pairs}
            for ab, x in component.items():
                projected[ab] += x/value
        return {triple: sum(x*projected[ab] for ab,x in self.bracket_triple(triple).items())
                for triple in self.triples[weight]}


def branches(colour, normalization):
    lift = ScalarLift(colour)
    metric = matrix(colour['metric'])
    highest = [(19,29,31,45,46,47), (19,20,29,45,46,47)]
    dimensions = [r['dimension'] for r in normalization['branches']]
    harmonic = normalization['weights']
    for u, dim, scalar in zip(highest, dimensions, harmonic, strict=True):
        v = tuple(next(j for x,j in metric if x == i) for i in u)
        wedge_metric = prod(metric[i,j] for i,j in zip(u,v))
        rows = Counter()
        for i, a in enumerate(u):
            for j, b in enumerate(sorted(v)):
                rest = (u[:i]+u[i+1:] + tuple(x for x in sorted(v) if x != b))
                factor = scalar * dim * sign(v) * (-1)**(i+j)
                for triple, coefficient in lift.scalar_values(a,b).items():
                    if coefficient:
                        rows[rest+triple] = factor*coefficient
        yield wedge_metric, rows


def read(z, name):
    return json.loads(contents(z, name))


def compute(z):
    colour = read(z, 'collected-O0/color-basis-native.json')
    normalization = read(z, 'independent-audit/target-normalization.json')
    wedge_rows = list(branches(colour, normalization))
    weights = read(z, 'exact-score/eta0-folded-exact.json')
    orbit = {tuple(k): (tuple(f), Fraction(s)) for k,(f,s) in weights['orbit']}
    folded = {tuple(k): Fraction(x) for k,x in weights['weights']}
    first_prime = 1000003
    unfolded = dict(graph(z, 'eta0-1000003-unfolded.in', first_prime))
    original = Counter()
    columns = []
    for metric, row in wedge_rows:
        folded_branch = Counter()
        for k,c in row.items():
            original[k] += c
            f,s = orbit.get(k, (k, 1))
            folded_branch[f] += s*c
        columns.append(folded_branch)
    assert {k: residue(v, first_prime) for k,v in original.items() if v} == unfolded
    assert {k:v for k in set().union(*columns) if (v := sum(r[k] for r in columns))} == folded

    boundary = Boundary(colour)
    b = {tuple(k): Fraction(x)/6 for k,x in read(z, 'exact-score/eta0-complete-exact-anchors.json')}
    b.update({tuple(k): Fraction(x)/6 for k,x in read(z,'cochains/inputs/exact-raw-anchors.json')['anchors']})
    y = {k:boundary.coefficient(k) for k in sorted(b)}
    assert set(y) == set(b)
    digest = sha256(json.dumps([[list(k), str(y[k])] for k in sorted(y)],separators=(',', ':')).encode()).hexdigest()
    assert digest == 'bd35f8393b86358cc4f329b3ee9a9b3f51b3267b620be85a31000a0c3cc96fbd'
    extra = sorted(set().union(*columns) - b.keys())
    assert len(extra) == 4 and all(original[k] == 0 for k in extra)
    tables, source, _ = inputs(z)
    action = FirstAction(tables, source)
    for k in extra:
        b[k] = -action.coefficient(k)
        y[k] = boundary.coefficient(k)
    target = read(z,'exact-score/exact_normalized_class.json')
    details = {}
    for name, cochain in (('Q0Y',y), ('archive_raw_anchor_over_6',b)):
        separate = [sum(x*cochain[k] for k,x in row.items()) for row in columns]
        scalar = sum(m*v for (m,_),v in zip(wedge_rows,separate))
        tests = []
        for prime in (1000003, 1000033):
            endpoint = graph(z, f'graph-{prime}.in', prime)
            k_value = sum(w*residue(cochain[k],prime) for k,w in endpoint) % prime
            tests.append(dict(prime=prime, etaK=k_value,
                              corrected_sum=(residue(scalar,prime)+k_value)%prime))
        entry = dict(branches=list(map(str,separate)), original_eta0=str(sum(separate)),
                     corrected_eta0=str(scalar), reductions=tests)
        if name == 'archive_raw_anchor_over_6':
            assert sum(separate) == -Fraction(target['eta0_corrected'])
            eta_k = -Fraction(target['etaK_corrected'])
            assert all(residue(eta_k, t['prime']) == t['etaK'] for t in tests)
            entry.update(etaK=str(eta_k), corrected_sum=str(scalar+eta_k))
        details[name] = entry
    mixed = list(map(Fraction, normalization['row']))
    harmonic = normalization['weights']
    physical = sum(c*x for c,x in zip(harmonic,mixed,strict=True))
    exported = sum(c*x/m for c,x,(m,_) in zip(harmonic,mixed,wedge_rows,strict=True))
    assert physical == 0 and exported == Fraction(38456,35)
    assert details['archive_raw_anchor_over_6']['corrected_sum']=='0'
    return dict(archive_sha256=ARCHIVE, wedge_metric=[str(m) for m,_ in wedge_rows],
                classical_mixed_row=list(map(str,mixed)),
                mixed_pairing_physical=str(physical),mixed_pairing_with_exported_weights=str(exported),
                unfolded_rows=len(unfolded), folded_rows=len(folded),
                newly_required_components=[[list(k), str(b[k]), str(y[k])] for k in extra],
                pairings=details)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    if not __debug__ or sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE:
        parser.error('use the pinned archive with Python assertions enabled')
    with ZipFile(args.archive) as archive:
        print(json.dumps(compute(archive), indent=2))
