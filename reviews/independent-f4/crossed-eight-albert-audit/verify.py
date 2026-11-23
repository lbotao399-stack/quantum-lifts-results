"""Replay the rational crossed-graph scores and the two order-704 determinants."""
from argparse import ArgumentParser
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NORMALIZATION = 2 * Fraction(28, 3)**4


def residue(value, prime):
    value = Fraction(value)
    return value.numerator * pow(value.denominator, -1, prime) % prime


def leaves(tree):
    if isinstance(tree, int):
        yield tree
    else:
        for child in tree:
            yield from leaves(child)


def order_sign(graph, bits, tree):
    plus, minus = graph['cp'], graph['cm']
    occupied = sorted(set(plus) | set(minus))
    owners = {edge: graph['E'][edge][bit] for edge, bit in zip(occupied, bits, strict=True)}
    visited = [edge for v in leaves(tree) for edge in occupied if owners[edge] == v]
    assert sorted(visited) == occupied
    inversions = 0
    for listed in (plus, minus):
        positions = [listed.index(edge) for edge in visited if edge in listed]
        inversions += sum(a > b for i, a in enumerate(positions) for b in positions[i+1:])
    return (-1)**inversions


def determinant(matrix, p):
    a = np.array(matrix, dtype=np.int64, copy=True) % p
    det = 1
    for i in range(len(a)):
        nonzero = np.flatnonzero(a[i:, i])
        if len(nonzero) == 0:
            return 0
        j = i + int(nonzero[0])
        if j != i:
            a[[i, j]] = a[[j, i]]
            det = -det
        pivot = int(a[i, i])
        det = det * pivot % p
        multiplier = a[i+1:, i] * pow(pivot, -1, p) % p
        a[i+1:, i+1:] = (a[i+1:, i+1:] - multiplier[:, None] * a[i, i+1:]) % p
    return det % p


def weighted_pair(graph, u, v):
    return graph['w'] * (-4*2392*7113106*u + 2*2299*1850212*v)


def read_poles(graphs):
    poles = {'left': {}, 'right': {}}
    totals = {}
    for line in (HERE / 'exact-poles.jsonl').read_text().splitlines():
        entry = json.loads(line)
        side = entry['side']
        if 'graph' in entry:
            i = entry['graph']
            assert i not in poles[side]
            poles[side][i] = list(map(Fraction, entry['poles']))
        else:
            assert side not in totals
            totals[side] = [Fraction(entry['raw']), Fraction(entry['unsubtracted_score'])]
    assert set(poles['left']) == set(range(18)) and set(poles['right']) == set(range(12,30))
    assert graphs['left'][:12] == graphs['right'][:12]
    poles['right'].update({i: poles['left'][i] for i in range(12)})
    return poles, totals


def rational_scores(graphs, poles, recorded):
    result = {}
    for side, bank in graphs.items():
        terms = [weighted_pair(graph, *poles[side][i]) for i, graph in enumerate(bank)]
        signs = [order_sign(graph, graph['bits'], graph['tree']) for graph in bank]
        raw = sum(terms)
        tree = sum(t*s for t, s in zip(terms, signs, strict=True)) / NORMALIZATION
        assert [raw, raw / NORMALIZATION] == recorded[side]
        assert raw == -92613596121446400 and tree == Fraction(-11362589638943400,49)
        # Only the two children of graph 0's right subtree are interchanged.
        g = bank[0]
        swapped = [g['tree'][0], list(reversed(g['tree'][1]))]
        assert signs[0] == -1 and order_sign(g, g['bits'], swapped) == 1
        difference = 2*terms[0]/NORMALIZATION
        assert difference == Fraction(7774403437171800,49)
        result[side] = (tree, raw / NORMALIZATION)
        print(f'{side}  tree-order score {tree}  graph-order score {result[side][1]}  '
              f'tree-order change on swapping graph 0 {difference}')
    return result


def check_observations(graphs, poles, observations, p):
    for side, entries in observations.items():
        for key, observation in entries.items():
            i = int(key)
            graph = graphs[side][i]
            sign = order_sign(graph, observation['bits'], observation['tree'])
            u, v = observation['poles']
            assert weighted_pair(graph, u, v) % p == observation['weight']
            assert [residue(t,p) for t in poles[side][i]] == [sign*u % p, sign*v % p]


def check_section(section):
    data = json.loads((HERE / 'section-graphs.json').read_bytes())
    rows = [json.loads(line) for line in (HERE / 'section-poles.jsonl').read_text().splitlines()]
    pairs = {r['graph']: tuple(map(Fraction,r['poles'])) for r in rows if 'graph' in r}
    expected = {r['pivot']: [r['left'],r['right']] for r in rows if 'pivot' in r}
    assert set(pairs) == set(range(len(data['graphs']))) and set(expected) == {700,701,702}
    weights = {i: -4*2392*7113106*u + 2*2299*1850212*v for i,(u,v) in pairs.items()}
    for row, banks in data['rows'].items():
        both = [sum((w*weights[i] for i,w in bank), Fraction()) / (2*Fraction(28,3)**3)
                for bank in banks]
        assert list(map(str,both)) == expected[int(row)]
        assert both == [Fraction(section.get(row,0))]*2
    print('six-cubic pivots 700, 701, 702:', *(section.get(str(i),'0') for i in (700,701,702)))


def main():
    cli = ArgumentParser(description=__doc__)
    cli.add_argument('--classical', type=Path, default=HERE.parent / 'sixth-glueball-classical-rank')
    args = cli.parse_args()
    data = json.loads((HERE / 'audit-data.json').read_bytes())
    graphs = data['graphs']
    poles, totals = read_poles(graphs)
    scores = rational_scores(graphs, poles, totals)
    check_section(data['section'])
    for p, base, old, projected in ((1009,1005,551,44), (1013,734,478,92)):
        check_observations(graphs,poles,data['observations'][str(p)],p)
        coord = data['coordinates'][str(p)]
        packed = (args.classical / f'selected-{p}.u16.gz').read_bytes()
        assert sha256(packed).hexdigest() == coord['chart_sha256']
        m = np.frombuffer(gzip.decompress(packed), dtype='<u2').reshape(703,703).astype(np.int64)
        c, y = np.asarray(coord['coefficients']), np.asarray(coord['outgoing'])
        assert determinant(m, p) == base and np.array_equal(c @ m % p, y)
        s = np.zeros(703, dtype=np.int64)
        for i, v in data['section'].items():
            s[int(i)] = residue(v,p)
        section = int(c @ s) % p
        assert section == projected
        for side, pair in scores.items():
            row = [residue(x,p) for x in pair]
            determinants = []
            for last in row:
                augmented = np.block([[m, s[:,None]], [y[None,:], np.array([[last]])]])
                d = determinant(augmented,p)
                assert d == base*(last-section) % p
                determinants.append(d)
            assert determinants == [old,0] and row[1] == section
            print(f'{p} {side}  classical {base}  tree {row[0]}  ordered {row[1]}  '
                  f'section {section}  minors {determinants}')


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with Python assertions enabled')
    main()
