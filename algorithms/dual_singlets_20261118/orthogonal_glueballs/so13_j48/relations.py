"""Verify the ordered-word minor and the rational second pairing in SO(13)."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from flint import fmpz, nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parent))
from orthogonal_glueball_relation_span import determinant, relation_dimension
from orthogonal_glueball_second_pair import second_dual
from so11_j36.verify import check_adjoints
from terminal_commutators import words


def initialize(probes, scores, sites, modulus):
    global ORDERED, SCORES, SITES, MODULUS
    ORDERED = [tuple(map(tuple, w)) for w in probes]
    SCORES = [(tuple(map(tuple, w)), Fraction(c)) for w, c in scores]
    SITES = sites
    MODULUS = modulus


def check_row(row):
    relation = determinant(8, row['spec'])
    assert relation.pf_degree > 13
    values = [relation.coefficient(w) for w in ORDERED]
    residues = [v.numerator * pow(v.denominator, -1, MODULUS) % MODULUS for v in values]
    assert residues == row['word_residues']
    sites = [relation.cartan(site) for site in SITES]
    assert sites == list(map(Fraction, row['sites']))
    score = sum(c * relation.coefficient(w) for w, c in SCORES)
    assert score == Fraction(row['score'])
    relation.clear()
    return row['index']


def verify(rebuild, workers):
    data = json.loads(gzip.decompress((HERE / 'relations.json.gz').read_bytes()))
    rows = data['rows']
    prime = data['prime']
    assert prime == 1000003 and fmpz(prime).is_prime()
    dimension = relation_dimension(13, 8, 2)
    assert dimension == 301
    assert len({r['index'] for r in rows}) == len(rows) == len(data['ordered_probes']) <= dimension
    minor = nmod_mat([r['word_residues'] for r in rows], prime).det()
    assert int(minor) == data['minor_determinant'] != 0
    for row in rows:
        first, second = map(Fraction, row['sites'])
        assert Fraction(row['score']) + Fraction(9, 1280)*first + Fraction(1, 80)*second == 0
    dual = second_dual(13, 8)
    assert dual['first_factorization_columns'] == 33 and dual['second_factorization_columns'] == 1266
    assert len(dual['scores']) == 251
    assert check_adjoints(13, 8, dual) == [39, 1362]
    if rebuild:
        sites = words((6, 6, 2, 2, 2))
        with ProcessPoolExecutor(workers, initializer=initialize,
                                 initargs=(data['ordered_probes'], dual['scores'],
                                           (sites[1], sites[2]), prime)) as pool:
            for i, row in zip(pool.map(check_row, rows), rows, strict=True):
                assert i == row['index']
    return dict(balanced_relations=dimension, independently_tested_relations=len(rows),
                minor_modulo_1000003=int(minor), nonzero_rational_residuals=0,
                all_relations_tested=(len(rows) == dimension))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild', action='store_true', help='regenerate the Pfaffian coefficients')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error('positive worker count and Python assertions are required')
    print(json.dumps(verify(args.rebuild, args.workers)))
