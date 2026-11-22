"""Recover the rational eta_K functional and evaluate it on independent inputs."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from sympy import Matrix, Rational

from endpoint_metric import read
from gram_data import export
from verify import ARCHIVE, HERE, contents
from ward_boundary import Boundary, residue


def matrix(rows):
    return Matrix([[Rational(x) for x in row] for row in rows])


def vector(values):
    return Matrix([Rational(x) for x in values])


def rebuild_gram(z, compiler):
    with TemporaryDirectory() as directory:
        folder = Path(directory)
        with (folder/'gram.in').open('w') as f:
            export(z, f)
        subprocess.run([compiler, '-std=c++17', '-O3', '-ftrapv', str(HERE/'gram11.cpp'),
                        '-o', str(folder/'gram11')], check=True)
        with (folder/'gram.in').open() as f:
            return subprocess.check_output([str(folder/'gram11')], stdin=f, text=True)


def weights(z, raw_gram):
    location = 'cochains/inputs/'
    lines = raw_gram.splitlines()
    pairs, nonzero, denominator, scale = map(int, lines[0].split())
    assert (pairs, nonzero, denominator, scale) == (573459, 563792, 757350, -157339)
    # Each of the eight rotated dual slots used twice the inverse metric.
    last = matrix([line.split() for line in lines[1:]]) * Rational(scale, 2**8 * denominator**2)
    assert last.shape == (9, 9) and last == last.T
    first_raw = contents(z, location+'block00-QQ.log').decode().splitlines()[1].split()
    assert first_raw[0] == 'GRAM' and len(first_raw) == 1+18**2
    first = vector(first_raw[1:]).reshape(18, 18)
    middle = matrix(read(z, location+'block10-exact-physical-Ht.json')['physical_Gram'])
    grams = first, middle, last
    reductions = {p: [json.loads(line) for line in
                      contents(z, location+f'convert-{p}.jsonl').splitlines()]
                  for p in (1000003, 1000033)}
    exact = read(z, location+'exact-restriction-and-lift-coordinates.json')['blocks']
    target_coeff = vector(read(z, location+'f4-j36-dphi-exact-solve.json')['coeff'])
    anchor = read(z, 'exact-score/positive-O1-Dphi-exact.json')['block_coordinates']
    covector = []
    for i, (gram, frame, test) in enumerate(zip(grams, exact, anchor, strict=True)):
        left, right = frame['a'], frame['b']
        restriction = matrix(frame['E'])
        assert gram == gram.T and gram * vector(test['t']) == vector(test['physical_Ht'])
        action = gram * target_coeff[left:right, :]
        dual = restriction.T.inv() * action
        assert restriction.T * dual == action
        for prime, blocks in reductions.items():
            row = blocks[i]
            assert row['block'] == frame['typ'] and row['dim'] == right-left
            assert gram.applyfunc(lambda x: residue(Fraction(x), prime)) == Matrix(row['H'])
            assert restriction.applyfunc(lambda x: residue(Fraction(x), prime)) == Matrix(row['E'])
            assert [residue(Fraction(x), prime) for x in dual] == row['etaK_weight']
        covector.extend(Fraction(x) for x in dual)
    assert len(covector) == 40
    return covector


def evaluate(z, raw_gram, bank):
    assert bank['archive_sha256'] == ARCHIVE and len(bank['sources']) == 4
    dual = weights(z, raw_gram)
    def apply(values):
        return sum(x*Fraction(y) for x, y in zip(dual, values, strict=True))

    for entry in bank['sources']:
        eta_k = apply(entry['coordinates'])
        total = Fraction(entry['eta0']) + eta_k
        for check in entry['reductions']:
            p = check['prime']
            assert residue(eta_k, p) == check['etaK']
            assert residue(total, p) == check['corrected_sum']
        entry.update(etaK=str(eta_k), corrected_sum=str(total))
    control = read(z, 'exact-score/exact_normalized_class.json')
    assert Fraction(bank['sources'][0]['etaK']) == -Fraction(control['etaK_corrected'])
    coords = read(z, 'cochains/inputs/exact-restriction-and-lift-coordinates.json')['blocks']
    assert bank['sources'][0]['coordinates'] == [str(s) for block in coords for s in block['coordinates']]
    boundary = Boundary(read(z, 'collected-O0/color-basis-native.json'))
    rows = read(z, 'cochains/inputs/exact-coordinate-rows.json')['rows']
    keys = {tuple(k) for row in rows for k, _ in row}
    values = {k: boundary.coefficient(k) for k in keys}
    coordinate = [sum(values[tuple(k)] * Fraction(c) for k, c in row) for row in rows]
    eta_k = apply(coordinate)
    check_y = json.loads((HERE/'endpoint-metric.json').read_text())['pairings']['Q0Y']
    assert Fraction(check_y['corrected_eta0']) + eta_k == 0
    bank.update(etaK_weights=list(map(str, dual)),
                gram11=dict(weight_pairs=573459, nonzero=563792,
                            common_denominator=757350, rotation_and_dimension_scale=-157339),
                Q0Y=dict(raw_coordinate_inputs=len(values), eta0=check_y['corrected_eta0'],
                         etaK=str(eta_k), corrected_sum='0'))
    return bank


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('bank', type=Path)
    parser.add_argument('--gram', type=Path, help='reuse the raw text output of gram11.cpp')
    parser.add_argument('--compiler', default='c++')
    args = parser.parse_args()
    if not __debug__ or sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE:
        parser.error('use the pinned archive and enabled assertions')
    with ZipFile(args.archive) as archive:
        gram = args.gram.read_text() if args.gram else rebuild_gram(archive, args.compiler)
        print(json.dumps(evaluate(archive, gram, json.loads(args.bank.read_text())), indent=2))
