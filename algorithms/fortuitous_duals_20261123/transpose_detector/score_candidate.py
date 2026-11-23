#!/usr/bin/env python3
"""Apply the certified source dual using only Python's standard library.

Input format: the same {"weight": [...], "length": 9, "terms": [[coefficient,
ordered_native_traces], ...]} as positive_sector/omega_rational.json.gz.
This routine scores coefficients. Closure and highest-weight membership are
independent prerequisites for interpreting the score as a fortuitous class.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def decode(word):
    return tuple(tuple((i, j, tuple(theta)) for i, j, theta in trace)
                 for trace in word)


def load_json(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == '.gz' else data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate', type=Path, nargs='?',
                        default=HERE.parent/'positive_sector/omega_rational.json.gz')
    args = parser.parse_args()
    detector_path = HERE/'omega_source_dual.json'
    detector = load_json(detector_path)
    receipt = load_json(HERE/'transpose_result.json')
    assert hashlib.sha256(detector_path.read_bytes()).hexdigest() == receipt['source_dual_sha256']
    candidate = load_json(args.candidate)
    if candidate['weight'] != detector['weight'] or candidate['length'] != detector['length']:
        raise ValueError('The candidate must use the certified weight and length.')
    source = {decode(word): Fraction(coefficient)
              for word, coefficient in detector['source_terms']}
    polynomial = defaultdict(Fraction)
    for coefficient, word in candidate['terms']:
        polynomial[decode(word)] += Fraction(coefficient)
    score = sum(value * polynomial.get(word, 0) for word, value in source.items())
    result = {
        'exact_coefficient_pairing': str(score),
        'source_coordinates': len(source),
        'input_terms': len(polynomial),
        'weight': detector['weight'],
        'length': detector['length'],
        'interpretation': (
            'For a finite Q0-closed cochain in the certified highest-weight block, '
            'a nonzero score proves a nonzero fortuitous component. In this '
            'one-dimensional cohomology block, zero score means the class is zero. '
            'This command does not check Q0-closure or highest-weight membership.'
        ),
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
