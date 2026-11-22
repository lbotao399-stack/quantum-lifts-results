"""Check the first page and the glueball spectator in thirteen colours."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parent))
import mixed_incoming as traces
from fortuity.full.letters import Q_letter
from orthogonal_glueball_first_pair import GHOSTS, PHI
from orthogonal_glueball_relation_span import relation_dimension
from orthogonal_glueball_second_pair import S, first_factorization, first_score, predecessors
from orthogonal_glueball_sources import physical_image_relations
from rank_four_lowered_replay import Q1d_multi
from q1d_exact import pair_terms
from verify_dot_coordinates import check_thirteen, read_bundle


def collect(terms):
    result = Counter()
    for word, value in terms:
        sign, canonical = traces.so_canon(word)
        if sign:
            result[canonical] += sign * value
    return {word: value for word, value in result.items() if value}


def check_spectator():
    first = PHI + GHOSTS
    second = {(x[0]+y[0], x[1]+y[1], tuple(sorted(x[2]+y[2])))
              for x, y in combinations_with_replacement(first, 2)
              if not set(x[2]).intersection(y[2])}
    assert len(second) == 12 and all(not Q_letter(x) for x in GHOSTS)
    for jet in set(first) | second:
        for ghost in GHOSTS:
            assert not pair_terms(jet, ghost) and not pair_terms(ghost, jet)
    nonzero_inputs = 0
    for i, scalar in enumerate(PHI):
        psi = (0, 0, tuple(j for j in range(3) if j != i))
        for word in (((psi, scalar),), ((psi, GHOSTS[0]), (GHOSTS[1], scalar))):
            image = Q1d_multi(word, 13, 'SO_odd')
            before = collect(Q1d_multi(word+(GHOSTS,), 13, 'SO_odd').items())
            after = collect((w+(GHOSTS,), value) for w, value in image.items())
            assert before == after and before
            nonzero_inputs += 1
    return dict(first_jets=len(first), second_jets=len(second), nonzero_trace_checks=nonzero_inputs)


def main(archives):
    pure = (S,)*8
    assert not traces.q({pure: 1}) and not Q1d_multi(pure, 13, 'SO_odd')
    lower = [relation_dimension(13, 8, k, False) for k in (0, 1)]
    assert lower == [0, 0]
    first = {w: Fraction(c) for w in predecessors([pure]) if (c := first_score(13, w))}
    columns = first_factorization(13, 8, first)
    assert len(first) == 13 and columns == 33
    balanced = relation_dimension(13, 8, 2)
    singlet = physical_image_relations(8, 2, 13)
    assert (balanced, singlet) == (301, 81)
    files = read_bundle(archives, 'so13-j42')
    coefficients = [check_thirteen(files, alt) for alt in (False, True)]
    assert coefficients == [-921600, -921600]
    print(json.dumps(dict(earlier_finite_relations=lower, first_dual_words=len(first),
                          first_dual_columns=columns, second_balanced_relations=balanced,
                          second_spin_singlet_relations=singlet, **check_spectator(),
                          source_coefficients=list(map(str, coefficients)), possible_orders=[2,3])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dot_archive_directory', type=Path)
    if not __debug__:
        parser.error('run with assertions enabled')
    main(parser.parse_args().dot_archive_directory)
