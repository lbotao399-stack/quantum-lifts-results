"""Complete the SO(11) third-order chain and query its stable Pfaffian derivative."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import chain
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import mixed_incoming as traces
from fortuity.full.wordenum import multitraces_with_content
from orthogonal_glueball_first_pair import GHOSTS, PHI, PfaffianCoefficient
from orthogonal_glueball_second_pair import tree_split
from orthogonal_odd_trace_primitive import primitive
from rank_four_lowered_replay import Q1d_multi


LETTERS = dict(zip('ABC+-', PHI + GHOSTS))


def leading_derivative(word, relation):
    """Coefficient of an ordered trace in the canonical stable Q0 O0 expansion.

    An input list gives the traces in their displayed order.  Return the
    coefficient multiplied by the sign relating that list to its canonical form.
    """
    if Counter(''.join(word)) != Counter('AAABBBCCC+++---') or any(len(t) < 2 for t in word):
        return Fraction()
    ordered = tuple(tuple(LETTERS[c] for c in t) for t in word)
    sign, canonical = traces.so_canon(ordered)
    if not sign:
        return Fraction()
    # Q0 is a derivation. Reversal of the first odd trace is a right inverse
    # on each such word, since its graded reverse represents its negative.
    assert traces.q(primitive(canonical)) == {canonical: 1}
    return sign * relation.coefficient(canonical)


def check_terminal(export=False):
    # Enumerate cyclic, graded multitraces first.  Orthogonal reversal is then
    # applied to every trace, before testing the classical equation.
    cyclic = multitraces_with_content(tuple(sorted(GHOSTS * 6)))
    basis = {w for s, w in map(traces.so_canon, cyclic) if s}
    sign, target = traces.so_canon((GHOSTS,) * 6)
    assert sign == 1
    primitive_support = set()
    sections = []
    columns = []
    symbols = {v: k for k, v in LETTERS.items()}
    symbols.update({(0, 2, ()): 'u', (1, 1, ()): 'v', (2, 0, ()): 'w'})
    def encode(word):
        return [''.join(symbols[letter] for letter in trace) for trace in word]
    for word in sorted(basis):
        cochain, section = tree_split(word)
        expected = Counter({word: 1})
        expected.subtract(section)
        assert traces.q(cochain) == {w: v for w, v in expected.items() if v}
        primitive_support.update(cochain)
        if section:
            assert section == {target: 1} and word == target
            sections.append(word)
        if export:
            columns.append(dict(word=encode(word), section=int(bool(section)),
                                h=[[str(v), encode(w)] for w, v in sorted(cochain.items())]))
    assert sections == [target]
    for word in primitive_support:
        letters = list(chain.from_iterable(word))
        assert len(letters) == 11 and all(not theta for _, _, theta in letters)
        assert sum(a for a, _, _ in letters) == sum(b for _, b, _ in letters) == 6
        assert not Q1d_multi(word, 11, 'SO_odd')
    summary = dict(cyclic_words=len(cyclic), orthogonal_words=len(basis),
                   primitive_words=len(primitive_support), surviving_sections=len(sections))
    if export:
        summary['columns'] = columns
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--derivative', nargs='+', metavar='TRACE',
                        help='return the signed coefficient on a list of traces using A,B,C,+,-')
    parser.add_argument('--terminal-map', action='store_true',
                        help='print all columns of the final primitive in divided jet coordinates')
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    if args.derivative:
        if any(c not in LETTERS for c in ''.join(args.derivative)):
            parser.error('use first-jet labels A, B, C, +, -')
        print(leading_derivative(args.derivative, PfaffianCoefficient(3, 0)))
    else:
        result = check_terminal(args.terminal_map)
        pf = PfaffianCoefficient(3, 0)
        assert leading_derivative('ABC A+ A- B+ B- C+ C-'.split(), pf) == 540
        assert leading_derivative('ABC ABC ABC +- +- +-'.split(), pf) == 360
        print(json.dumps(result))
