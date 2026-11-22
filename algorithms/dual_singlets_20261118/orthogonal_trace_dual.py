"""Pull a functional on orthogonal first-jet traces through two supercharges."""
from collections import Counter
from fractions import Fraction

import mixed_incoming as traces
from orthogonal_glueball_adjoint_replay import parents
from orthogonal_glueball_second_pair import predecessors, quantum_columns, tree_columns, tree_split
from rank_four_lowered_replay import Q1d_multi


def pullback(colours, functional, sign):
    """Return a row on classical images and verify its transposed equation."""
    quantum = {}
    sources = set(parents(functional, active_even=False)) | quantum_columns(functional)
    for source in sources:
        value = Fraction()
        for word, coeff in Q1d_multi(source, colours, 'SO_odd').items():
            orientation, key = traces.so_canon(word)
            value += orientation * coeff * functional.get(key, 0)
        if value:
            quantum[source] = sign * value
    row = {}
    candidates = set()
    for raw in predecessors(functional):
        orientation, word = traces.so_canon(raw)
        if orientation:
            candidates.add(word)
    for word in candidates:
        primitive, section = tree_split(word)
        expected = Counter({word: 1})
        expected.subtract(section)
        assert traces.q(primitive) == {w:c for w,c in expected.items() if c}
        value = sum(c * quantum.get(w, 0) for w, c in primitive.items())
        if value:
            row[word] = value
    checked = sources | tree_columns(row)
    for source in checked:
        value = sum(c*row.get(w, 0) for w,c in traces.q({source: 1}).items())
        if value != quantum.get(source, 0):
            raise ArithmeticError(('stable adjoint equation failed',source,value,quantum.get(source,0)))
    return row, dict(ordered_candidates=len(candidates), nonzero_words=len(row),
                     quantum_sources=len(sources), checked_columns=len(checked))
