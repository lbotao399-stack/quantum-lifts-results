"""Enumerate the unchanged multilinear orthogonal trace coordinates exactly.

Replacement for mixed_colour_relations.trace_words.  All labels are distinct,
all colour letters are even, single-letter traces are absent, and reversing a
length-k trace carries the existing sign (-1)**k.  This module only constructs
the canonical *basis labels*.  It never applies finite-colour identities or a
projection.  Signed substitutions must still use the caller's colour_trace.
"""

from itertools import combinations, permutations
from math import comb, factorial


def trace_count(degree):
    """Number of canonical labelled multi-traces, using an exact recurrence."""
    if not isinstance(degree, int) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    counts = [1, 0]
    for n in range(2, degree + 1):
        counts.append(sum(comb(n - 1, k - 1)
                          * (1 if k == 2 else factorial(k - 1) // 2)
                          * counts[n - k] for k in range(2, n + 1)))
    return counts[degree]


def iter_trace_words(labels):
    """Generate every canonical multi-trace on distinct sorted labels once.

    The least remaining label starts the first trace.  Each trace has that
    smallest label first, and for length at least three its second label is
    smaller than its last.  This fixes cyclic rotation and reversal uniquely.
    Recursive use of the least remaining label fixes the ordering of traces.
    """
    labels = tuple(labels)
    if labels != tuple(sorted(set(labels))):
        raise ValueError("labels must be sorted and distinct")
    if not labels:
        yield ()
        return
    if len(labels) == 1:
        return
    first, tail = labels[0], labels[1:]
    for length in range(2, len(labels) + 1):
        if len(labels) - length == 1:
            continue
        for members in combinations(tail, length - 1):
            members_set = frozenset(members)
            remaining = tuple(x for x in tail if x not in members_set)
            if length == 2:
                for rest in iter_trace_words(remaining):
                    yield ((first,) + members,) + rest
            else:
                # The endpoints are unordered before orientation; selecting
                # them as an increasing pair avoids generating both reversals.
                for left, right in combinations(members, 2):
                    middle = tuple(x for x in members if x not in (left, right))
                    for order in permutations(middle):
                        trace = (first, left) + order + (right,)
                        for rest in iter_trace_words(remaining):
                            yield (trace,) + rest


def trace_words(degree):
    """Return exactly the old sorted list; suitable for unchanged column IDs."""
    if not isinstance(degree, int) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    return sorted(iter_trace_words(tuple(range(degree))))


def necklaces_fixed_content(counts):
    """Canonical cyclic words on an even alphabet with the given multiplicity.

    FKM recursion generates one necklace per cyclic orbit; it does not build
    all labelled words or permute equal labels.  Internal support for the
    repeated-content front end, not a claim about arbitrary superletters.
    """
    counts = tuple(counts)
    n = sum(counts)
    if n < 2 or not counts or min(counts) < 0:
        return
    active = tuple(i for i, c in enumerate(counts) if c)
    remaining = [counts[i] for i in active]
    word = [0] * (n + 1)
    remaining[0] -= 1

    def visit(position, period):
        if position > n:
            if n % period == 0:
                yield tuple(active[i] for i in word[1:])
            return
        lower = word[position - period]
        for letter in range(lower, len(active)):
            if remaining[letter]:
                word[position] = letter
                remaining[letter] -= 1
                yield from visit(position + 1, period if letter == lower else position)
                remaining[letter] += 1

    yield from visit(2, 1)


def _bracelets(counts):
    """Nonvanishing trace representatives with the old transpose convention."""
    for word in necklaces_fixed_content(counts):
        reverse = word[::-1]
        backwards = min(reverse[i:] + reverse[:i] for i in range(len(word)))
        if word > backwards or (word == backwards and len(word) % 2):
            continue
        yield word


def trace_words_fixed_content(counts, maximum_trace_length=None):
    """All formal nonzero multi-traces for ordinary even Lie-algebra letters.

    Returns exactly the list obtained by identifying distinct labelled traces
    with repeated labels and discarding zero reversal orbits.  Coefficients
    and signed relation specialization remain the responsibility of the caller.
    """
    from itertools import product

    counts = tuple(counts)
    if any(not isinstance(c, int) or c < 0 for c in counts):
        raise ValueError("counts must be nonnegative integers")
    n = sum(counts)
    maximum = n if maximum_trace_length is None else maximum_trace_length
    if not isinstance(maximum, int) or maximum < 0:
        raise ValueError("maximum trace length must be a nonnegative integer")
    generators = sorted((word, content)
                        for content in product(*(range(c + 1) for c in counts))
                        if 2 <= sum(content) <= maximum
                        for word in _bracelets(content))

    def visit(remaining, start):
        if not any(remaining):
            yield ()
            return
        for i in range(start, len(generators)):
            word, content = generators[i]
            if all(a <= b for a, b in zip(content, remaining)):
                leftover = tuple(b - a for a, b in zip(content, remaining))
                if sum(leftover) == 1:
                    continue
                for rest in visit(leftover, i):
                    yield (word,) + rest

    return list(visit(counts, 0))
