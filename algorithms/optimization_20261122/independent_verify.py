"""Independent coordinate and sign audit of dihedral_traces.py.

Run with --baseline pointing to the read-only pinned research checkout.
This verifier does not infer finite-rank identities or quantum obstructions.
"""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import factorial, gcd
from pathlib import Path
import subprocess
import sys

import dihedral_traces as fast


def canonical_trace(word):
    """Literal complete group orbit, independent of the fast construction."""
    word = tuple(word)
    reverse = word[::-1]
    forward = min(word[j:] + word[:j] for j in range(len(word)))
    backward = min(reverse[j:] + reverse[:j] for j in range(len(word)))
    if forward == backward and len(word) % 2:
        return forward, 0
    return (forward, 1) if forward <= backward else (backward, (-1) ** len(word))


def identify(word, labels):
    parts, sign = [], 1
    for trace in word:
        canonical, factor = canonical_trace(tuple(labels[j] for j in trace))
        parts.append(canonical)
        sign *= factor
    return tuple(sorted(parts)), sign


def independent_counts(maximum):
    """Coefficient recurrence from (2-2z) A'(z)=(2z-z*z) A(z)."""
    values = [1, 0, 1]
    for n in range(2, maximum):
        values.append(n * values[n] + n * values[n - 1]
                      - (n * (n - 1) // 2) * values[n - 2])
    return values[:maximum + 1]


def digest(value):
    return sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def even_bracelet_burnside(content):
    """Independent fixed-content count for even trace length."""
    degree = sum(content)
    assert degree >= 2 and degree % 2 == 0

    def multinomial(counts):
        value = factorial(sum(counts))
        for count in counts:
            value //= factorial(count)
        return value

    rotations = []
    for shift in range(degree):
        length = degree // gcd(degree, shift)
        rotations.append(multinomial(tuple(c // length for c in content))
                         if all(c % length == 0 for c in content) else 0)
    edge = multinomial(tuple(c // 2 for c in content)) if all(c % 2 == 0 for c in content) else 0
    vertex = 0
    for left, right in product(range(len(content)), repeat=2):
        remaining = tuple(c - int(i == left) - int(i == right)
                          for i, c in enumerate(content))
        if min(remaining) >= 0 and all(c % 2 == 0 for c in remaining):
            vertex += multinomial(tuple(c // 2 for c in remaining))
    count, remainder = divmod(sum(rotations) + degree // 2 * (edge + vertex), 2 * degree)
    assert remainder == 0
    return dict(content=content, raw_words=multinomial(content),
                rotation_fixed_counts=rotations, edge_reflection_fixed=edge,
                vertex_reflection_fixed=vertex, bracelets=count)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--native", action="store_true", help="also compile and audit the native catalogue")
    args = parser.parse_args()
    if not __debug__:
        parser.error("assertions must be enabled")
    source = args.baseline / "algorithms/dual_singlets_20261118/f4/j36_source"
    sys.path.insert(0, str(source))
    import mixed_colour_relations as old

    counts = independent_counts(14)
    assert counts == [fast.trace_count(n) for n in range(15)]
    records = []
    lists = {}
    for n in range(9):
        old.colour_trace.cache_clear()
        expected = old.trace_words(n)
        actual = fast.trace_words(n)
        assert actual == expected
        assert len(actual) == counts[n] == len(set(actual))
        for multi in actual:
            assert sorted(j for word in multi for j in word) == list(range(n))
            assert multi == tuple(sorted(multi))
            for word in multi:
                assert len(word) >= 2
                assert canonical_trace(word) == (word, 1)
        lists[n] = actual
        records.append(dict(degree=n, count=len(actual), sha256=digest(actual)))

    action_records = []
    for n in range(2, 9):
        words = lists[n]
        lookup = {word: j for j, word in enumerate(words)}
        actions = []
        for i, (old_indices, old_signs) in enumerate(old.signed_transpositions(words)):
            labels = list(range(n))
            labels[i], labels[i + 1] = labels[i + 1], labels[i]
            values = [identify(word, labels) for word in words]
            indices = [lookup[word] for word, sign in values]
            signs = [sign for word, sign in values]
            assert indices == old_indices and signs == old_signs.tolist()
            assert all(indices[indices[j]] == j
                       and signs[j] * signs[indices[j]] == 1 for j in range(len(words)))
            actions.append((indices, signs))
        # Relations of the full signed symmetric-group action.
        def compose(j, sequence):
            sign = 1
            for i in sequence:
                indices, signs = actions[i]
                sign *= signs[j]
                j = indices[j]
            return j, sign
        for i in range(n - 2):
            assert all(compose(j, (i, i + 1, i)) == compose(j, (i + 1, i, i + 1))
                       for j in range(len(words)))
        for i in range(n - 1):
            for k in range(i + 2, n - 1):
                assert all(compose(j, (i, k)) == compose(j, (k, i))
                           for j in range(len(words)))
        action_records.append(dict(degree=n, coordinates=len(words),
                                   adjacent_actions=len(actions), coxeter_relations=True))

    content_cases = set()
    for n in range(9):
        for a in range(n + 1):
            content_cases.add((a, n - a))
    content_cases.update([(5, 3), (4, 4), (2, 2, 2, 2), (3, 2, 2), (0, 3, 0, 4),
                          (), (0, 0), (1, 1, 1), (1, 1, 2, 2)])
    content_records = []
    for counts_case in sorted(content_cases):
        n = sum(counts_case)
        labels = tuple(i for i, count in enumerate(counts_case) for _ in range(count))
        expected = sorted({word for multi in lists[n]
                           for word, sign in [identify(multi, labels)] if sign})
        actual = fast.trace_words_fixed_content(counts_case)
        assert actual == expected, counts_case
        for maximum in (0, 2, 3, 4, n):
            restricted = [multi for multi in expected
                          if all(len(word) <= maximum for word in multi)]
            assert fast.trace_words_fixed_content(counts_case, maximum) == restricted
        content_records.append(dict(content=counts_case, count=len(actual),
                                    sha256=digest(actual)))

    # Literal brute-force necklace orbits, including periodic and zero-sign words.
    necklace_cases = []
    for content in ((0, 0), (1,), (2,), (3,), (4,), (2, 2), (3, 3),
                    (0, 2, 0, 2), (1, 1, 1), (2, 2, 2), (3, 1, 2)):
        n = sum(content)
        raw = []
        if n >= 2:
            for word in product(range(len(content)), repeat=n):
                if tuple(word.count(i) for i in range(len(content))) == content:
                    raw.append(word)
        cyclic = sorted({min(word[j:] + word[:j] for j in range(n)) for word in raw})
        assert sorted(fast.necklaces_fixed_content(content)) == cyclic
        bracelets = sorted({canonical for word in raw
                            for canonical, sign in [canonical_trace(word)] if sign})
        assert sorted(fast._bracelets(content)) == bracelets
        necklace_cases.append(dict(content=content, necklaces=len(cyclic),
                                   nonzero_bracelets=len(bracelets)))

    native_result = None
    if args.native:
        import numpy as np
        import native_catalogue

        contents = [(2,), (3,), (4,), (2, 2), (3, 3), (0, 2, 0, 2),
                    (1, 1, 1), (2, 2, 2), (3, 1, 2)]
        for content in contents:
            values = native_catalogue.fixed_content_catalogue(content)
            decoded = [tuple((int(code) >> (2 * j)) & 3
                             for j in range(sum(content) - 1, -1, -1)) for code in values]
            assert decoded == sorted(fast._bracelets(content))
        content = (4, 4, 4, 2)
        burnside = even_bracelet_burnside(content)
        values = native_catalogue.catalogue()
        assert len(values) == burnside["bracelets"] == 112950
        assert values.dtype == np.int32 and np.all(values[1:] > values[:-1])
        for code in values:
            word = tuple((int(code) >> (2 * j)) & 3 for j in range(13, -1, -1))
            assert tuple(word.count(i) for i in range(4)) == content
            assert canonical_trace(word) == (word, 1)
        # Test the native canonicalizer on the complete dihedral orbit of
        # every independently certified representative, including stabilizers.
        forward = values.astype(np.uint32)
        reverse = np.zeros_like(forward)
        for j in range(14):
            reverse = (reverse << 2) | ((forward >> (2 * j)) & 3)
        for shift in range(14):
            assert np.array_equal(native_catalogue.canonical14(forward), values)
            assert np.array_equal(native_catalogue.canonical14(reverse), values)
            forward = ((forward << 2) & ((1 << 28) - 1)) | (forward >> 26)
            reverse = ((reverse << 2) & ((1 << 28) - 1)) | (reverse >> 26)
        native_result = dict(burnside=burnside, small_content_cases=len(contents),
                             complete_orbit_checks=28 * len(values),
                             catalogue_sha256=sha256(values.astype("<i4").tobytes()).hexdigest())

    result = dict(status="PASS", baseline_commit=subprocess.check_output(
                  ["git", "-C", str(args.baseline), "rev-parse", "HEAD"], text=True).strip(),
                  baseline_sha256={name: sha256((source / name).read_bytes()).hexdigest()
                                   for name in ("mixed_colour_relations.py", "dressed_tree.py")},
                  independent_count_recurrence=counts,
                  baseline_coordinate_lists=records, signed_actions=action_records,
                  repeated_content=content_records, brute_orbits=necklace_cases,
                  native_catalogue=native_result,
                  scope="Exact coordinate/sign equivalence; no quantum-lift or completeness claim.")
    serialized = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(json.dumps(dict(status="PASS", labelled_degrees=len(records),
                          signed_action_degrees=len(action_records),
                          repeated_content_cases=len(content_records),
                          brute_orbit_cases=len(necklace_cases))))


if __name__ == "__main__":
    main()
