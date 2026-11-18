"""Replay the two scalar restrictions and one mixed restriction of the SU(3) level-34 operator."""

import argparse
import gzip
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE_HASH = "ce7fdd0bb3cd6a22b9d294f143553f9f6a71fdde21fba2aa0eb0a4d0bdcc462b"
WORDS_HASH = "e2ecc349c6485c0d2063d60ba955d820d7f5ba01032a84645caba66045de5593"
VENDORED_WORDS_HASH = "f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e"
PUBLIC_PRIME = 67108859


def load_words(path):
    source = path.read_bytes()
    if hashlib.sha256(source).hexdigest() == VENDORED_WORDS_HASH:
        # The copy used by the published solver has one extra trailing newline.
        assert source.endswith(b"\n\n")
        source = source[:-1]
    assert hashlib.sha256(source).hexdigest() == WORDS_HASH
    spec = importlib.util.spec_from_file_location("bmn_words", path)
    words = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(words)
    return words


def add_scaled(total, row, scale):
    for key, value in row.items():
        new_value = total.get(key, 0) + scale * value
        if new_value:
            total[key] = new_value
        else:
            total.pop(key, None)


def load_source(path, words):
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == SOURCE_HASH
    rows = json.loads(gzip.decompress(raw))
    source = {}
    for traces, coefficient in rows:
        # Field strength is 7 in the literal paper packet and 6 in bmn.words.
        word = tuple(tuple(6 if letter == 7 else letter for letter in trace) for trace in traces)
        assert words.charges_of_content(words.content(word)) == (6, 8, 5, 4)
        assert word not in source
        source[word] = Fraction(coefficient)
    assert len(source) == 5815
    return source


def restrict_block(expr, counts, words, fermion=None):
    """Extract a block coefficient, with the single abelian fermion ordered on the left."""
    @lru_cache(None)
    def options(trace, removed_fermion):
        positions = [i for i, letter in enumerate(trace) if letter < 3 and counts[letter]]
        out = Counter()
        for size in range(min(sum(counts), len(positions)) + 1):
            for chosen in combinations(positions, size):
                degree = tuple(sum(trace[i] == label for i in chosen) for label in range(3))
                if any(n > limit for n, limit in zip(degree, counts)):
                    continue
                remaining = tuple(letter for i, letter in enumerate(trace)
                                  if i not in chosen and i != removed_fermion)
                if not remaining:
                    out[degree, ()] += 2 + (-2) ** len(trace)
                else:
                    sign, canonical = words.canon_word(remaining)
                    if sign:
                        out[degree, canonical] += sign
        return [(degree, word, value) for (degree, word), value in out.items() if value]

    result = Counter()
    for multi, coefficient in expr.items():
        placements = [(None, None, 1)]
        if fermion is not None:
            placements = []
            odd_before = 0
            for t, trace in enumerate(multi):
                for i, letter in enumerate(trace):
                    if letter == fermion:
                        placements.append((t, i, (-1) ** odd_before))
                    odd_before += words.is_ferm(letter)
        for selected_trace, position, outer_sign in placements:
            factors = [options(trace, position if t == selected_trace else None)
                       for t, trace in enumerate(multi)]
            for choices in product(*factors):
                if tuple(map(sum, zip(*(choice[0] for choice in choices)))) != counts:
                    continue
                sign, canonical = words.canon_multi([word for _, word, _ in choices if word])
                for _, _, value in choices:
                    sign *= value
                if sign:
                    result[canonical] += coefficient * sign * outer_sign
    return {word: value for word, value in result.items() if value}


def even_part(expr, words):
    result = Counter()
    for multi, coefficient in expr.items():
        for reverse in (False, True):
            sign, traces = 1, []
            for trace in multi:
                if reverse:
                    odd = sum(map(words.is_ferm, trace))
                    sign *= (-1) ** (len(trace) + odd * (odd - 1) // 2)
                cyclic, word = words.canon_word(trace[::-1] if reverse else trace)
                sign *= cyclic
                traces.append(word)
            assert sign
            exchange, word = words.canon_multi(traces)
            result[word] += coefficient * sign * exchange / 2
    return {word: value for word, value in result.items() if value}


@lru_cache(None)
def exterior_sign(left, right):
    if left & right:
        return 0
    inversions = 0
    while right:
        bit = right & -right
        inversions += (left >> bit.bit_length()).bit_count()
        right ^= bit
    return (-1) ** inversions


def multiply(left, right):
    result = Counter()
    for (even_left, odd_left), x in left.items():
        for (even_right, odd_right), y in right.items():
            sign = exterior_sign(odd_left, odd_right)
            if sign:
                result[even_left + even_right, odd_left | odd_right] += x * y * sign
    return {key: value for key, value in result.items() if value}


@lru_cache(None)
def trace_polynomial(trace):
    """Generic traceless 2 by 2 matrices [[u,v],[w,-u]] in a polynomial/exterior ring."""
    assert len(trace) < 16
    result = Counter()
    for path in product((0, 1), repeat=len(trace)):
        even = odd = 0
        sign = 1
        for label, i, j in zip(trace, path, path[1:] + path[:1]):
            entry = 0 if i == j else 1 + i
            if i == j == 1:
                sign = -sign
            if 3 <= label <= 5:
                bit = 1 << (3 * (label - 3) + entry)
                sign *= exterior_sign(odd, bit)
                if not sign:
                    break
                odd |= bit
            else:
                species = 3 if label == 6 else label
                # Each four-bit digit holds the exponent of one even matrix entry.
                even += 1 << (4 * (3 * species + entry))
        if sign:
            result[even, odd] += sign
    return {key: value for key, value in result.items() if value}


@lru_cache(maxsize=4096)
def polynomial_of_word(multi):
    assert sum(map(len, multi)) < 16
    result = {(0, 0): 1}
    for trace in multi:
        result = multiply(result, trace_polynomial(trace))
    return result


def polynomial(expr):
    result = {}
    for multi, coefficient in expr.items():
        add_scaled(result, polynomial_of_word(multi), coefficient)
    return result


def reduce_row(row, basis):
    residual, used = row.copy(), {}
    for pivot, vector, combination in basis:
        coefficient = residual.get(pivot)
        if coefficient:
            add_scaled(residual, vector, -coefficient)
            add_scaled(used, combination, coefficient)
    return residual, used


def boundary_basis(predecessors, words):
    basis = []
    for index, word in enumerate(predecessors):
        residual, used = reduce_row(polynomial(words.Q_multi(word)), basis)
        if residual:
            pivot = min(residual)
            scale = Fraction(1, residual[pivot])
            coefficients = {index: 1}
            add_scaled(coefficients, used, -1)
            basis.append((pivot, {key: scale * v for key, v in residual.items()},
                          {key: scale * v for key, v in coefficients.items()}))
    return basis


def pairing(functional, row):
    return sum(value * row.get(key, 0) for key, value in functional.items())


def residue_functional(basis, coordinate):
    functional = {coordinate: Fraction(1)}
    for pivot, vector, _ in reversed(basis):
        add_scaled(functional, {pivot: 1}, -pairing(functional, vector))
    return functional


def residue_mod_prime(rational):
    rational = Fraction(rational)
    return rational.numerator * pow(rational.denominator, -1, PUBLIC_PRIME) % PUBLIC_PRIME


def test_channel(source, counts, words, fermion=None):
    target = tuple(n - m - int(fermion is not None and i != fermion - 3)
                   for i, (n, m) in enumerate(zip((8, 5, 4), counts)))
    two_j = 6 - int(fermion is not None)
    image = restrict_block(source, counts, words, fermion)
    for word in image:
        assert words.charges_of_content(words.content(word)) == (two_j,) + target
    values = polynomial(image)
    assert not polynomial(words.Q_expr(image))
    predecessors = words.words_of_sector((two_j + 1,) + target, 3)
    basis = boundary_basis(predecessors, words)
    residual, used = reduce_row(values, basis)
    primitive = {predecessors[index]: value for index, value in used.items()}
    replay = values.copy()
    add_scaled(replay, polynomial(words.Q_expr(primitive)), -1)
    assert replay == residual
    witness = None
    if residual:
        coordinate = min(residual)
        functional = residue_functional(basis, coordinate)
        answer = pairing(functional, values)
        assert answer == residual[coordinate] and residue_mod_prime(answer)
        modular = {key: residue_mod_prime(value) for key, value in functional.items()}
        for word in predecessors:
            boundary = polynomial(words.Q_multi(word))
            assert pairing(functional, boundary) == 0
            assert pairing(modular, boundary) % PUBLIC_PRIME == 0
        assert sum(modular.get(key, 0) * residue_mod_prime(value)
                   for key, value in values.items()) % PUBLIC_PRIME == residue_mod_prime(answer)
        witness = dict(monomial_coefficients=[[e, o, str(v)] for (e, o), v in sorted(functional.items())],
                       boundary_pairings_zero=len(predecessors), rational_pairing=str(answer),
                       public_prime_pairing=residue_mod_prime(answer))
    return dict(scalar_powers=counts, abelian_fermion_label=fermion, target=(two_j,) + target, restricted_trace_terms=len(image),
                matrix_entry_monomials=len(values), closed=True, predecessor_trace_words=len(predecessors),
                rational_boundary_rank=len(basis), residual_monomials=len(residual),
                decomposition_primitive=[[multi, str(value)] for multi, value in sorted(primitive.items())],
                witness=witness)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, nargs="?", default=HERE / "published_O34_input.json.gz")
    parser.add_argument("--words", type=Path, default=HERE.parent / "su3_independent_prime_20261118"
                        / "runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/bmn/words.py")
    parser.add_argument("--emit-source", type=Path, help="Save the rational C-even source in BMN letter order")
    args = parser.parse_args()
    words = load_words(args.words)
    source = load_source(args.source, words)
    assert exterior_sign(1, 2) == 1 and exterior_sign(2, 1) == -1 and exterior_sign(1, 1) == 0
    # The lower colour block contributes to the coefficient in a pure scalar trace.
    assert restrict_block({((0, 0),): 1}, (2, 0, 0), words) == {(): 6}
    assert restrict_block({((0, 0, 0),): 1}, (3, 0, 0), words) == {(): -6}
    # Odd abelian extraction is an odd derivative on products of colour traces.
    mixed = {((0, 3), (0, 5)): 1}
    assert restrict_block(mixed, (1, 0, 0), words, 5) == {((0, 3),): -6}
    assert restrict_block(mixed, (1, 0, 0), words, 3) == {((0, 5),): 6}
    projected = even_part(source, words)
    assert even_part(projected, words) == projected
    assert even_part(words.Q_expr(source), words) == words.Q_expr(projected)
    raw_counts = [len(restrict_block(source, counts, words)) for counts in ((3, 1, 0), (4, 0, 0))]
    assert raw_counts == [0, 0]
    packet = dict(source="Kim, Lee, Lee and Song, arXiv:2610.20340, SU(3) O34 highest component",
                  projection="(1+C)/2", field_convention="paper classical fields",
                  letters=["phi1", "phi2", "phi3", "psi1", "psi2", "psi3", "f"],
                  n_terms=len(projected), terms=[[str(v), w] for w, v in sorted(projected.items())])
    packed = json.dumps(packet, separators=(",", ":")).encode()
    assert gzip.decompress((HERE / "O34_C_even.json.gz").read_bytes()) == packed
    if args.emit_source is not None:
        args.emit_source.write_bytes(gzip.compress(packed, mtime=0))
    records = [test_channel(projected, counts, words) for counts in ((3, 1, 0), (4, 0, 0))]
    records.append(test_channel(projected, (3, 0, 0), words, 5))
    assert [r["residual_monomials"] == 0 for r in records] == [True, True, False]
    record = dict(source_sha256=SOURCE_HASH, words_sha256=WORDS_HASH,
                          source_terms=len(source), C_even_terms=len(projected),
                          C_even_json_sha256=hashlib.sha256(packed).hexdigest(),
                          Q0_commutes_with_C_projection_on_source=True,
                          unprojected_scalar_channel_trace_terms=raw_counts,
                          coefficient_field="Q", public_prime=PUBLIC_PRIME,
                          channels=records)
    assert json.loads(json.dumps(record)) == json.loads((HERE / "certificate.json").read_text())
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
