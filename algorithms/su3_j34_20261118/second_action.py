"""Apply the existing full-letter SU(3) Q2 to the projected level-34 source."""

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from math import lcm
from pathlib import Path

from verify import add_scaled, exterior_sign, multiply


HERE = Path(__file__).resolve().parent
BMN_LETTERS = [(0, 0, (i,)) for i in range(3)] + [
    (0, 0, (1, 2)), (0, 0, (0, 2)), (0, 0, (0, 1)), (0, 0, (0, 1, 2))]
BMN_SIGNS = (1, 1, 1, 1, -1, 1, -1)
SOURCE_SHA = "673bbcbcf9bafba63b66dfab9a19c306178c5839269abe35da538505166b66c1"


def read_source(path):
    import su_trace_adapter_v1 as algebra
    from fortuity.bmn.words import QLET
    from fortuity.full.letters import Q_letter

    for label, letter in enumerate(BMN_LETTERS):
        left, right = Counter(), Counter()
        for c, first, second in Q_letter(letter):
            left[first, second] += BMN_SIGNS[label] * c
        for c, (first, second) in QLET.get(label, []):
            right[BMN_LETTERS[first], BMN_LETTERS[second]] += c * BMN_SIGNS[first] * BMN_SIGNS[second]
        assert left == right

    raw = gzip.decompress(path.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA
    packet = json.loads(raw)
    source = defaultdict(Fraction)
    for c, traces in packet["terms"]:
        sign = 1
        for trace in traces:
            for label in trace:
                sign *= BMN_SIGNS[label]
        word = tuple(tuple(BMN_LETTERS[label] for label in trace) for trace in traces)
        source[word] += sign * Fraction(c)
    result = algebra.canonical(source)
    assert len(result) == packet["n_terms"] == 6270
    assert algebra.conjugate(result) == result
    assert {algebra.fw.sector_of(w) for w in result} == {(6, 6, 16, 10, 8)}
    return result


def native_image(chunk):
    import su_quantum_integer_v1 as quantum
    return quantum.act(dict(chunk), 2, 3)


def cartan(image):
    """Evaluate all traces on generic traceless diagonal supermatrices over Q."""
    from fortuity.full import words
    from fortuity.full.letters import is_ferm

    def sort_sign(odd):
        if len(set(odd)) != len(odd):
            return 0
        return (-1) ** sum(x > y for i, x in enumerate(odd) for y in odd[i + 1:])

    # An odd colour-diagonal species has just two independent Grassmann coordinates.
    reduced = defaultdict(Fraction)
    for word, coefficient in image.items():
        odd = [letter for trace in word for letter in trace if is_ferm(letter)]
        if any(n > 2 for n in Counter(odd).values()):
            continue
        sign = 1
        for trace in word:
            sign *= sort_sign([letter for letter in trace if is_ferm(letter)])
        if sign:
            exchange, monomial = words.canon_multi(tuple(tuple(sorted(trace)) for trace in word))
            if exchange:
                reduced[monomial] += sign * exchange * coefficient
    reduced = {w: c for w, c in reduced.items() if c}
    alphabet = sorted({letter for word in reduced for trace in word for letter in trace})
    labels = {letter: i for i, letter in enumerate(alphabet)}

    @lru_cache(maxsize=4096)
    def trace_value(trace):
        result = {}
        for colour, indices in enumerate(((0,), (1,), (0, 1))):
            value = {(0, 0): 1}
            for letter in trace:
                factors = {}
                for i in indices:
                    slot = 2 * labels[letter] + i
                    key = (0, 1 << slot) if is_ferm(letter) else (1 << (4 * slot), 0)
                    factors[key] = -1 if colour == 2 else 1
                value = multiply(value, factors)
            add_scaled(result, value, 1)
        return result

    total = {}
    for word, coefficient in reduced.items():
        assert sum(map(len, word)) < 16
        value = {(0, 0): 1}
        for trace in word:
            value = multiply(value, trace_value(trace))
        add_scaled(total, value, coefficient)
    return dict(formal_diagonal_trace_terms=len(reduced), alphabet=alphabet,
                terms=[[str(e), str(o), str(c)] for (e, o), c in sorted(total.items())])


def save_rows(path, rows):
    """Hash the uncompressed, deterministic JSON lines as they are written."""
    digest = hashlib.sha256()
    temporary = path.with_suffix(".partial")
    with temporary.open("wb") as file, gzip.GzipFile(filename="", fileobj=file, mode="wb", mtime=0) as stream:
        for item in rows:
            line = (json.dumps(item, separators=(",", ":")) + "\n").encode()
            digest.update(line)
            stream.write(line)
    temporary.replace(path)
    return digest.hexdigest()


def run(runtime, source_path, output, workers):
    import multiprocessing
    sys.path.insert(0, str(runtime.resolve()))
    import su_trace_adapter_v1 as algebra
    import su_quantum_integer_v1 as quantum

    source = read_source(source_path)
    entries = sorted(source.items())
    # Check nonzero examples from both field-strength contents with the Python trace rule.
    controls = []
    for has_field_strength in (False, True):
        for i, (word, coefficient) in enumerate(entries):
            if any(BMN_LETTERS[6] in trace for trace in word) != has_field_strength:
                continue
            example = {word: coefficient}
            native = quantum.act(example, 2, 3)
            if native:
                assert native == quantum.reference(example, 2, 3)
                controls.append(dict(source_index=i, field_strength=has_field_strength, terms=len(native)))
                break
    assert len(controls) == 2
    # SU(3), tr(A^2) = 2(u^2+uv+v^2). For two odd species, reversing the trace flips the sign.
    scalar = BMN_LETTERS[0]
    diagonal_test = cartan({((scalar, scalar),): 1})
    assert diagonal_test["terms"] == [["2", "0", "2"], ["17", "0", "2"], ["32", "0", "2"]]
    psi, other = BMN_LETTERS[3:5]
    assert cartan({((psi, other),): 1, ((other, psi),): 1})["terms"] == []
    assert cartan({((psi, psi),): 1})["terms"] == []
    assert exterior_sign(1, 2) == 1 and exterior_sign(2, 1) == -1

    chunks = [entries[start:start + 384] for start in range(0, len(entries), 384)]
    image = {}
    if workers == 1:
        for chunk in chunks:
            add_scaled(image, native_image(chunk), 1)
    else:
        with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("fork")) as pool:
            for part in pool.map(native_image, chunks):
                add_scaled(image, part, 1)
    assert image and algebra.conjugate(image) == image
    assert {sum(map(len, w)) for w in image} == {10}
    assert {algebra.fw.sector_of(w) for w in image} == {(5, 5, 16, 10, 8)}
    output.mkdir(parents=True, exist_ok=True)
    digest = save_rows(output / "Q2.jsonl.gz", ([str(image[w]), w] for w in sorted(image)))
    restricted = cartan(image)
    with (output / "Q2_Cartan.json").open("w") as file:
        json.dump(restricted, file, separators=(",", ":"))
        file.write("\n")
    record = dict(source_sha256=SOURCE_SHA, source_terms=len(entries), charge_conjugation=1,
                  group="SU(3)", order=2, coefficient_field="Q", dictionary_letters_checked=7,
                  python_rule_controls=controls, output_length=10, doubled_charges=[5, 5, 16, 10, 8],
                  raw_trace_terms=len(image), common_denominator=lcm(*(c.denominator for c in image.values())),
                  output_jsonl_sha256=digest, cartan_trace_terms=restricted["formal_diagonal_trace_terms"],
                  cartan_monomials=len(restricted["terms"]),
                  cartan_json_sha256=hashlib.sha256((output / "Q2_Cartan.json").read_bytes()).hexdigest(),
                  first_correction_solved=False, second_obstruction_computed=False)
    path = output / "second_action.json"
    path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", type=Path, default=HERE.parent / "su3_independent_prime_20261118"
                        / "runtime/work/sp3-native.nosync")
    parser.add_argument("--source", type=Path, default=HERE / "O34_C_even.json.gz")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    assert args.workers > 0
    run(args.runtime, args.source, args.output, args.workers)
