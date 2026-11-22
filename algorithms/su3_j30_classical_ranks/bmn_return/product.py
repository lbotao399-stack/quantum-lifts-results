"""Reconstruct the SU(3) BMN product and apply the eight signed Taylor substitutions."""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ALG = HERE.parents[1]
sys.path[:0] = [str(ALG/'dual_singlets_20261118'), str(ALG/'su3_j34_20261118')]
from dual633_pairing import integral
from gaugino_primary_order import charge
from mixed_incoming import sector
from so7_bmn_633 import export_trace
from so8_level24_descendant import check_letter, substitution
from sp3_omega_product_replay import trace_factor
import su_trace_adapter_v1 as su
from fortuity.full.letters import is_ferm
from fortuity.full.words import canon_multi, canon_word, word_parity


FACTORS = ((0,0,0,0,2,0,0,0), (0,0,1,0,1,0,0,2),
           (0,0,1,1,0,0,0,2), (0,0,2,0,1,0,0,3))
STEPS = ((0,1), (0,2), (1,1), (1,2))


def bmn_product():
    result = {(): 1}
    for description in FACTORS:
        factor = trace_factor(description, unitary=True)
        assert factor and not su.differential(factor)
        expanded = Counter()
        for (trace,), coefficient in factor.items():
            for word, value in result.items():
                sign, key = canon_multi((*word, trace))
                if sign:
                    expanded[key] += value * coefficient * sign
        result = {word: n for word,n in expanded.items() if n}
    return su.projection.project(result, sector((0,0,5,4,4), 9))


@lru_cache(maxsize=4096)
def trace_action(trace, spin, flavour, upward):
    result = Counter()
    phase = 1
    for i,letter in enumerate(trace):
        coefficient, changed = substitution(letter, spin, flavour, upward)
        if coefficient:
            sign, word = canon_word(trace[:i] + (changed,) + trace[i+1:])
            if sign:
                result[word] += phase * coefficient * sign
        if is_ferm(letter):
            phase = -phase
    return tuple((k,v) for k,v in result.items() if v)


def action(row, spin, flavour, upward):
    result = Counter()
    for word,value in row.items():
        phase = 1
        for i,trace in enumerate(word):
            for changed, coefficient in trace_action(trace, spin, flavour, upward):
                sign, key = canon_multi(word[:i] + (changed,) + word[i+1:])
                if sign:
                    result[key] += phase * coefficient * value * sign
            if word_parity(trace):
                phase = -phase
    return {k:v for k,v in result.items() if v}


def digest(row):
    encoded = json.dumps(sorted(row.items()), separators=(',', ':')).encode()
    return sha256(encoded).hexdigest()


def verify(native_path):
    polynomial = bmn_product()
    initial, denominator = integral(polynomial)
    assert len(initial) == 477 and not su.differential(initial)
    assert su.conjugate(initial) == {w:-n for w,n in initial.items()}
    assert charge(initial) == ((0,0,5,4,4),9)
    terms, file_denominator, native_digest = export_trace(polynomial, native_path)
    assert terms == len(initial) and denominator == file_denominator
    current = initial
    sizes = []
    for upward, steps in ((True, STEPS), (False, reversed(STEPS))):
        for spin,flavour in steps:
            letters = {l for w in current for trace in w for l in trace}
            for letter in letters:
                check_letter(letter, spin, flavour, upward)
            current = action(current, spin, flavour, upward)
            sizes.append(len(current))
        if upward:
            assert charge(current) == ((2,2,5,2,2),9)
            for root in ((0,1),(1,2),('spin',1)):
                assert not su.projection.action(current, *root)
            raised_digest = digest(current)
    assert current == {w:240*v for w,v in initial.items()}
    return dict(factors=FACTORS, common_denominator=denominator, steps=STEPS,
                bmn_terms=terms, terms_after_substitutions=sizes, return_coefficient=240,
                native_bmn_sha256=native_digest, raised_numerator_sha256=raised_digest)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export-bmn', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    with TemporaryDirectory() as tmp:
        print(json.dumps(verify(args.export_bmn or Path(tmp)/'bmn.txt'), indent=2))
