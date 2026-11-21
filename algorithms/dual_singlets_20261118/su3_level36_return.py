"""Lower and raise the two even SU(3) BMN classes at level 36 in exact trace coordinates."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
from pathlib import Path

from su3_bmn_descendants import (LEVEL34_SHA, WORDS, add_scaled, flavour_root, highest, load_words,
                                 pairing, polynomial, read_operator, restrict_block)


HERE = Path(__file__).resolve().parent
DATA = HERE/"taylor_hooks"
LEVEL32_SHA = "897fbf27bced5562684f9537ed1bd9e1b8a352f7149cc3b99ddb22c89b6c823e"


def lower_letter(letter, axis, words):
    if letter == words.F:
        return ((2, (words.PSI[axis],)),)
    if letter in words.PSI:
        flavour = letter-3
        return tuple((sign, (k,)) for k in range(3) if (sign := words._eps(axis, k, flavour)))
    return ()


def check_letters(words):
    def on_word(word, rules):
        result = Counter()
        sign = 1
        for i, letter in enumerate(word):
            for n, replacements in rules(letter):
                result[word[:i]+replacements+word[i+1:]] += sign*n
            if words.is_ferm(letter):
                sign = -sign
        return {w: n for w, n in result.items() if n}

    def compose(letter, left, right):
        result = Counter()
        for w, n in on_word((letter,), right).items():
            for v, k in on_word(w, left).items():
                result[v] += n*k
        return result

    for axis in range(3):
        def down(letter):
            return lower_letter(letter, axis, words)

        def up(letter):
            if letter in words.PHI:
                return tuple((c, (3+k,)) for k in range(3) if (c := words._eps(axis, letter, k)))
            return ((1, (words.F,)),) if letter == words.PSI[axis] else ()

        def classical(letter):
            return words.QLET.get(letter, ())

        for letter in range(7):
            mixed = compose(letter, down, classical)
            mixed.update(compose(letter, classical, down))
            assert not any(mixed.values())
            charge = compose(letter, down, up)
            charge.update(compose(letter, up, down))
            flavours = ((letter,) if letter < 3 else
                        tuple(i for i in range(3) if i != letter-3) if letter < 6 else (0, 1, 2))
            assert +charge == Counter({(letter,): sum(i != axis for i in flavours)}) and min(charge.values(), default=0) >= 0
    return 21


def lower(expr, axis, words):
    @lru_cache(maxsize=None)
    def on_trace(trace):
        result = Counter()
        sign = 1
        for i, letter in enumerate(trace):
            for n, replacements in lower_letter(letter, axis, words):
                phase, canonical = words.canon_word(trace[:i]+replacements+trace[i+1:])
                if phase:
                    result[canonical] += sign*n*phase
            if words.is_ferm(letter):
                sign = -sign
        return tuple((w, n) for w, n in result.items() if n)

    result = Counter()
    for multi, coefficient in expr.items():
        sign = 1
        for k, trace in enumerate(multi):
            for w, n in on_trace(trace):
                words.add_term(result, multi[:k]+(w,)+multi[k+1:], coefficient*sign*n)
            if words.word_parity(trace):
                sign = -sign
    return dict(result)


def even_descendant(level, source, words):
    if level == 34:
        return highest(words.Qplus_expr(source, 2), words)
    result = {}
    for j, k, l in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        component = source if j == 0 else flavour_root(source, j, 0, words)
        add_scaled(result, words.Qplus_expr(words.Qplus_expr(component, l), k),
                   Fraction(1) if j == 0 else Fraction(1, 4))
    return result


def read_functional(path, field):
    rows = json.loads(path.read_bytes())[field]["functional"]
    functional = {(e, o): Fraction(c) for e, o, c in rows}
    assert len(rows) == len(functional) and all(functional.values())
    return functional


def calculate(args):
    words = load_words(args.words)
    checked = check_letters(words)
    source = read_operator(args.source, LEVEL32_SHA if args.level == 32 else LEVEL34_SHA)
    original = even_descendant(args.level, source, words)
    assert len(original) == {32: 82088, 34: 48451}[args.level]
    lowered = lower(original, 2, words)
    assert words.expr_sector(lowered) == (6, 8, 5, 4)
    assert not flavour_root(lowered, 0, 1, words) and not flavour_root(lowered, 1, 2, words)
    returned = highest(words.Qplus_expr(lowered, 2), words)
    assert words.expr_sector(returned) == (7, 8, 5, 5)
    pure = read_functional(args.pure, "quartic_descendant")
    mixed = read_functional(args.mixed, "C_even_length11_raised")

    def coordinates(expr):
        result = []
        for counts, fermion, functional in (((4, 0, 0), None, pure), ((3, 0, 0), 5, mixed)):
            active = restrict_block(expr, counts, words, fermion)
            assert not polynomial(words.Q_expr(active))
            result.append(str(pairing(functional, polynomial(active))))
        return result

    original_values = coordinates(original)
    assert original_values == {32:["-163296000","816480000"],34:["0","400982400"]}[args.level]
    return dict(parent_level=args.level, checked_letter_brackets=checked, raising_axis=3,
                original_terms=len(original), lowered_terms=len(lowered), returned_terms=len(returned),
                original_pairings=original_values, returned_pairings=coordinates(returned))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", required=True, choices=(32, 34), type=int)
    parser.add_argument("--words", type=Path, default=WORDS)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--pure", type=Path, default=DATA/"su3_level32/descendant.json")
    parser.add_argument("--mixed", type=Path, default=DATA/"su3_bmn_sources/descendants/level34.json")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable checks on the letter brackets and pairings")
    if args.source is None:
        args.source = (DATA/"su3_level32/O32_C_even.json.gz" if args.level == 32
                       else HERE.parent/"su3_j34_20261118/O34_C_even.json.gz")
    args.out.write_text(json.dumps(calculate(args), indent=2)+"\n")
