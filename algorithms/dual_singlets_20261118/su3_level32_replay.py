"""Verify the SU(3) level-32 source, its pure scalar restriction and its level-36 descendant."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from su3_bmn_descendants import (WORDS, even_part, flavour_root, load_words, polynomial,
                                 read_operator, restrict_block, verify_classical, witness)
from su3_level36_scalar_channel import published_core
from verify import add_scaled, pairing


LEVEL32_SHA = "897fbf27bced5562684f9537ed1bd9e1b8a352f7149cc3b99ddb22c89b6c823e"
DATA = Path(__file__).resolve().parent/"taylor_hooks/su3_level32"


def restrictions(parent, words):
    assert even_part(parent, words) == parent and words.expr_sector(parent) == (5,8,4,4)
    pure = polynomial(restrict_block(parent, (4,0,0), words))
    core = polynomial(published_core(words))
    assert pure and pure == {key: -1020600*c for key,c in core.items()}
    descendant = {}
    for j,k,l in ((0,1,2), (1,2,0), (2,0,1)):
        component = parent if j == 0 else flavour_root(parent, j, 0, words)
        raised = words.Qplus_expr(words.Qplus_expr(component, l), k)
        add_scaled(descendant, raised, Fraction(1) if j == 0 else Fraction(1,4))
    assert words.expr_sector(descendant) == (7,8,5,5)
    assert not flavour_root(descendant, 0, 1, words) and not flavour_root(descendant, 1, 2, words)
    assert even_part(descendant, words) == descendant
    active = restrict_block(descendant, (4,0,0), words)
    quartic = witness(active, (7,4,5,5), words)
    mixed = restrict_block(descendant, (3,0,0), words, 5)
    mixed_value = polynomial(mixed)
    assert not polynomial(words.Q_expr(mixed))
    return (dict(source_terms=len(parent), quartic_scalar_coordinates=len(pure), core_coefficient=-1020600,
                 descendant_terms=len(descendant), descendant_weight=[8,5,5], descendant_length=11,
                 quartic_descendant=quartic, mixed_active_coordinates=len(mixed_value)),
            polynomial(active), mixed_value)


def pair_descendants(parent, words, level34_record, core_record):
    data, pure, mixed = restrictions(parent, words)
    def evaluate(functional, expression):
        return pairing({(i,j):Fraction(c) for i,j,c in functional}, expression)
    first = core_record["twice_raised_core"]["functional"]
    second = level34_record["C_even_length11_raised"]["functional"]
    row1 = [evaluate(first, pure), Fraction(0)]
    row2 = [evaluate(second, mixed), Fraction(level34_record["C_even_length11_raised"]["pairing"])]
    assert core_record["level34_pure_four_scalar"]["active_coordinates"] == 0
    determinant = row1[0]*row2[1]-row1[1]*row2[0]
    assert str(determinant) == core_record["determinant"]
    data.update(pairings=[[str(v) for v in row] for row in (row1,row2)], determinant=str(determinant))
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("classical","descendant"))
    parser.add_argument("--source", type=Path, default=DATA/"O32_C_even.json.gz")
    parser.add_argument("--words", type=Path, default=WORDS)
    parser.add_argument("--engine", type=Path)
    parser.add_argument("--quartic", type=Path, default=DATA.parent/"su3_bmn_sources/quartic.json")
    parser.add_argument("--level34", type=Path, default=DATA.parent/"su3_bmn_sources/descendants/level34.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or (args.mode=="classical" and not args.engine):
        parser.error("enable assertions and supply --engine for the classical variation")
    parent = read_operator(args.source, LEVEL32_SHA)
    words = load_words(args.words)
    if args.mode == "classical":
        record = verify_classical(parent, words, args.engine)
    else:
        record = pair_descendants(parent, words, json.loads(args.level34.read_bytes()),
                                  json.loads(args.quartic.read_bytes()))
    args.out.write_text(json.dumps(record, indent=2) + "\n")
