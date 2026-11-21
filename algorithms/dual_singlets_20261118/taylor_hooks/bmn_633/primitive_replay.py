"""Verify the rational symplectic primitive on a character-sized integer frame."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import mixed_incoming as m
from dual633_pairing import integral
from gaugino_primary_order import charge
from sp3_omega_actions import signed_frame
from sp3_omega_replay import pair
from symplectic_primitive import HERE, PAGE, WEIGHT, source_and_target
from trace_descent import action, evaluate, settings


def load_primitive(path):
    page, basis, polynomial = source_and_target()
    receipt = json.loads(gzip.decompress(path.read_bytes()))
    assert receipt["page_sha256"] == hashlib.sha256(PAGE.read_bytes()).hexdigest()
    assert page["outgoing"]["determinant"] == receipt["modular_minor_determinant"]
    entries = receipt["coefficients_in_outgoing_pivots"]
    assert len(entries) == len(set(i for i,c in entries)) == receipt["modular_coefficients"]
    primitive = Counter()
    for i, text in entries:
        assert 0 <= i < len(basis)
        coefficient = Fraction(text)
        for w,c in basis[i].items():
            primitive[w] += coefficient*c
    primitive = {w:c for w,c in primitive.items() if c}
    assert charge(primitive) == (WEIGHT, 9)
    assert (len(primitive), integral(primitive)[1]) == (receipt["primitive_terms"],receipt["primitive_denominator"])
    return primitive, polynomial


def replay(path, seed, work, workers, project_first):
    settings()
    primitive, polynomial = load_primitive(path)
    difference = Counter(m.q(primitive))
    difference.subtract(polynomial)
    residual, den = integral({w:c for w,c in difference.items() if c})
    frames = signed_frame(WEIGHT, 10, seed)
    values = evaluate([residual], frames, 6, "Sp")[0]
    assert len(values) == 1832 and not any(values)
    result = dict(primitive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), seed=seed,
                  primitive_terms=len(primitive), product_terms=len(polynomial),
                  formal_residual_terms=len(residual), formal_residual_denominator=den,
                  full_frame={degree:dict(dimension=data["dimension"], minor_mod_1000003=data["determinant_mod_prime"])
                              for degree,data in frames.items()}, exact_zero_coordinates=len(values))
    if project_first:
        image = action(primitive, 1, 6, "Sp")
        assert charge(image) == ((1,1,5,2,2),9)
        lowered = HERE.parents[1] / "sp3_omega/lowered.json.gz"
        page = json.loads(gzip.decompress(lowered.read_bytes()))
        library = m.library(work,6,m.P)
        result.update(lowered_page_sha256=hashlib.sha256(lowered.read_bytes()).hexdigest(),
                      first_action_terms=len(image), first_action_denominator=integral(image)[1],
                      prime=m.P, first_action_pairing=pair(image,page,library,workers))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--primitive", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--seed", type=int, default=1738401)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--project-first", action="store_true")
    args = p.parse_args()
    if not __debug__: p.error("enable assertions")
    args.work.mkdir(parents=True, exist_ok=True)
    result = replay(args.primitive,args.seed,args.work,args.workers,args.project_first)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
