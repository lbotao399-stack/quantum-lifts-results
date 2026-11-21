"""Transport the quantum relation to omega and compare its coefficient over the rationals."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from flint import nmod_mat
import mixed_incoming as m
from dual633_pairing import integral
from formal_product_correction import LOWER, grouped_basis
from full_extremal import project
from gaugino_primary_order import charge
from primitive_replay import load_primitive
from so7_bmn_633 import export_trace
from so8_level24_chain import solve_image
from so8_level24_descendant import check_letter, raising_images, substitute
from sp3_omega_actions import centre_basis, load_page
from sp3_omega_actions import signed_frame
from sp3_omega_substitutions import STEPS
from symplectic_primitive import HERE
from terminal_interval_ranks import sample
from trace_descent import action, evaluate, highest_weights, multiply_rows, settings


UPPER = (3,3,5,0,0)
OMEGA = HERE.parents[1] / "sp3_omega"


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def linear_combination(basis, entries):
    result = Counter()
    for i,text in entries:
        c = Fraction(text)
        for w,n in basis[i].items():
            result[w] += c*n
    return {w:c for w,c in result.items() if c}


def correction_from_archive(archive):
    data = read(archive)
    basis = grouped_basis(9)
    row = Counter()
    for block in data["blocks"]:
        label = tuple(map(tuple,block["trace_degrees"]))
        row.update(linear_combination(basis[label], block["first_coefficients"]))
    return {w:c for w,c in row.items() if c}


def return_polynomial(row):
    scaled, denominator = integral(row)
    for spin,flavour in reversed(STEPS):
        for letter in {l for w in scaled for tr in w for l in tr}:
            check_letter(letter,spin,flavour,True)
        scaled = substitute(scaled,spin,flavour,True)
    projected = project(scaled,m.sector(UPPER,9))
    assert charge(projected) == (UPPER,9)
    assert all(not any(z.values()) for z in raising_images(projected))
    return {w:Fraction(c,denominator) for w,c in projected.items()}


def returned_rows(primitive,archive):
    k,p = load_primitive(primitive)
    first = correction_from_archive(archive)
    assert m.q(first) == {w:-c for w,c in action(p,1,6,"Sp").items()}
    relation = Counter(first)
    relation.subtract(action(k,1,6,"Sp"))
    source = load_page(OMEGA/"target.json.gz")
    cycle = read(OMEGA/"cycle.json.gz")
    omega = linear_combination(centre_basis(source),cycle["cycle_basis_coefficients"])
    scaled, denominator = integral(omega)
    for spin,flavour in STEPS:
        scaled = substitute(scaled,spin,flavour,False)
    low = project(scaled,m.sector(LOWER,9))
    assert len(low)==63508
    return (return_polynomial({w:c for w,c in relation.items() if c}),
            return_polynomial({w:Fraction(c,denominator) for w,c in low.items()}))


def export(rows,folder):
    relations, omega = rows
    difference = Counter(relations)
    difference.update({w:112*c for w,c in omega.items()})
    records = {}
    for name,row in (("returned",relations),("omega",omega),("sum112",{w:c for w,c in difference.items() if c})):
        n,den,digest = export_trace(row,folder/(name+".txt"))
        records[name] = dict(terms=n,denominator=den,sha256=digest)
    (folder/"exports.json").write_text(json.dumps(records,indent=2)+"\n")


def rational_coefficient(returned, reference, out):
    page = load_page(OMEGA/"target.json.gz")
    rows = highest_weights(UPPER,8)
    parents = [r for degree in sorted(rows) for r in rows[degree]]
    basis = [m.q(parents[i]) for i in page["incoming"]["rows"]]
    assert len(basis)==1607
    ref, ref_den = integral(reference)
    right, right_den = integral(returned)
    rows = basis+[ref,right]
    words = sorted({w for row in rows for w in row})
    assert {m.nferm(w) for w in words} == {4}
    alphabet = sorted({l for content in m.contents(UPPER,9) for l in content})
    blocks = []
    for i in range(48):
        values = sample(words,4,6,20,974327+4099*i,"Sp",alphabet)
        blocks.append(multiply_rows(rows,words,values))
        if (i+1)*100 < len(basis)+1:
            continue
        table = np.column_stack(blocks)
        rank = int(nmod_mat((table[:-1] % m.P).tolist(),m.P).rank())
        print("upper cycle evaluations", table.shape, rank, flush=True)
        if rank==len(basis)+1:
            break
    else:
        raise ValueError("integer evaluations do not separate upper cycles")
    coeff, info = solve_image(table[:-1],table[-1].tolist())
    relative = coeff[len(basis)]*ref_den/right_den
    result = dict(weight=UPPER,lower_weight=LOWER,upper_cohomology=1,
                  upper_cochain_dimension=2689,upper_outgoing_rank=1081,
                  boundary_rank=len(basis),signed_cycle_space_rank=rank,
                  integer_coordinates=table.shape[1],field_for_minor=m.P,
                  returned_terms=len(right),reference_terms=len(ref),returned_denominator=right_den,
                  reference_denominator=ref_den,relative_coefficient=str(relative),
                  boundary_coefficients=[[page["incoming"]["rows"][i],str(v/right_den)]
                                         for i,v in sorted(coeff.items()) if i<len(basis)],
                  **info)
    m.write(out,result)
    return result


def replay_boundary(returned, reference, coefficients, seed):
    record = read(coefficients)
    rows = highest_weights(UPPER,8)
    parents = [row for n in sorted(rows) for row in rows[n]]
    primitive = linear_combination(parents,record["boundary_coefficients"])
    relative = Fraction(record["relative_coefficient"])
    residual = Counter(returned)
    residual.subtract({w:relative*c for w,c in reference.items()})
    residual.subtract(m.q(primitive))
    numerator,den = integral({w:c for w,c in residual.items() if c})
    frame = signed_frame(UPPER,9,seed)
    values = evaluate([numerator],frame,6,"Sp")[0]
    assert len(values)==2689 and not any(values)
    return dict(coefficients_sha256=hashlib.sha256(coefficients.read_bytes()).hexdigest(),
                relative_coefficient=str(relative),upper_primitive_terms=len(primitive),
                upper_primitive_denominator=integral(primitive)[1],formal_residual_terms=len(numerator),
                formal_residual_denominator=den,seed=seed,exact_zero_coordinates=len(values),
                integer_frame_minor_mod_1000003=frame[4]["determinant_mod_prime"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("primitive","formal","work","out"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--check",type=Path)
    parser.add_argument("--seed",type=int,default=2719409)
    args = parser.parse_args()
    if not __debug__: parser.error("enable assertions")
    settings()
    args.work.mkdir(parents=True,exist_ok=True)
    rows = returned_rows(args.primitive,args.formal)
    if args.check:
        result = replay_boundary(*rows,args.check,args.seed)
        args.out.write_text(json.dumps(result,indent=2)+"\n")
    else:
        export(rows,args.work)
        result = rational_coefficient(*rows,args.out)
    print("relative to the projected lowered omega:", result["relative_coefficient"],flush=True)
