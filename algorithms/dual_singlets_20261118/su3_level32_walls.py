"""Find one three-colour correction satisfying the level-32 equation on three scalar loci."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys

from su3_level32_replay import DATA

sys.path.insert(0, str(DATA.parents[2] / "su3_independent_prime_20261118/runtime/work/sp3-native.nosync"))
import sp_environment
import block_first
import block_reduce
import collision_first
import collision_solve
import su_trace_adapter_v1 as algebra
from g2_full_classical_screen_v1 import contents
from replay_collision import derivative


def right_hand_sides(folder):
    manifest = json.loads((folder / "extraction.json").read_bytes())
    assert manifest["rhs_sha256"] == "f7810267d100f7cc9c9b4c93e045540ab2291df8b8bb35eeb4ad77feddd47a8a"
    assert manifest["denominator"] == 3 and manifest["sign_applied_to_input"] == -1
    rhs = defaultdict(lambda: defaultdict(dict))
    for case in manifest["rows"]:
        data = gzip.decompress((folder / case["file"]).read_bytes())
        assert hashlib.sha256(data).hexdigest() == case["trace_sha256"]
        assert len(data.splitlines()) == case["restricted_trace_words"]
        for value, raw in map(json.loads, data.splitlines()):
            word = tuple(map(tuple, raw))
            rhs[case["scalar"], case["power"]][collision_solve.content(word)][word] = value
    return rhs


def local_quartic(scalar, power, word):
    # A derivative letter is active on all three scalar loci. The four-trace
    # insertion then vanishes on the other two by Cayley-Hamilton in size two.
    i = next(i for i,tr in enumerate(word) if any(block_first.letter(k)[:2] != (0,0) for k in tr))
    pair = (scalar,scalar)
    rest = (pair,)*(power//2-1)
    den = 6**(power//2-1)
    return [(word[:i] + (pair+word[i],) + word[i+1:] + rest, Fraction(-1,2*den)),
            (word + (pair,) + rest, Fraction(1,4*den))]


def coupled(folder, quartic=False):
    rhs = right_hand_sides(folder)
    domain = sorted(w for raw in contents((1, 1, 7, 3, 3), 10)
                    for w in collision_first.pair_words(tuple(sorted(map(block_first.code, raw)))))
    assert len(set(domain)) == len(domain)
    extra = {}
    equations, labels, blocks = [], [], []
    for scalar in (1, 2, 4):
        projected = defaultdict(dict)
        for i, w in enumerate(domain):
            power, factor, reduced = collision_first.restrict_word(w, scalar)
            if factor:
                assert power % 2 == 0 and reduced not in projected[power] and factor == 6**(power//2)
                projected[power][reduced] = (i, factor)
        for power in sorted(set(projected) | {n for s,n in rhs if s==scalar}):
            images = defaultdict(lambda: defaultdict(dict))
            for w, (i, factor) in projected[power].items():
                added = None
                if quartic and power > 0:
                    added = len(domain) + len(extra)
                    extra[added] = (scalar,power,w)
                for target, c in collision_solve.differential(w, scalar).items():
                    images[collision_solve.content(target)][i][target] = factor*c
                    if added is not None:
                        images[collision_solve.content(target)][added][target] = c
            before = len(equations)
            for t in sorted(set(images) | set(rhs[scalar,power])):
                diagonal, rows, target = collision_solve.block_data(rhs[scalar,power][t], images[t], t)
                for k in collision_solve.compress(rows, target):
                    equations.append(({i:row[k] for i,row in rows.items() if k in row}, target.get(k,0)))
                    labels.append((scalar,power,t,diagonal,*k))
            blocks.append([scalar,power,len(projected[power]),len(equations)-before])
    kind, coordinates, value = collision_solve.eliminate(equations)
    result = dict(source="SU(3) level 32", full_pair_words=len(domain), equations=len(equations),
                  blocks=blocks, kind=kind)
    if quartic:
        result["quartic_combinations"] = len(extra)
    if kind == "solution":
        full = Counter()
        for i,c in coordinates.items():
            terms = [(domain[i],1)] if i < len(domain) else local_quartic(*extra[i])
            for w, coefficient in terms:
                full[block_reduce.to_fort(w)] += c*coefficient/3
        # The negative-transpose colour involution commutes with restriction
        # and with Q0. Keep the charge conjugation of the input representative.
        full = {w:v/2 for w,v in algebra.parity_project(full,1).items()}
        assert algebra.conjugate(full) == full
        result.update(selected_rank=value, correction=[[block_reduce.from_fort(w),str(c)] for w,c in sorted(full.items())])
    else:
        result.update(pairing=str(value/3), functional=[[labels[i],str(c)] for i,c in sorted(coordinates.items())])
    return result


def verify(folder, certificate):
    if certificate["kind"] != "solution":
        raise ValueError("a three-colour correction is required")
    rhs = right_hand_sides(folder)
    solution = {tuple(map(tuple, raw)):Fraction(c) for raw,c in certificate["correction"]}
    cases = []
    for scalar in (1, 2, 4):
        residual = defaultdict(lambda: defaultdict(Counter))
        for (s,n), groups in rhs.items():
            if s == scalar:
                for t,row in groups.items():
                    residual[n][t].update({w:-Fraction(c,3) for w,c in row.items()})
        for word, c in solution.items():
            power, factor, active = collision_first.restrict_word(word, scalar)
            if factor:
                for w, v in derivative(active, scalar).items():
                    residual[power][collision_solve.content(w)][w] += factor*c*v
        for power, groups in sorted(residual.items()):
            terms = 0
            for group in groups.values():
                row = {w:c for w,c in group.items() if c}
                assert not block_reduce.polynomial(row)
                terms += len(row)
            cases.append(dict(scalar=scalar, power=power, contents=len(groups), formal_residual_words=terms,
                              nonzero_polynomial_coordinates=0))
    return dict(correction_words=len(solution), generic_traceless_upper_block=True, cases=cases)


def verify_pairs(folder, certificate):
    rhs = right_hand_sides(folder)
    tests = defaultdict(lambda: defaultdict(dict))
    for (scalar,power,t,diagonal,even,odd), coefficient in certificate["functional"]:
        assert diagonal == collision_solve.diagonal_for(t)
        tests[scalar,power][tuple(t)][even,odd] = Fraction(coefficient)
    boundary, target = Counter(), Fraction(0)
    domains = set()
    for (scalar,power), groups in sorted(tests.items()):
        local, _ = collision_first.parent_words(scalar,power,(1,1,7,3,3))
        factor = 6**(power//2)
        for t,ell in groups.items():
            coords = collision_solve.sum_poly(rhs[scalar,power][t],collision_solve.diagonal_for(t))
            target += sum((ell.get(k,0)*v for k,v in coords.items()), Fraction())/3
        for parent in local:
            _, lifted = block_first.canonical_multi(parent+((scalar,scalar),)*(power//2))
            domains.add(lifted)
            for word,c in collision_solve.differential(parent,scalar).items():
                t = collision_solve.content(word)
                if t in groups:
                    ell = groups[t]
                    coords = collision_solve.polynomial(word,collision_solve.diagonal_for(t))
                    boundary[lifted] += factor*c*sum(ell.get(k,0)*v for k,v in coords.items())
    assert not any(boundary.values()) and target == Fraction(certificate["pairing"]) and target
    return dict(pair_words=certificate["full_pair_words"], pair_words_on_supported_loci=len(domains),
                functional_terms=len(certificate["functional"]), nonzero_boundary_pairings=0,
                negative_first_image_pairing=str(target))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("pairs", "quartic", "verify", "verify-pairs"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--certificate", type=Path, required=True)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    if not __debug__ or (a.mode.startswith("verify") and not a.out):
        p.error("enable assertions and pass --out for verify")
    if not a.mode.startswith("verify"):
        result = coupled(a.work, a.mode=="quartic")
        a.certificate.write_text(json.dumps(result, separators=(",", ":")) + "\n")
        print(result["kind"], result["equations"], result.get("selected_rank"), flush=True)
    else:
        replay = verify if a.mode == "verify" else verify_pairs
        result = replay(a.work, json.loads(a.certificate.read_bytes()))
        a.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result) if a.mode == "verify-pairs" else (result["correction_words"],len(result["cases"])), flush=True)
