"""Identify both common eight-letter singlets and test the physical currents on two Cartan jets."""
import argparse
import json
import sys
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path

import mixed_incoming as m
from mixed33_exact import apply, read, sha

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import verify_g2_superfield_sign as sf
from verify_g2_s3_currents import add, components, exterior, scale, trace
from fortuity.full.letters import Q_letter, letters_up_to


def superfields():
    currents = {(s,i): sf.graviton(s,i) for s,i in product(range(3),range(2))}
    bilinear = sf.canonical({(((1,0,()),(0,1,())),): 1})
    lg = add(*(scale(sf.mult(currents[s,i][0],currents[s,1-i][1]),Fraction((-1)**i,6))
               for s,i in currents))
    lll = add(*(scale(sf.mult(currents[s,i][0],currents[t,j][0],
                              sf.even_derivative(currents[u,1-j][0],1-i)),
                      sf.epsilon(s,t,u)*(-1)**(i+j))
                for s,t,u in permutations(range(3)) for i,j in product(range(2), repeat=2)))
    return sf.mult(bilinear,bilinear,lg), sf.mult(bilinear,lll)


def cartan(group):
    if group == "SO(7)":
        return [{(2*i,2*i+1): 1, (2*i+1,2*i): -1} for i in range(3)]
    assert group == "USp(6)"
    return [{(i,i): 1, (i+3,i+3): -1} for i in range(3)]


def field(matrices, masks):
    return {k: {mask: a} for matrix,mask in zip(matrices,masks) for k,a in matrix.items()}


def matmul(left, right):
    out = {}
    for (i,j), u in left.items():
        for (j2,k), v in right.items():
            if j == j2:
                out[i,k] = add(out.get((i,k),{}), exterior(u,v))
    return {k: v for k,v in out.items() if v}


def evaluate(expr, jets):
    result = {}
    for w,c in expr.items():
        value = {0: c}
        for a,b in w:
            value = exterior(value, trace(jets.get(a,{}), jets.get(b,{})))
        result = add(result, value)
    return result


def row(group, matter, expressions):
    cols = cartan(group)
    plus, minus = (1,0,()), (0,1,())
    jets = {plus:field(cols,(1,2,4)), minus:field(cols,(8,16,32)),
            (0,0,(0,)):field(cols[:1],(0,))}
    if matter:
        jets[0,0,(1,2)] = field(cols[:1],(64,))
        mask = 127
    else:
        jets[0,0,(1,)] = field(cols[1:2],(0,))
        jets[1,0,(2,)] = field(cols[:1],(0,))
        mask = 61
    letters = letters_up_to(3,3)
    for l in letters:
        result = {}
        for c,a,b in Q_letter(l):
            for pos,v in matmul(jets.get(a,{}),jets.get(b,{})).items():
                result[pos] = add(result.get(pos,{}), scale(v,c))
        assert not any(result.values())
    s = trace(jets[plus],jets[minus])
    def physical(sign):
        return [exterior(s,p).get(mask,0) for p in components(jets, sign)]
    check = [evaluate(e,jets).get(mask,0) for e in expressions]
    assert check == physical(1)
    return [list(map(str,physical(-1))), list(map(str,check))], len(letters), mask


def run():
    folder = Path(__file__).resolve().parent/"mixed33"
    cert = read(folder/"formal_certificate.json.gz")
    plan = read(folder/"outgoing.json.gz")
    assert sha(folder/"outgoing.json.gz") == cert["outgoing_plan"]
    raw = list(map(m.decode, plan["source_dictionary"]))
    sources = list(map(dict, plan["sources"]))
    archived = [{raw[j]:v for j,v in apply(dict(rec),sources).items()}
                for rec in cert["formal_class_recipes"]]
    a,b = superfields()
    assert a == scale(archived[1], Fraction(-1,12)) and b == scale(archived[0], -1)
    assert not m.q(a) and not m.q(b)
    output = {}
    for group, d1rank in (("SO(7)",1), ("USp(6)",0)):
        rows = [row(group,matter,(a,b)) for matter in (True,False)]
        physical, classical = ([list(map(Fraction,v[0][k])) for v in rows] for k in (0,1))
        det = lambda xs: xs[0][0]*xs[1][1] - xs[0][1]*xs[1][0]
        assert det(physical) and det(physical) == -det(classical)
        tag = "SO_odd" if group == "SO(7)" else "Sp"
        dimensions = {r["group"]:r for r in read(folder.parent/"charge_maps.json")
                      if r["mu"] == [3,3,1,1,1]}
        ds = dimensions[tag]["by_length_and_fermion_degree"]
        middle, terminal = (sum(ds[str(i)].values()) for i in (8,9))
        for p in (1000003,1000033):
            test = read(folder/f"classes_{tag}_{p}.json")
            assert test["rank_prefixes"][-1] == [227+d1rank]*2
            assert read(folder/f"result_{tag}_{p}.json")["rank"] == 225
        dim = middle - terminal - 225
        assert dim - d1rank == 2
        output[group] = dict(physical_matrix=[list(map(str,r)) for r in physical],
                             physical_determinant=str(det(physical)),
                             superfield_matrix=[list(map(str,r)) for r in classical],
                             exterior_masks=[r[2] for r in rows], Q0_letter_checks=[r[1] for r in rows],
                             classical_dimension=dim, protected_rank=2,
                             previously_certified_outgoing_d1_rank=d1rank,
                             remaining_unaccounted_outgoing_sources=0)
    return dict(mu=plan["mu"], length=8, doubled_charges=plan["source_sector"],
                formal_class_trace_terms=list(map(len,archived)),
                current_coordinates_in_formal_class_order=[["0","-1/12"],["-1","0"]],
                results=output, instanton_amplitude_evaluated=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert __debug__
    certificate = run()
    args.out.write_text(json.dumps(certificate, indent=2)+"\n")
    print({g:z["physical_determinant"] for g,z in certificate["results"].items()})
