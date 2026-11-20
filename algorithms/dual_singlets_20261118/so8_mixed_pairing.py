"""Project the eight-colour first actions onto the complete reflection-even terminal quotient."""
import argparse
from fractions import Fraction
from pathlib import Path

from flint import fmpq_mat

import mixed_incoming as m
from mixed33_exact import read
from rank_four_lowered_pair import run as first_action
from rank_four_orthogonal_paths import run as trace_relations


def run(source, classes, symplectic, work, workers):
    work.mkdir(parents=True, exist_ok=True)
    trace_relations(work, workers, colours=8)
    finite = read(work/"orthogonal_paths.json")
    formal = read(classes)
    assert finite["weight"] == formal["weight"] and finite["group"] == "SO(8)"
    assert finite["all_templates_checked"] or finite["classical"] == finite["cartan"]
    relations = [row["singlet"] for row in finite["relations"]]
    assert fmpq_mat(relations+formal["orthogonal_relations"]).rank() == finite["rank"]
    projections = [dict((i,Fraction(n)) for i,n in column)
                   for column in formal["weight_functionals"]]
    restrictions = [[sum(c*col.get(i,0) for i,c in row["raw"]) for col in projections]
                    for row in finite["relations"]]
    # Each of the five nine-colour relations restricts to an eight-colour relation.
    # The other two duals were constructed independently of the three Cartan coordinates.
    matrix = fmpq_mat([[str(x) for x in row] for row in restrictions])
    reduced, rank = matrix.rref()
    assert rank == finite["rank"]-len(formal["orthogonal_relations"])
    pivots = [next(i for i in range(2) if reduced[j,i]) for j in range(rank)]
    duals = []
    for i in range(2):
        if i in pivots:
            continue
        vec = [Fraction(0),Fraction(0)]
        vec[i] = 1
        for j,pivot in enumerate(pivots):
            vec[pivot] = -Fraction(str(reduced[j,i]))
        assert all(sum(x*y for x,y in zip(row,vec)) == 0 for row in restrictions)
        duals.append(vec)
    assert len(duals)+3 == finite["classical"]
    actions, projected = {}, []
    for channel in ("gauginos","scalars"):
        target = work/(channel+".json")
        first_action(source, classes, target, workers, channel=channel, cartan=symplectic, colours=8)
        rec = read(target)
        values = list(map(Fraction,rec["projected_first_action"]))
        assert rec["group"]=="SO(8)" and values[2:]==[0,0,0]
        surviving = [sum(x*y for x,y in zip(values[:2],vec)) for vec in duals]
        projected.append(surviving)
        actions[channel] = dict(five_formal_coordinates=list(map(str,values)),
                                eight_colour_coordinates=list(map(str,surviving)))
    rank = fmpq_mat([[str(x) for x in row] for row in projected]).rank() if duals else 0
    out = dict(group="SO(8)", parity="reflection-even", weight=finite["weight"], level=30,
               classical_dimension=finite["classical"], cartan_rank=finite["cartan"],
               extra_relations_on_nine_colour_duals=[[str(x) for x in row] for row in restrictions],
               surviving_duals=[[str(x) for x in v] for v in duals],
               actions=actions, first_differential_rank_on_tested_sources=rank,
               residual_cartan_kernel_upper=len(duals)-rank)
    m.write(work/"so8_mixed_pairing.json",out)
    print("eight colours",out,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).with_name("rank_four")
    p.add_argument("--source",type=Path,default=base/"source.json.gz")
    p.add_argument("--classes",type=Path,default=base/"classes.json.gz")
    p.add_argument("--symplectic",type=Path,default=base/"mixed_Sp.json")
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    args=p.parse_args()
    assert __debug__ and args.workers>0
    run(**vars(args))
