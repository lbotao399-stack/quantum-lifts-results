"""Rank the pure gaugino singlet at bidegree (5,5), including its nine-letter sources."""
import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpz_mat

import mixed_incoming as m
from f4_formal_cycles import probe
from f4_singlet_characters import singlets
from rank_four_s5 import decode_rows
from mixed33_exact import read


def run(work, workers, prime, pure):
    work.mkdir(exist_ok=True, parents=True)
    lib = m.library(work, 27, prime)
    weight = (5, 5, 0, 0, 0)
    incoming = []
    records = []
    for length in range(2, 10):
        formal = m.prepare(weight, length)
        source = decode_rows(formal, "sources")
        dimension = singlets("F4", weight, length)
        words, rows = m.flatten(source + incoming)
        evaluated, frames = probe(rows, words, dimension, prime, lib, workers, 629500+111*length)
        source_values, boundary_values = np.split(evaluated,[len(source)])
        basis = m.minor(source_values,prime)
        boundary = m.minor(boundary_values,prime) if incoming and np.any(boundary_values) else {"rank":0}
        incoming = decode_rows(formal, "images")
        records.append(dict(length=length, cochains=dimension[length], basis=basis, incoming=boundary,
                            formal_sources=len(source), frames=frames))
        print("length",length,"cochains",dimension[length],"evaluated",basis["rank"],
              "incoming rank",boundary["rank"],flush=True)
    target = read(pure)
    assert target["weight"]==list(weight) and len(incoming)==target["formal_lower"]
    cochains, boundaries = target["formal_terminal"], target["formal_lower"]
    values=target["terminal_evaluations"]
    assert len(values)==cochains+boundaries+1
    rank=fmpz_mat(values[:cochains]).rank()
    image=values[cochains:-1]
    image_rank=fmpz_mat(image).rank()
    assert target["ranks"]==[fmpz_mat(target["lower_evaluations"]).rank(),rank,
                             image_rank,fmpz_mat(values[cochains:]).rank()]==[20,12,11,12]
    functional=[(i,Fraction(q)) for i,q in target["normalized_target_covector"]]
    assert [sum(row[i]*q for i,q in functional) for row in values[cochains:]]==[0]*boundaries+[1]
    assert rank==singlets("F4",weight,10)[10] and not singlets("F4",weight,11)
    records.append(dict(length=10,cochains=rank,basis={"rank":rank},incoming={"rank":image_rank},
                        exact_terminal_file=pure.name))
    for index, row in enumerate(records):
        outgoing = records[index+1]["incoming"]["rank"] if index+1<len(records) else 0
        row["classical_cohomology_upper"] = row["cochains"]-row["incoming"]["rank"]-outgoing
        assert row["classical_cohomology_upper"]>=0
    for previous, current in zip(records,records[1:]):
        assert current["classical_cohomology_upper"]==0 or previous["classical_cohomology_upper"]==0
    assert all(r["basis"]["rank"]==r["cochains"] for r in records)
    # Every differential has at least one endpoint with saturated consecutive ranks.
    for row in records:
        row["rational_cohomology"]=row["classical_cohomology_upper"]
    m.write(work/"pure_chain.json", dict(weight=weight, prime=prime, lengths=records))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--prime",type=int,default=1000003)
    p.add_argument("--pure",type=Path,default=Path(__file__).with_name("f4")/"pure.json.gz")
    assert __debug__
    run(**vars(p.parse_args()))
