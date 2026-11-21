"""Independent character, supercharge and integer-minor checks for the shortest SU(N) sectors."""
import argparse
from collections import Counter
from itertools import product
from pathlib import Path

import trace_descent as descent
from exact_dims import _multiplicity, partitions
from fortuity.full.letters import theta_sign
from fortuity.full.words import canon_word
from minimum_length_sources import source_bound
from saturated_flavour import highest_weight
from su3_charge_dimensions_v1 import parity_parts as three_colour_parts
from unitary_characters import parity_parts


def casimir_component(power, plus, minus, theta):
    # Supercoordinates stand to the left of coefficient matrices in C.
    terms=Counter()
    def append(word, used, used_plus, used_minus, coefficient):
        if len(word)==power:
            if (used_plus,used_minus,used)==(plus,minus,theta):
                sign,canonical=canon_word(word)
                if sign:
                    terms[(canonical,)]+=coefficient*sign
            return
        for i in range(plus-used_plus+1):
            for j in range(minus-used_minus+1):
                for choice in product((0,1),repeat=len(theta)):
                    additional=tuple(t for t,k in zip(theta,choice) if k)
                    sign=theta_sign(used,additional)
                    if sign:
                        append(word+((i,j,(0,)+additional),),tuple(sorted(used+additional)),
                               used_plus+i,used_minus+j,
                               coefficient*sign*(-1)**(len(used)*len(additional)))
    append((),(),0,0,1)
    return {w:c for w,c in terms.items() if c}


def checks():
    shapes=0
    for length in range(9):
        for bosons in range(length+1):
            for b,f in product(partitions(bosons),partitions(length-bosons)):
                assert parity_parts(b,f,3)==three_colour_parts(b,f)
                assert parity_parts(b,f,4)==tuple(_multiplicity(b,f,p,6) for p in (0,1))
                assert parity_parts(b,f,2)[1]==0
                shapes+=1
    descent.settings()
    classes=[]
    for size in (4,5):
        for power in range(2,size+1):
            for plus,minus,theta in ((0,0,()),(0,0,(1,2)),(1,0,(1,)),(1,1,(1,2))):
                expr=casimir_component(power,plus,minus,theta)
                assert expr
                group="SU-" if power%2 else "SU+"
                images=[descent.action(expr,n,size,group) for n in range(4)]
                assert not any(images)
                classes.append(dict(size=size,power=power,plus=plus,minus=minus,theta=theta,
                                    trace_terms=len(expr)))
    minors=[]
    for size,group,flavours,grade,seed in ((4,"SU-",(8,3,3),0,81173),
                                          (5,"SU+",(10,2,2),0,247599),
                                          (4,"SU-",(6,3,3),1,81173),
                                          (5,"SU-",(6,3,3),1,81173)):
        row=source_bound(flavours,0,grade,size,group,seed)
        predicted=highest_weight(size,tuple(row["weight"]),group)
        assert row["rational_cohomology_upper"]==predicted
        row["casimir_jet_multiplicity"]=predicted
        minors.append(row)
    return dict(species_shapes_per_independent_character_check=shapes,
                independent_character_checks=3*shapes,single_trace_components=classes,
                supercharge_images_tested=4*len(classes),source_minors=minors)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    assert __debug__
    data=checks()
    descent.m.write(args.out,data)
    print(data["independent_character_checks"],data["supercharge_images_tested"],
          [(r["source_dimension"],r["classical_rank_lower"],r["casimir_jet_multiplicity"])
           for r in data["source_minors"]])
