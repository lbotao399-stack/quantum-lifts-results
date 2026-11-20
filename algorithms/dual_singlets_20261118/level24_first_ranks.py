"""Replay the level-24 order-one singlet map in eight and nine orthogonal colours."""
import argparse
import hashlib
from pathlib import Path

from flint import fmpz_mat, nmod_mat
import mixed_incoming as m
import terminal_interval_ranks as terminal
from rank_four_s5 import decode_rows


def blocks(colours):
    records, left, right = [], [], []
    previous, previous_scales = [], []
    for weight in ((1,1,3,3,3),(2,2,2,2,2)):
        length=sum(weight)
        formal=m.prepare(weight,length-1)
        sources=decode_rows(formal,"sources")
        classical=decode_rows(formal,"images")
        targets=decode_rows(m.prepare(weight,length),"sources")
        assert all(not m.q(row) for row in targets)
        quantum=[terminal.first(row,colours,"SO") for row in sources] if not records else []
        scales=terminal.row_scales(quantum) if quantum else [1]*len(sources)
        a,b,frame=terminal.coordinates(targets,classical,previous,scales,previous_scales,
                                       weight,colours,"SO")
        left.append(a)
        right.append(b)
        records.append(dict(**frame,weight=weight,formal_parents=len(sources)))
        previous,previous_scales=quantum,scales
    return records,left[0],right[1],left[1]


def ranks(first,advanced,second):
    classical=[fmpz_mat(first).rank(),fmpz_mat(second).rank()]
    full=[a+b for a,b in zip(first,advanced,strict=True)]
    empty=[0]*len(first[0])
    full.extend(empty+b for b in second)
    coupled=fmpz_mat(full).rank()
    modular=[]
    for prime in (1000003,1000033):
        values=[nmod_mat(matrix,prime).rank() for matrix in (first,second,full)]
        assert values==classical+[coupled]
        modular.append(dict(prime=prime,ranks=values))
    return dict(classical=classical,coupled=coupled,first_pairings=coupled-sum(classical),
                finite_field_controls=modular)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--colours",required=True,type=int,choices=(8,9))
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    assert __debug__
    records,first,advanced,second=blocks(args.colours)
    result=ranks(first,advanced,second)
    assert result["first_pairings"]=={8:1,9:2}[args.colours]
    result.update(colours=args.colours,source_charges=[1,1,8,8,8],target_charges=[0,0,8,8,8],
                  level=24,source_grade=1,target_grade=2,frames=records,
                  integer_blocks_sha256=[hashlib.sha256(m.compact(a)).hexdigest()
                                         for a in (first,advanced,second)])
    m.write(args.out,result)
    print(result["classical"],result["coupled"],result["first_pairings"],flush=True)


if __name__=="__main__":
    main()
