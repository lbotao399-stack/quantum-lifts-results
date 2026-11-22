"""Rebuild all leading-relation contractions and compare the exact result."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
import json
from multiprocessing import get_context
from pathlib import Path
import sys

from flint import fmpq_mat

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from orthogonal_glueball_cartan_basis import matrix
from orthogonal_glueball_first_pair import PfaffianCoefficient
from orthogonal_glueball_third_pair import decode, predecessors, rows, targets


def contract(bounds):
    start,end=bounds
    data=rows(11,SUPPORT[start:end],LOOKUP,SCORES)
    pf=PfaffianCoefficient(3,0)
    hits=0
    first=Counter()
    third=Fraction()
    for encoded,entries,last in data:
        coefficient=pf.coefficient(decode(encoded))
        if coefficient:
            hits+=1
            third+=coefficient*Fraction(last)
            for i,n in entries:
                first[i]+=coefficient*Fraction(n)
    for method in (pf.trace,pf.aux_trace,pf.traces_product,pf.coefficient,pf.tilings,pf.exterior_sign):
        method.cache_clear()
    return [start,end,len(data),hits],first,third


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('assertions and a positive worker count are required')
    case=json.loads(gzip.decompress((HERE/'second.json.gz').read_bytes()))['cases']['11']
    expected=json.loads(gzip.decompress((HERE/'forward.json.gz').read_bytes()))
    basis,LOOKUP,SCORES=targets(6,case['dual'])
    SUPPORT=predecessors(set(LOOKUP)|set(SCORES))
    assert len(SUPPORT)==109525
    splits=[(i,min(i+70,len(SUPPORT))) for i in range(0,len(SUPPORT),70)]
    total=Counter()
    score=Fraction()
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        for row,(counts,first,third) in zip(expected['chunks'],pool.map(contract,splits),strict=True):
            assert counts==row[:4]
            total.update(first)
            score+=third
    trace=[str(total[i]) for i in range(len(basis))]
    assert trace==expected['first_trace'] and score==Fraction(expected['third_before_finite_correction'])
    orbit=(fmpq_mat([trace])*matrix(6)).tolist()[0]
    assert [str(v) for v in orbit]==expected['first_Cartan']
    print(len(SUPPORT),'leading words,',len(trace),'commuting coordinates, uncorrected score',score)
