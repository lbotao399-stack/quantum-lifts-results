"""Calculate the classical determinant-pair homology and its rational multiplicity bounds."""
import argparse
import json
from pathlib import Path

from currents import PairComplex
from determinant import calculate as rational_weight

HERE=Path(__file__).resolve().parent

def run(colour, method):
    from colour_quotient import calculate, dominant, weyl
    data=json.loads((HERE/'colour.json').read_bytes())
    pair=PairComplex(5,8,data,5)
    if method=='augmented':
        other=PairComplex(6,8,data,4)
        return calculate(pair,other,colour,[127,131])
    weights={}
    multiplicity=0
    for (a,b),sign in weyl().items():
        w=dominant((colour[0]+1-a,colour[1]+1-b))
        if w not in weights:
            weights[w]=rational_weight(pair,w)
        multiplicity+=sign*weights[w]['homology']
    assert multiplicity>=0
    return dict(colour=colour,weight_spaces=list(weights.values()),multiplicity=multiplicity)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('method',choices=('rational','augmented'))
    parser.add_argument('colour',nargs=2,type=int)
    args=parser.parse_args()
    if not __debug__:parser.error('Enable assertions')
    print(json.dumps(run(tuple(args.colour),args.method),indent=2))
