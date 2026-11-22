"""Replay the full second-page exclusions for SO(12), S^7 and SO(14), S^8."""
import argparse
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from orthogonal_glueball_even_span import span
from orthogonal_glueball_relation_span import relation_dimension
from orthogonal_glueball_second_pair import second_dual


def load(name):
    return json.loads(gzip.decompress((HERE/name).read_bytes()))


def check(colours, power, weights, args):
    original = load(f'second{colours}.json.gz')
    complete = load(f'span{colours}.json.gz')
    assert original['job']=='5693393' and complete['job']=='5693483'
    assert (original['colours'],original['power'])==(colours,power)
    rows = complete['basis_rows']
    dim = relation_dimension(colours,power,2)
    assert len(rows)==dim==complete['rank']==complete['balanced_flavour_relations']
    assert all(row == original['rows'][row['index']] for row in rows)
    assert complete['lower_full_weight_relations']==[relation_dimension(colours,power,k,False) for k in (0,1)]==[0,0]
    matrix = fmpq_mat(complete['minor'])
    assert matrix.nrows()==matrix.ncols()==dim and len(complete['probes'])==dim
    assert matrix.det() and str(matrix.det())==complete['determinant']
    coords = fmpq_mat([row['cartan'] for row in rows])
    augmented = fmpq_mat([row['cartan']+[row['second_score']] for row in rows])
    assert coords.rank()==augmented.rank()==complete['cartan_rank']==complete['augmented_rank']==len(weights)
    for row in rows:
        assert Fraction(row['second_score'])==sum(c*Fraction(value) for c,value in zip(weights,row['cartan']))
    assert original['dual']==second_dual(colours,power)
    if args.rebuild:
        rerun = span(original,args.workers,384 if colours==12 else 96,176 if colours==12 else 48)
        assert rerun=={k:v for k,v in complete.items() if k!='job'}
    print(f'SO({colours}), power {power}: {dim} independent relations; ranks {coords.rank()}, {augmented.rank()}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('positive workers and assertions required')
    check(12,7,[Fraction(-1,320),Fraction(-3,64),Fraction(-5,64)],args)
    check(14,8,[Fraction(1,1280)],args)
