"""Verify the lower-page quotient, the source correction, and both third pairings."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from flint import fmpq_mat

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
import mixed_incoming as traces
from orthogonal_glueball_adjoint_replay import parents
from orthogonal_glueball_cartan_basis import matrix
from orthogonal_glueball_relation_span import evaluate_relation, quotient_functional, relation_dimension
from orthogonal_glueball_second_pair import S, first_score, predecessors, second_dual, tree_columns
from orthogonal_glueball_second_sources import Relation, evaluate
from orthogonal_glueball_third_pair import decode
from rank_four_lowered_replay import Q1d_multi, terminal
from terminal_commutators import words


def read(name):
    return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))


def check_adjoints(colours,power,dual):
    pure=(S,)*power
    orientation,monomial=terminal(pure)
    first={w:Fraction(first_score(colours,w)) for w in predecessors([pure]) if first_score(colours,w)}
    second={decode(w):Fraction(v) for w,v in dual['scores']}
    sizes=[]
    for order,covector,lower in ((1,first,{pure:1}),(2,second,first)):
        columns=tree_columns(covector)|set(parents(lower,active_even=False))
        for v in columns:
            image=Fraction()
            for w,c in Q1d_multi(v,colours,'SO_odd').items():
                if order==1:
                    s,t=terminal(w)
                    if t==monomial:
                        image+=orientation*s*c
                else:
                    s,t=traces.so_canon(w)
                    image-=s*c*first.get(t,0)
            classical=sum(c*covector.get(w,0) for w,c in traces.q({v:1}).items())
            assert image==classical
        sizes.append(len(columns))
    assert sizes==[39,1362]
    return sizes


def rebuild_row(args):
    return evaluate_relation(*args)


def check_lower(case,full,workers):
    colours,power=case['dual']['colours'],case['dual']['power']
    span=read('relations'+str(colours))
    assert span['job']=='5693152'
    dim=relation_dimension(colours,power,2)
    assert span['grade_two_flavour_singlet_relations']==span['rank']==dim
    assert [relation_dimension(colours,power,k,False) for k in (0,1)]==span['lower_full_weight_relations']==[0,0]
    assert not span['nonzero_residuals']
    residual,definition=quotient_functional(case)
    assert definition==span['site_functional']
    cols=span['minor_columns']
    rows=span['independent']
    minor=fmpq_mat([[r['probe_values'][k] for k in cols] for r in rows]).det()
    assert minor and str(minor)==span['minor_determinant'] and len(rows)==len(cols)==dim
    assert all(Fraction(r['residual'])==0 for r in rows)
    if full:
        probes=[tuple(map(tuple,w)) for w in span['probes']]
        args=[(r['spec'],power,colours,probes,residual) for r in rows]
        with ProcessPoolExecutor(workers) as pool:
            for row,(values,value) in zip(rows,pool.map(rebuild_row,args),strict=True):
                assert values==row['probe_values'] and value==row['residual']
    dual=second_dual(colours,power)
    assert all(dual[k]==v for k,v in case['dual'].items())
    assert dual['first_factorization_columns']==33 and dual['second_factorization_columns']==1266
    return dict(colours=colours,power=power,finite_relations=dim,nonzero_minor=str(minor),
                adjoint_columns=check_adjoints(colours,power,dual))


def intervals(rows,n):
    assert rows[0][0]==0 and rows[-1][1]==n
    assert all(a<b for a,b,*_ in rows) and all(x[1]==y[0] for x,y in zip(rows,rows[1:]))


def check_third(case, rebuild):
    forward,reverse,first=read('forward'),read('adjoint'),read('first')
    assert (forward['job'],reverse['job'],first['job'])==('5693064','5693222','5693264')
    assert (first['parents'],first['nonzero'])==(181596,130488)
    assert len(first['first_trace'])==257 and first['first_trace']==forward['first_trace']
    intervals(forward['chunks'],forward['support'])
    intervals(reverse['blocks'],reverse['parents'])
    assert len(forward['chunks'])==1565 and forward['support']==109525
    case_rows=case['rows']
    coords=fmpq_mat([r['cartan'] for r in case_rows])
    scores=fmpq_mat([[r['second_score']] for r in case_rows])
    trace=fmpq_mat([forward['first_trace']])
    plane=trace*matrix(6)
    assert plane==fmpq_mat([forward['first_Cartan']])
    pivot_rows=[1,9,18]
    if rebuild:
        basis=words((4,4,2,2,2))
        for i in pivot_rows:
            saved=case_rows[i]
            spec=saved['spec']
            actual=evaluate((spec,case['dual']['scores'],11,6))
            assert actual==saved
            relation=Relation(6,spec)
            ordered=fmpq_mat([[str(relation.commuting_coefficient(w)) for w in basis]])
            assert ordered*matrix(6)==fmpq_mat([saved['cartan']])
            relation.clear()
    selected=fmpq_mat([case_rows[i]['cartan'] for i in pivot_rows])
    weights=fmpq_mat([[135,-90,-90]])
    assert weights*selected==plane and selected.rank()==coords.rank()==3
    correction=(weights*fmpq_mat([[scores[i,0]] for i in pivot_rows]))[0,0]
    augmented=fmpq_mat([r['cartan']+[r['second_score']] for r in case_rows])
    last=Fraction(forward['third_before_finite_correction'])
    last_row=forward['first_Cartan']+[str(-last)]
    assert [augmented.rank(),fmpq_mat(augmented.tolist()+[last_row]).rank()]==[3,4]
    assert str(correction)=='25920' and last==-17280
    value=last+Fraction(str(correction))
    assert reverse['nonzero']==sum(r[2] for r in reverse['blocks'])==12362
    assert value==Fraction(reverse['coordinate'])==sum(Fraction(r[3]) for r in reverse['blocks'])==8640
    return dict(group='SO(11)',level=36,first_coordinates=len(first['first_trace']),pivot_rows=pivot_rows,finite_correction=str(correction),
                uncorrected=str(last),third_coordinate=str(value))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-relations',action='store_true')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('assertions and a positive worker count are required')
    data=read('second')
    assert data['job']=='5693010'
    for n in ('11','13'):
        print(json.dumps(check_lower(data['cases'][n],args.rebuild_relations,args.workers)),flush=True)
    print(json.dumps(check_third(data['cases']['11'],args.rebuild_relations)),flush=True)
