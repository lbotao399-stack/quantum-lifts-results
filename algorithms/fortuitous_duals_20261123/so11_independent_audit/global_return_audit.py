#!/usr/bin/env python3
"""Independent global Fischer-functional audit of the SO11 factor25200.

This does not use the candidate's odd-return or determinant routines.
SymPy supplies determinant coefficients; the calculation sums physical
odd-marker derivatives on every row-degree(3,3,3) marker monomial.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,product
from math import factorial,prod
import ast,hashlib,json,time,argparse
import sympy as sp

HERE=Path(__file__).resolve().parent
SOURCE_COMMIT='9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9'
SOURCE_FILE='orthogonal_glueball_first_pair.py'
SOURCE_SHA256='eeb6d6d3acb94df435bbac3424c695abcd82156c569e2c9aa37be5676f553033'

def default_source():
    # Prefer the other reproducible package's complete source checkout; the
    # bundled unchanged source makes this audit independently portable.
    candidates=(
        HERE.parent/'positive_sector/inputs/algorithms/dual_singlets_20261118'/SOURCE_FILE,
        HERE/'inputs'/SOURCE_FILE,
    )
    return next((path for path in candidates if path.is_file()),candidates[-1])

def original_coefficient_class(source):
    raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256:
        raise ValueError('The Pfaffian source differs from the pinned, audited file')
    tree=ast.parse(raw,filename=str(source))
    selected=[node for node in tree.body
              if isinstance(node,(ast.FunctionDef,ast.ClassDef))
              and node.name in {'scalar_code','PfaffianCoefficient'}]
    assert len(selected)==2
    # No census, native evaluator, trace canonicalizer, or candidate return
    # routine is imported. These two upstream definitions are unchanged.
    namespace=dict(Counter=Counter,Fraction=Fraction,lru_cache=lru_cache,
                   permutations=permutations,product=product,
                   factorial=factorial,prod=prod,sign=permutation_sign)
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(source),'exec'),namespace)
    return namespace['PfaffianCoefficient']

def integer_partitions(total,minimum=1):
    if not total:
        yield ()
    for first in range(minimum,total+1):
        for rest in integer_partitions(total-first,first):
            yield (first,)+rest

def permutation_sign(order):
    return -1 if sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))%2 else 1

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pfaffian-source',type=Path,default=default_source(),
                        help='Unchanged pinned orthogonal_glueball_first_pair.py')
    parser.add_argument('--out-dir',type=Path,default=HERE)
    args=parser.parse_args();start=time.monotonic()
    if not __debug__:
        parser.error('Run with Python assertions enabled; do not use -O')
    variables=sp.symbols('x00 x01 x02 x10 x11 x12 x20 x21 x22')
    determinant=sp.Matrix(3,3,variables).det()
    coefficients={k:{degree:int(c) for degree,c in sp.Poly(determinant**k,*variables).terms()} for k in (3,5)}
    # The physical return order is A-,A+,B-,B+,C-,C+. For left odd
    # derivatives, the sign on eta_+0 eta_+1 eta_+2 eta_-0 eta_-1 eta_-2
    # is the permutation sign of the order in which variables are removed.
    assignments=[]
    for plus in permutations(range(3)):
        for minus in permutations(range(3)):
            increments=[0]*9;removal=[]
            for flavour in range(3):
                increments[3*flavour+minus[flavour]]+=1
                increments[3*flavour+plus[flavour]]+=1
                removal.extend((3+minus[flavour],plus[flavour]))
            s=permutation_sign(removal)
            assert s==permutation_sign(plus)*permutation_sign(minus)
            assignments.append((tuple(increments),s,removal))
    row_compositions=[v for v in product(range(4),repeat=3) if sum(v)==3]
    checks=[]
    for rows in product(row_compositions,repeat=3):
        degree=tuple(v for row in rows for v in row)
        terms=[]
        for increments,s,removal in assignments:
            original=tuple(a+b for a,b in zip(degree,increments))
            c=coefficients[5].get(original,0)
            terms.append(s*c*prod(factorial(n) for n in original))
        lhs=sum(terms)
        rhs=25200*coefficients[3].get(degree,0)*prod(factorial(n) for n in degree)
        assert lhs==rhs
        checks.append({'scalar_marker_exponents':degree,'physical_fischer_coefficient':lhs,'rhs_coefficient':rhs})
    assert len(checks)==1000 and sum(bool(r['rhs_coefficient']) for r in checks)==54
    # Audit the original implementation's factorial convention, rather than
    # assuming that coefficients of a determinant power equal its differential.
    PfaffianCoefficient=original_coefficient_class(args.pfaffian_source)
    factorial_checks=[]
    for power in (3,5):
        observed=PfaffianCoefficient(power,0).weights
        expected={sum(n<<(3*i) for i,n in enumerate(degree)):c*prod(factorial(n) for n in degree) for degree,c in coefficients[power].items()}
        assert observed==expected
        factorial_checks.append({'power':power,'all_weight_coefficients':len(expected),'matched':True})
    # Exercise the original coefficient method itself, bypassing its cache.
    # A degree15 nonzero trace monomial has at most seven trace factors.
    # These partitions cover every possible multiplicity pattern and a few
    # patterns that homogeneity would exclude, without changing the oracle.
    pf=PfaffianCoefficient(3,0)
    multiplicities=[partition for count in range(1,8) for partition in integer_partitions(count)]
    coefficient_method_checks=0
    for x,c in pf.weights.items():
        for partition in multiplicities:
            words=tuple((i,) for i,n in enumerate(partition) for _ in range(n))
            pf.traces_product=lambda words,x=x:{(63,x,0):1}
            actual=PfaffianCoefficient.coefficient.__wrapped__(pf,words)
            assert actual==-Fraction(c,prod(factorial(n) for n in partition))
            coefficient_method_checks+=1
    x=next(iter(pf.weights))
    for mask in range(63):
        pf.traces_product=lambda words,mask=mask:{(mask,x,0):1}
        assert PfaffianCoefficient.coefficient.__wrapped__(pf,((0,),))==0
    pf.traces_product=lambda words:{(63,x,1):1}
    assert PfaffianCoefficient.coefficient.__wrapped__(pf,((0,),))==0
    # Moving each odd marker past the preceding odd colour matrices creates
    # 15 exchanges. Grouping the six colour letters into any trace factors
    # separates that count into within-trace and between-trace exchanges.
    koszul=[]
    for cuts in product((False,True),repeat=5):
        lengths=[1]
        for cut in cuts:
            if cut:lengths.append(1)
            else:lengths[-1]+=1
        within=sum(n*(n-1)//2 for n in lengths)
        between=sum(a*b for i,a in enumerate(lengths) for b in lengths[i+1:])
        assert within+between==15
        koszul.append({'odd_letters_per_trace':lengths,'exchanges':within+between,'sign':-1})
    raw=(json.dumps(checks,separators=(',',':'))+'\n').encode()
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/'all_marker_functional_checks.json').write_bytes(raw)
    result={'passed':True,'method':'Independent SymPy determinant expansion and complete Fischer-functional transfer, not a selected-trace probe test','complete_scalar_marker_monomials':len(checks),'nonzero_monomials':54,'admissible_plus_minus_marker_assignments':len(assignments),'physical_return_order':['A-','A+','B-','B+','C-','C+'],'marker_order':['+0','+1','+2','-0','-1','-2'],'base_removal_order':[3,0,4,1,5,2],'base_removal_sign':1,'physical_fischer_factor':25200,'left_marker_insertion_factor':-25200,'marker_colour_koszul_sign':-1,'trace_factorizations_checked':len(koszul),'source_factorial_weights':factorial_checks,'source_coefficient_method_checks':coefficient_method_checks,'all_trace_multiplicity_patterns_up_to_seven_factors':len(multiplicities),'unwanted_masks_or_double_grade_rejected':64,'universal_total_letter_degree':15,'scalar_ancestor_flavour_degrees':[5,5,5],'returned_scalar_flavour_degrees':[3,3,3],'returned_ghost_degrees':[3,3],'conclusion':'P6 G_C = 25200 R_C globally in the formal graded trace algebra for PfaffianCoefficient(3,0) and the unchanged scalar order5 ancestor','coefficient_checks_sha256':hashlib.sha256(raw).hexdigest(),'source_commit':SOURCE_COMMIT,'source_path':'algorithms/dual_singlets_20261118/'+SOURCE_FILE,'upstream_file_sha256':SOURCE_SHA256,'source_definitions_loaded_unchanged_via_ast':['scalar_code','PfaffianCoefficient'],'seconds':time.monotonic()-start}
    (args.out_dir/'global_return_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
