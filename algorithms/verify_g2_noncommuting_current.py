"""Replay Chang--Lin (2.6) and (2.7)--(2.10) independently on six noncommuting jets.

The imaginary unit is SymPy I. A single exterior coefficient eta occurs once
in every surviving trace; the integer matrix records its coefficient. No two
odd input letters are simultaneously nonzero. No physical current is assumed
in computing the two variations.
"""
import argparse
import json
import sys
from pathlib import Path
from itertools import combinations, permutations
from sympy import I, Matrix, Rational, diag, zeros

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor'))
from fortuity.full.letters import Q_letter


def seven(x):
    assert x.shape==(3,3) and x.trace()==0
    return diag(x,-x.T,0)


def comm(x,y):return x*y-y*x


def verify_embedding(matrices):
    # Volume forms on 3 and 3*, and their evaluation two-form wedged with
    # the fixed seventh covector. Verify the infinitesimal stabilizer on all
    # triples independently of the block Lie-bracket construction.
    form={}
    for (i,j,k),v in [((0,1,2),1),((3,4,5),1),*(( (r,3+r,6),1) for r in range(3))]:
        for p in permutations(range(3)):
            form[tuple((i,j,k)[s] for s in p)]=v*(-1)**sum(p[u]>p[t] for u,t in combinations(range(3),2))
    for x in matrices:
        for i,j,k in combinations(range(7),3):
            assert sum(x[l,i]*form.get((l,j,k),0)+x[l,j]*form.get((i,l,k),0)+x[l,k]*form.get((i,j,l),0) for l in range(7))==0
    return len(matrices)*35


def run():
    a=zeros(3);a[0,1]=1;b=a.T;h=comm(a,b);A,B,H=map(seven,(a,b,h));z=zeros(7)
    assert comm(A,B)==H and (H*H).trace()==4
    embedding=verify_embedding((A,B,H));rows=[]
    for s in range(3):
        j,k=(s+1)%3,(s+2)%3;orientation=1 if j<k else -1
        for spin,axis in zip(('+','-'),((1,0),(0,1))):
            scalar_j=(0,0,(j,));scalar_k=(0,0,(k,));gauge=(*axis,());matter=(0,0,tuple(sorted((j,k))))
            # (2.7)-(2.8) in increasing theta order, applied before Q_letter.
            coeff={scalar_j:-2*I*A,scalar_k:-2*I*B,gauge:-I*H}
            def q(l):return sum((c*coeff.get(u,z)*coeff.get(v,z) for c,u,v in Q_letter(l)),z)
            assert q(scalar_j)==q(scalar_k)==q(gauge)==z
            qj=q((*axis,(j,)))/(-2*I);qk=q((*axis,(k,)))/(-2*I);qm=q(matter)/(-2*I*orientation)
            assert qj==-I*comm(H,A) and qk==-I*comm(H,B) and qm==-2*I*comm(A,B)
            scalar=((A*qk-B*qj)/2).trace();fermion=(qm*H).trace()
            # Check again without the Taylor routines, using (2.6) on physical fields.
            assert scalar==I*(H*comm(A,B)).trace()==4*I
            assert fermion==-2*I*(comm(A,B)*H).trace()==-8*I
            # Both archived trace terms have unit coefficient in their Taylor variables;
            # the sorted theta-pair orientation is already in its archived coefficient.
            arch_scalar=(coeff[scalar_j]*q((*axis,(k,)))-coeff[scalar_k]*q((*axis,(j,)))).trace()/2
            arch_fermion=orientation*(q(matter)*coeff[gauge]).trace()
            assert arch_scalar==-4*scalar and arch_fermion==-2*fermion and arch_scalar+arch_fermion==0
            assert arch_scalar-arch_fermion==-32*I
            rows.append(dict(flavour=s+1,dotted_spin=spin,physical_scalar_half_variation_over_i_eta=int(scalar/I),
                             physical_ordered_psi_lambda_variation_over_i_eta=int(fermion/I),
                             closed_matter_to_scalar_half_coefficient=str(-scalar/fermion),
                             archive_scalar_variation_over_i_eta=int(arch_scalar/I),archive_matter_variation_over_i_eta=int(arch_fermion/I),
                             archive_sum_over_i_eta=0,opposite_relative_sign_over_i_eta=-32))
    assert len(rows)==6 and {r['closed_matter_to_scalar_half_coefficient'] for r in rows}=={str(Rational(1,2))}
    return dict(representation='sl3 inside g2 in 3 + 3* + 1; rho(X)=diag(X,-X.T,0)',probe='cyclic scalars E12,E21; gaugino eta*[E12,E21]; remaining fields zero',
                trace_commutator_squared=4,invariant_three_form_component_checks=embedding,six_components=rows,
                field_order='matter fermion before gaugino',scope='Necessary and sufficient relative sign in the two-term physical ansatz with the Chang-Lin component transformations. The last column tests the opposite relative sign in the same Taylor dictionary, not the differential of the Bomans et al. current in its own conventions.',instanton_matrix_element=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);args=p.parse_args();assert __debug__
    text=json.dumps(run(),indent=2)+'\n';print(text,end='')
    if args.output:args.output.write_text(text)
