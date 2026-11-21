"""Entrywise tests of d tr(phi dc) and the complete singlet minor, in archive Taylor variables."""
from fractions import Fraction as F
from itertools import product
import argparse,json
from pathlib import Path
from sympy import zeros
from verify_g2_noncommuting_current import seven,verify_embedding
import verify_g2_s3_currents as comp
import verify_g2_superfield_sign as sup


def run():
    assert __debug__;x=zeros(3);x[0,1]=1;X=seven(x);Y=seven(x.T)
    def sparse(m,coeff):return {(i,j):comp.scale(coeff,int(m[i,j])) for i,j in product(range(7),repeat=2) if m[i,j]}
    def evaluate(word,jets):
        # Independent coefficient-of-eta matrix polynomial; every probe word
        # is boson times fermion. Compare both word orders by an ordinary trace.
        value=0
        for traces,c in word.items():
            assert len(traces)==1 and len(traces[0])==2;u,v=traces[0]
            if u in jets and v in jets:value+=c*(jets[u]*jets[v]).trace()
        return value
    checks=[]
    for k,h in product(range(2),repeat=2):
        ak,ah=comp.AXES[k],comp.AXES[h];second=(ak[0]+ah[0],ak[1]+ah[1],());factor=1+int(k==h)
        scalar=(0,0,(0,));gaugino=(*ah,());diff_scalar=(*ak,(0,))
        l,_=sup.graviton(0,h);word=sup.even_derivative(l,k)
        for ds,dg in ((0,1),(1,0),(1,1)):
            ordinary={scalar:X,gaugino:Y,diff_scalar:ds*X,second:F(dg,factor)*Y};jet={key:sparse(X if key in (scalar,diff_scalar) else Y,{0:ds if key==diff_scalar else 1} if key in (scalar,diff_scalar) else {1:F(dg,factor) if key==second else 1}) for key in ordinary}
            val=comp.differentiated_jet(jet,0,k,h);assert val=={1:2*(ds+dg)} and evaluate(word,ordinary)==val[1]
            checks.append(dict(derivative=k,current_spin=h,scalar_derivative=ds,gaugino_derivative=dg,coefficient=int(val[1])))
    # The two rows are independently assembled as trace components and as
    # differentiated canonical superfield words; settle the full signed minor
    # in one convention without assuming the physical matter dictionary.
    minor=[]
    for matter,mask in ((True,31),(False,13)):
        values=comp.components(comp.diagonal_probe(matter),+1);row=[v.get(mask,0) for v in values];assert all(F(x).denominator==1 for x in row);minor.append(list(map(int,row)))
    old=sup.run();assert minor==old['superfield_six_letter_minor']==[[64,0],[16,576]]
    return dict(field='rationals with one exterior generator eta for the matrix probes',embedding_three_form_checks=verify_embedding((X,Y)),derivative_probes=checks,same_spin_only_gaugino_coefficient=checks[0]['coefficient'],archive_component_rows=minor,determinant=minor[0][0]*minor[1][1]-minor[0][1]*minor[1][0],classical_incoming_cohomology_dimension=2,instanton_matrix_element_computed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args();text=json.dumps(run(),indent=2)+'\n';print(text,end='');
    if a.output:a.output.write_text(text)
