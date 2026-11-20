"""Necessary root-subalgebra conditions on the F4 level-30 one-fermion Cartan singlets."""
import argparse
from itertools import product
from pathlib import Path

from sympy import Matrix, Rational, eye, symbols

import mixed_incoming as m
from f4_cartan_sources import dimensions


def normal(vectors):
    h=Matrix.hstack(*vectors)
    assert h.shape==(4,3)
    return Matrix([
        (-1)**(i+1)*h.extract([j for j in range(4) if j!=i],range(3)).det() for i in range(4)])


def run():
    e=[eye(4)[:,i] for i in range(4)]
    roots=[sg*x for x in e for sg in (-1,1)]
    roots.extend(Matrix(v)/2 for v in product((-1,1),repeat=4))
    assert all(a.dot(a)==1 for a in roots)
    # The 24 roots form a Weyl-stable set in these Euclidean root coordinates.
    root_set={tuple(a) for a in roots}
    simple=(e[1]-e[2],e[2]-e[3],e[3],(e[0]-e[1]-e[2]-e[3])/2)
    for a in simple:
        assert {tuple(x-2*a.dot(x)*a/a.dot(a)) for x in roots}==root_set

    weight=(0,0,5,5,5)
    space=dimensions(weight,14)
    assert space=={1:2}
    # The degree-three classical closure test for a root sl2 has equal values
    # on its three diagonal coefficients and zero on all six off-diagonals.
    H=Matrix.diag(1,-1);E=Matrix([[0,1],[0,0]]);F=Matrix([[0,0],[1,0]])
    x=[H,E,F]
    jac=Matrix(3,3,lambda i,j:(x[j]*(x[(i+1)%3]*x[(i+2)%3]-x[(i+2)%3]*x[(i+1)%3])).trace())
    assert jac==2*eye(3)

    s=symbols("s")
    values=[]
    for name,root,h in (("short",e[3],e[:3]),
                        ("long",e[2]+e[3],[e[0],e[1],e[2]-e[3]])):
        v=normal(h)
        assert v==root and all(root.dot(z)==0 for z in h)
        trace=0;hessian=0;individual=[]
        for i in range(3):
            moved=h.copy();moved[i]=h[i]+s*root
            vs=normal(moved)
            a=(vs.dot(vs)**2*moved[i].dot(root)).expand().coeff(s,1)
            b=sum(30*(r.dot(vs))**4*r.dot(moved[i])*r.dot(root) for r in roots).expand().coeff(s,1)
            trace+=a;hessian+=b
            individual.append([str(a),str(b)])
        # A second evaluation uses homogeneity and the orthogonal projector
        # along the root. The root norm remains explicit for the long root.
        p4=sum(r.dot(root)**4 for r in roots)
        p6=sum(r.dot(root)**6 for r in roots)
        assert trace==3*root.dot(root)**3 and hessian==210*p6-120*root.dot(root)*p4
        values.append(dict(root=name,normal=list(map(int,root)),fourth_moment=str(p4),sixth_moment=str(p6),
                           probes=individual,trace_test=[str(trace),str(hessian)]))
    matrix=Matrix([[Rational(y) for y in v['trace_test']] for v in values])
    assert matrix==Matrix([[3,Rational(225,2)],[24,-360]]) and matrix.det()==-3780
    return dict(weight=weight,length=14,cartan_character_by_fermion_degree=space,
                short_root_count=len(roots),simple_reflections_checked=len(simple),
                sl2_trace_matrix=[[int(y) for y in row] for row in jac.tolist()],
                probes=values,cartan_conditions=[[str(v) for v in row] for row in matrix.tolist()],
                determinant=str(matrix.det()),classically_closed_cartan_image_dimension=2-matrix.rank())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    assert __debug__
    data=run()
    m.write(args.out,data)
    print(data["cartan_conditions"],data["determinant"],data["classically_closed_cartan_image_dimension"])
