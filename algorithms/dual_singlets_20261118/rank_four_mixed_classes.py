"""Rational singlet representatives and the two terminal projectors on the Cartan kernel."""
import argparse
from math import gcd
from pathlib import Path

from flint import fmpq_mat, fmpz_mat

import mixed_incoming as m
from mixed33_exact import read
from orthogonal_relations import lower
from rank_four_terminal import bounds, quotient
from terminal_commutators import words


def primitive(vector):
    values = list(map(int, vector))
    div = gcd(*values)
    return [v//div for v in values]


def run(sp, so, functional, output):
    sym, orth, test = map(read,(sp,so,functional))
    case = bounds(3)
    weight = tuple(case["weight"])
    basis, cov = quotient(weight,case["formal"])
    # The highest-weight kernel realizes the semisimple singlet quotient.
    raising = []
    for src,dst in ((1,0),(3,2),(4,3)):
        target = list(weight)
        target[src] -= 1
        target[dst] += 1
        target_words = words(tuple(target))
        column = [lower(w,src,dst) for w in basis]
        raising.extend([[c.get(w,0) for c in column] for w in target_words])
    matrix = fmpz_mat(raising)
    hw, n = matrix.nullspace()
    assert n == 10
    highest = fmpq_mat([primitive([hw[i,j] for i in range(len(basis))]) for j in range(n)])
    assert not any((fmpq_mat(matrix)*highest.transpose()).entries())
    transition = highest*fmpq_mat(cov)
    assert transition.rank() == n
    gram = [r["singlet"] for r in orth["relations"]]
    chosen, complement = [], []
    for i,rel in enumerate(sym["relations"]):
        v = rel["singlet"]
        if fmpz_mat(gram+complement+[v]).rank() > len(gram)+len(complement):
            chosen.append(i)
            complement.append(v)
    assert len(gram) == 5 and len(complement) == 2
    polynomials = fmpq_mat(complement)*transition.inv()*highest
    assert polynomials*fmpq_mat(cov) == fmpq_mat(complement)
    assert not any((fmpq_mat(matrix)*polynomials.transpose()).entries())
    coeff = dict(test["weight_coefficients"])
    poles = []
    for v,index in zip(polynomials.tolist(),chosen):
        blocks = [(i,str(c)) for i,c in enumerate(v) if c]
        value = sum(c*coeff.get(i,0) for i,c in enumerate(v))
        assert value == test["symplectic_relations"][index]
        poles.append(dict(symplectic_relation=index,singlet=blocks,pure_substitution=str(value)))
    # Solve for duals vanishing on every SO(9) trace relation.
    frame = fmpq_mat(gram+complement)
    transpose = frame.transpose().tolist()
    rows = []
    for i,row in enumerate(transpose):
        if fmpq_mat([transpose[j] for j in rows]+[row]).rank() > len(rows):
            rows.append(i)
    rhs = fmpq_mat([[0,0]]*len(gram)+[[1,0],[0,1]])
    sol = fmpq_mat([[frame[i,j] for j in rows] for i in range(frame.nrows())]).solve(rhs)
    extended = fmpq_mat(frame.ncols(),2)
    for i,j in enumerate(rows):
        for k in range(2):
            extended[j,k] = sol[i,k]
    assert frame*extended == rhs
    dual = fmpq_mat(cov)*extended
    assert polynomials*dual == fmpq_mat([[1,0],[0,1]])
    m.write(output,dict(weight=weight,formal_weight_basis=basis,orthogonal_relations=gram,
                        quotient_representatives=complement,classes=poles,
                        weight_functionals=[[(i,str(dual[i,k])) for i in range(dual.nrows()) if dual[i,k]]
                                            for k in range(2)],
                        projector_values=[[1,0],[0,1]],classical_SO=5,classical_Sp=3))
    print("two classes",[(len(p["singlet"]),p["pure_substitution"]) for p in poles],flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("sp","so","functional","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    assert __debug__
    run(**vars(p.parse_args()))
