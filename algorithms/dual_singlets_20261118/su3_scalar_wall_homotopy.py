"""Check the scalar split on all off-diagonal entries and evaluate a trial wall correction transversely."""
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path

from su3_level32_levi import SCALAR, TRANSLATION, extensions
from su3_level32_trace_factors import evaluator
from su3_level32_scalar_descent import restored
from su3_levi_equation import DATA, allowed_targets, first_residual, formal, joint


def derivation(word):
    result = Counter()
    sign = 1
    for i, letter in enumerate(word):
        for c, two in formal.SPLITS[letter]:
            result[word[:i]+two+word[i+1:]] += sign*c
        if formal.ODD[letter]:
            sign = -sign
    return result


def check_pairs():
    pairs = extensions()
    marked = {m for m, _ in pairs.values()}
    assert not (set(pairs) & marked) and set(pairs) | marked == set(range(1,len(formal.ALPH)))
    assert not formal.SPLITS[0]
    for k, splits in enumerate(formal.SPLITS):
        squared = Counter()
        for c, word in splits:
            for v, n in derivation(word).items():
                squared[v] += c*n
        assert not any(squared.values())
        projected = [(c,w) for c,w in splits if 0 in w]
        if k in marked:
            bare, sign = next((a,s) for a,(b,s) in pairs.items() if b==k)
            assert Counter(projected) == Counter([(sign,bytes((0,bare))),(-sign,bytes((bare,0)))])
            assert formal.ODD[bare] != formal.ODD[k]
        else:
            assert not projected
    diagonal = (1,1,-2)
    contractions = []
    for i in range(3):
        for j in range(3):
            gap = diagonal[i]-diagonal[j]
            if gap:
                for bare, (m, sign) in sorted(pairs.items()):
                    coefficient = sign*gap
                    inverse = Fraction(1,coefficient)
                    assert coefficient*inverse == 1
                    contractions.append([i,j,bare,m,coefficient,str(inverse)])
    assert len(contractions) == 60
    assert all(sum(degrees[:2])+sum(degrees[3:]) >= 1 for degrees in formal.DEGREES[1:])
    return dict(classical_letters=len(formal.ALPH), ordered_splits=sum(map(len,formal.SPLITS)),
                matrix_entries=contractions)


def normal_value(t, row, parameter, seed):
    diagonal = [(1,0,0,0,0,0,0,0,-1),(0,0,0,0,1,0,0,0,-1)]
    mixed = {2,5,6,7}
    basis = diagonal + [tuple((parameter if i in mixed else 1)*int(i==j) for j in range(9))
                        for i in (1,2,3,5,6,7)]
    evaluate, assignment = evaluator(t.translate(TRANSLATION), basis, seed, {TRANSLATION[0]: SCALAR})
    value = sum((c*evaluate(tuple(tr.translate(TRANSLATION) for tr in w)) for w,c in row.items()), Fraction())
    return value, assignment


def interpolate(values):
    coefficients = [Fraction()] * len(values)
    for i, value in enumerate(values):
        polynomial = [Fraction(value)]
        for j in range(len(values)):
            if j == i:
                continue
            next_poly = [Fraction()] * (len(polynomial)+1)
            for k,c in enumerate(polynomial):
                next_poly[k] -= c*j*j/(i*i-j*j)
                next_poly[k+1] += c/(i*i-j*j)
            polynomial = next_poly
        coefficients = [a+b for a,b in zip(coefficients,polynomial,strict=True)]
    return coefficients


def check(level, data, correction, seed):
    pairs = check_pairs()
    targets = allowed_targets(level)
    limit = max(len(t)-t.count(0) for t in targets)
    vector = restored(joint.s.read(correction))
    rows = first_residual(level, data, vector, lambda _: True)
    for t,row in sorted(rows.items(),key=lambda item:(len(item[1]),item[0])):
        if not row:
            continue
        at_one, assignment = normal_value(t, row, 1, seed)
        if not at_one:
            continue
        degree = (len(t)-t.count(0))//2
        values = [normal_value(t, row, x, seed)[0] for x in range(degree+1)]
        assert values[0] == 0 and values[1] == at_one
        coefficients = interpolate(values)
        extra, _ = normal_value(t, row, -degree-1, seed)
        assert extra == sum(c * (degree+1)**(2*i) for i,c in enumerate(coefficients))
        if not coefficients[1]:
            continue
        return dict(level=level, scalar_eigenvalues=[1,1,-2], maximum_non_scalar_letters=limit,
                    minimum_scalar_power_in_primitive=10-limit, scalar_split=pairs,
                    trial_normal_witness=dict(content=t.hex(), seed=seed, assignment_at_one=assignment,
                        coefficients_in_t_squared=list(map(str,coefficients)), value_at_one=str(at_one),
                        extra_point=-degree-1, extra_value=str(extra)))
    raise ValueError("no quadratic normal coefficient found")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level",type=int,required=True,choices=(32,34))
    p.add_argument("--data",type=Path)
    p.add_argument("--correction",type=Path,required=True)
    p.add_argument("--seed",type=int,default=894611)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    data=args.data or (DATA if args.level==32 else Path(__file__).resolve().parents[1]/"su3_j34_20261118")
    report=check(args.level,data,args.correction,args.seed)
    joint.s.write(args.out,report)
    print(args.level,report["trial_normal_witness"]["coefficients_in_t_squared"],flush=True)
