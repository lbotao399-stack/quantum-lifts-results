"""Rebuild the five-scalar correction on the SU(2) x U(1) tests at level 32."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from flint import fmpz_mat
from su3_level32_layers import DATA, content, formal, images, joint, residual
from su3_level32_trace_factors import code, evaluator, multiply
from even_content_frames import projected
from verify import add_scaled


BASIS = ((1,0,0,0,0,0,0,0,-1), (0,0,0,0,1,0,0,0,-1),
         (0,1,0,0,0,0,0,0,0), (0,0,0,1,0,0,0,0,0))
SCALAR = (1,0,0,0,1,0,0,0,-2)
TRANSLATION = bytes(map(code,formal.ALPH)) + bytes(256-len(formal.ALPH))


def projected_correction(words):
    even = {}
    for w,c in words.items():
        add_scaled(even,projected(w),c)
    result = {}
    for variant,sign in zip(images(even),(1,-1,-1,1),strict=True):
        add_scaled(result,variant,sign)
    return result


def extensions():
    result = {}
    for parent,splits in enumerate(formal.SPLITS):
        for sign,pair in splits:
            if pair[0] == 0:
                assert pair[1] and sign in (-1,1)
                assert (-sign,bytes((pair[1],0))) in splits and pair[1] not in result
                result[pair[1]] = (parent,sign)
    return result


@lru_cache(maxsize=10000)
def insert(trace):
    parents = extensions()
    terms = Counter()
    phase = 1
    for i,letter in enumerate(trace):
        if letter in parents:
            parent,sign = parents[letter]
            rotate,word = formal.canon(trace[:i]+bytes((parent,))+trace[i+1:])
            terms[word] += phase*sign*rotate
        if formal.ODD[letter]:
            phase = -phase
    return {w:c for w,c in terms.items() if c}


def expand(recipe):
    word = joint.s.unpack_hex(recipe["word"])
    if "insert" not in recipe:
        return projected_correction({word:1})
    i,j = recipe["insert"]
    assert 0 <= i < j < len(word)
    terms = Counter()
    for left,x in insert(word[i]).items():
        for right,y in insert(word[j]).items():
            factors = list(word)
            factors[i],factors[j] = left,right
            sign,canon = formal.multi(factors)
            terms[canon] += sign*x*y
    return projected_correction({w:c for w,c in terms.items() if c})


def variation(words):
    result = {}
    for word,c in words.items():
        add_scaled(result,formal.q(word),c)
    return result


def test_point(target, seed):
    evaluate, assignment = evaluator(target.translate(TRANSLATION),BASIS,seed,{TRANSLATION[0]:SCALAR})
    return lambda word:evaluate(tuple(t.translate(TRANSLATION) for t in word)), assignment


def evaluate_row(row, evaluate):
    return sum((v*evaluate(w) for w,v in row.items()),Fraction())


def verify(data, folder, out, alternate_seed):
    assert all(multiply(SCALAR,x)==multiply(x,SCALAR) for x in BASIS)
    probes = sorted(joint.s.read(folder/"gl2_probes.json"),key=lambda p:p["content"])
    saved = joint.s.read(folder/"levi_solved.json")
    targets = [bytes.fromhex(p["content"]) for p in probes]
    assert [v.hex() for v in targets] == saved["targets"] and len(set(targets)) == len(probes)
    positions = {v:i for i,v in enumerate(targets)}
    tests = []
    for t,p in zip(targets,probes,strict=True):
        ev,assigned = test_point(t,p["seed"])
        assert assigned["exterior_degree"] == p["assignment"]["exterior_degree"]
        assert assigned["assignments"] == [[k,i,tuple(m)] for k,i,m in p["assignment"]["assignments"]]
        tests.append(ev)
    old = {joint.s.unpack_hex(w):Fraction(c) for w,c in joint.s.read(folder/"correction.json.gz")}
    assert all(content(w).count(0)==5 for w in old)
    before = residual(data,old,lambda t:t.count(0)==5)
    assert {t for t,r in before.items() if r} == set(targets)
    rhs = [-evaluate_row(before[t],ev) for t,ev in zip(targets,tests,strict=True)]
    assert rhs == list(map(Fraction,saved["rhs"])) == [-Fraction(p["coefficient"]) for p in probes]
    delta = {}
    rows = []
    for recipe,weight,archived in zip(saved["recipes"],saved["coefficients"],saved["values"],strict=True):
        terms = expand(recipe)
        assert all(content(w).count(0)==5 for w in terms)
        image = variation(terms)
        assert all(content(w).count(0)==5 for w in image)
        values = [0]*len(targets)
        for w,c in image.items():
            t = content(w)
            if t in positions:
                j = positions[t]
                values[j] += c*tests[j](w)
        assert values == archived
        rows.append(values)
        add_scaled(delta,terms,Fraction(weight))
    pivots = saved["pivot_columns"]
    assert len(set(pivots))==len(rows)==len(pivots)
    minor = fmpz_mat([[r[i] for i in pivots] for r in rows]).det()
    assert minor
    improved = old.copy()
    add_scaled(improved,delta,1)
    for w,c in variation(delta).items():
        add_scaled(before[content(w)],{w:1},c)
    assert all(evaluate_row(before[t],ev)==0 for t,ev in zip(targets,tests,strict=True))
    holdout = {}
    for t,row in sorted(before.items()):
        if row:
            ev,_ = test_point(t,alternate_seed)
            value = evaluate_row(row,ev)
            if value:
                holdout[t.hex()] = str(value)
    out.mkdir(parents=True,exist_ok=True)
    joint.s.write(out/"levi_correction.json.gz",[[joint.s.hexes(w),str(c)] for w,c in sorted(improved.items())])
    result = dict(base_terms=len(old),change_terms=len(delta),combined_terms=len(improved),tests=len(tests),
                  image_rank=len(rows),determinant=str(minor),nonzero_fit_residuals=0,
                  five_scalar_contents=len(before),five_scalar_formal_zero=sum(not r for r in before.values()),
                  independent_seed=alternate_seed,independent_nonzero=holdout)
    joint.s.write(out/"levi_replay.json",result)
    print({k:v for k,v in result.items() if k not in ("determinant","independent_nonzero")},
          "nonzero independent tests",len(holdout),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--folder",type=Path,default=DATA/"scalar_layers")
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--alternate-seed",type=int,default=620003)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    verify(args.data,args.folder,args.out,args.alternate_seed)
