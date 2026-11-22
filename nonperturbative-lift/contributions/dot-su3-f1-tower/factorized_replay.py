"""Expand the level-32 SU(3) operator from the published arXiv:2412.08695v2 TeX."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations, product
import json
from pathlib import Path
import re
import sys


def additional_preimages(dmk):
    # App. A.1. Only the fundamental [1,0] and [0,2] tensors are missing
    # from the level-30 implementation. Their two lower labels are symmetrized.
    def tr(*letters):
        return dmk.tr(*({i: 1} for i in letters))

    @lru_cache(None)
    def vector(i):
        terms = []
        for a, b in product(range(3), repeat=2):
            A, B = 3+a, 3+b
            terms.extend([
                (3, tr(i,A,a,b,B)),
                (1, dmk.mul(tr(i,a),tr(b,A,B))), (1, dmk.mul(tr(i,a),tr(b,B,A))),
                (-3, dmk.mul(tr(i,A),tr(a,b,B))), (3, dmk.mul(tr(i,A),tr(b,a,B))),
                (-1, dmk.mul(tr(i,A,B),tr(a,b)))])
        return dmk.add(*terms)

    @lru_cache(None)
    def tensor(i, j):
        terms = []
        for left, right in ((i,j), (j,i)):
            J = 3+right
            for a,b,c in permutations(range(3)):
                e = Fraction(5*dmk.eps(a,b,c), 18)
                terms.extend([(2*e, tr(3+left,J,a,b,c)), (e, tr(3+left,a,J,b,c))])
            for a,b in permutations(range(3),2):
                e = Fraction(dmk.eps(a,b,left), 2)
                if not e:
                    continue
                for c in range(3):
                    C = 3+c
                    for scalars in ((a,b,c), (a,c,b), (c,a,b)):
                        for k in range(4):
                            terms.append((e, tr(J,*scalars[:k],C,*scalars[k:])))
                    terms.extend([
                        (-5*e/3, dmk.mul(tr(J,a,b),tr(C,c))),
                        (-e/3, dmk.mul(tr(J,a,c),tr(C,b))),
                        (-e/3, dmk.mul(tr(J,c,a),tr(C,b))),
                        (e/3, dmk.mul(tr(J,b),tr(C,a,c))),
                        (e/3, dmk.mul(tr(J,b),tr(C,c,a)))])
        return dmk.add(*terms)
    return {"r14_10":vector, "r14_02":tensor}


def build(dmk, terms):
    dmk._RAW.update(additional_preimages(dmk))
    total = Counter()
    for k, (coefficient, factors) in enumerate(terms, 1):
        # This printed operator includes the signs of its odd prefactors.
        # Applying Q_0 restores those signs in the relation among the R tensors.
        expression = dmk.term_expr(factors, dict.fromkeys("ijkl", 0))
        for word, number in expression.items():
            total[word] += coefficient * number
        print(f"Expression {k}/47, {len(expression)} terms", flush=True)
    return {w:c for w,c in total.items() if c}



from os import environ
from pathlib import Path
R007=Path(environ['R007_EXTRACT']).resolve()
for rel,expected in [
  ('q0-certificate/vendor/fortuity/bmn/dmk.py','535ec58677878ffe5b613ec93bc4539c7e58eeefc8e532b98e791ccff0fc86b8'),
  ('q0-certificate/vendor/fortuity/bmn/words.py','f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e'),
  ('q0-certificate/q0/finite_colour.py','4d5d70b7ead3676a35bef4417686ba378e0b51829ead1d880ef943dd9ea42f6b'),
]:
 assert hashlib.sha256((R007/rel).read_bytes()).hexdigest()==expected
sys.path.insert(0,str(R007/'q0-certificate/vendor'))
sys.path.insert(0,str(R007/'q0-certificate/q0'))
assert len(sys.argv)==2, 'supply the frozen O32_C_even.json.gz'
from fortuity.bmn import dmk,words
from finite_colour import evaluate,trace
R3=range(3)
def R10(i):
 return dmk.add(*[(dmk.eps(a,b,c),dmk.mul(dmk.u2(i,a),dmk.u2(d,b),dmk.v2(c,d))) for a,b,c in permutations(R3) for d in R3])
def R02(i,j):
 terms=[]
 for I,J in ((i,j),(j,i)):
  for a,b in product(R3,repeat=2):
   e=Fraction(dmk.eps(a,b,I),2)
   if not e: continue
   for c,d,f in permutations(R3):
    terms.append((e*dmk.eps(c,d,f),dmk.mul(dmk.u2(a,c),dmk.u2(b,d),dmk.v2(f,J))))
   for c,d,f in product(R3,repeat=3):
    k=dmk.eps(J,c,d)
    if k: terms.append((-2*e*k,dmk.mul(dmk.u3(a,c,f),dmk.v3(b,d,f))))
 return dmk.add(*terms)

terms=[(c,[(name,' '.join(ids)) for name,ids in fs]) for c,fs in json.loads(Path(__file__).with_name('factors.json').read_text())]
assert len(terms)==47
print(Counter(f[-1][0] for c,f in terms),flush=True)
for name,I in [('r12_03',(2,2,2)),('r12_11',(0,2)),('r12_22',(0,0,2,2)),('r14_13',(0,2,2,2)),('r16_23',(0,0,2,2,2))]:
 r=dmk.tensor(name)(*I); R=dmk.tensor(dmk.R_OF_r[name])(*I)
 assert not evaluate(dmk.add((1,words.Q_expr(r)),(-1,R))),name
 print('existing projected identity',name,I,'PASS',flush=True)
 trace.cache_clear()
def even_part(expr, words):
    result = Counter()
    for multi, coefficient in expr.items():
        for reverse in (False, True):
            sign, traces = 1, []
            for trace in multi:
                if reverse:
                    odd = sum(map(words.is_ferm, trace))
                    sign *= (-1) ** (len(trace) + odd * (odd - 1) // 2)
                cyclic, word = words.canon_word(trace[::-1] if reverse else trace)
                sign *= cyclic
                traces.append(word)
            assert sign
            exchange, word = words.canon_multi(traces)
            result[word] += coefficient * sign * exchange / 2
    return {word: value for word, value in result.items() if value}


rr=additional_preimages(dmk)
for family,ind,R in [('r14_10',(0,),R10),('r14_02',(2,2),R02),('r14_02',(1,2),R02)]:
    r=rr[family](*ind)
    assert not evaluate(dmk.add((1,words.Q_expr(r)),(-1,R(*ind))))
    trace.cache_clear()
    print('exact new r/R',family,ind,'PASS',flush=True)

from contextlib import redirect_stdout
from io import StringIO
with redirect_stdout(StringIO()):
 polynomial=even_part(build(dmk,terms),words)
raw_source=Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(raw_source).hexdigest()=='40c4cbb876c01811d0dd83035ab9080ad71ebb603331304b309f4d0d0335e623'
source=gzip.decompress(raw_source)
new=json.dumps(dict(level=32,weight=[8,4,4],colour_parity='even',length=11,n_terms=len(polynomial),terms=[[str(c),w] for w,c in sorted(polynomial.items())]),separators=(',',':')).encode()+b'\n'
assert new==source, (len(new),len(source),len(polynomial))
print('FROZEN SOURCE BYTE MATCH',len(polynomial),hashlib.sha256(new).hexdigest(),flush=True)

def mult(*es):
 out={():Fraction(1)}
 for e in es:
  o=Counter()
  for x,c in out.items():
   w=[i for i in x if i[0] in ('V2','V3')]
   for y,d in e.items():
    z=[i for i in y if i[0] in ('V2','V3')]
    if set(w)&set(z): continue
    o[tuple(sorted(x+y))]+=(-1)**sum(i>j for i in w for j in z)*c*d
  out={k:v for k,v in o.items() if v}
 return out
def gen(name,*i):return {((name,*i),):1}
def u2(i,j):return gen('U2',*sorted((i,j)))
def u3(i,j,k):return gen('U3',*sorted((i,j,k)))
def raw2(i,j):return gen('V2',i,j)
def raw3(i,j,k):return gen('V3',*sorted((i,j)),k)
def v2(i,j): return dmk.add((1,raw2(i,j)),*([(Fraction(-1,3),raw2(a,a)) for a in R3] if i==j else []))
def v3(i,j,k):return dmk.add((1,raw3(i,j,k)),*([(Fraction(-1,4),raw3(j,a,a)) for a in R3] if i==k else []),*([(Fraction(-1,4),raw3(i,a,a)) for a in R3] if j==k else []))
for f in [u2,u3,v2,v3,mult]: setattr(dmk,f.__name__ if f!=mult else 'mul',f)
dmk._RAW.update(u2=u2,u3=u3,v2=v2,v3=v3,R14_10=R10,R14_02=R02)
dmk.R_OF_r.update(r14_10='R14_10',r14_02='R14_02')
dmk.tensor.cache_clear()
for name in dmk.R_OF_r.values():
 fn=dmk._RAW[name]
 if hasattr(fn,"cache_clear"): fn.cache_clear()
rr=[];raw=[]
for c,f in terms:
 p=dmk.term_expr(f,dict.fromkeys('ijkl',0),use_R=True)
 rr.append((c*(-1)**dmk.prefactor_parity(f),p));raw.append((c,p))
assert len(dmk.add(*rr))==0 and len(dmk.add(*raw))==716
print("graded",len(dmk.add(*rr)),"ungraded",len(dmk.add(*raw)),flush=True)
