"""Evaluate restricted full-letter traces as exact matrix-entry polynomials, and solve the first equation over Q."""
import argparse, gzip, hashlib, json, sys
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from block_first import HERE, RHS_SHA, canonical_multi, canonical_trace, code, letter, odd
from verify import add_scaled, exterior_sign, multiply

@lru_cache(None)
def trace(t):
    result=Counter(); assert len(t)<16
    for path in product((0,1),repeat=len(t)):
        even=bits=0;sign=1
        for x,i,j in zip(t,path,path[1:]+path[:1]):
            entry=0 if i==j else 1+i; slot=3*x+entry
            if i==j==1:sign=-sign
            if odd(x):
                bit=1<<slot;sign*=exterior_sign(bits,bit);bits|=bit
                if not sign:break
            else:even+=1<<(4*slot)
        if sign:result[even,bits]+=sign
    return {m:c for m,c in result.items() if c}
@lru_cache(512)
def value(w):
    assert sum(map(len,w))<16
    r={(0,0):1}
    for t in w:r=multiply(r,trace(t))
    return r

def polynomial(expr):
    r={}
    for w,c in expr.items():add_scaled(r,value(w),c)
    return r

def words(raw):return tuple(tuple(t) for t in raw)
def to_fort(w):return tuple(tuple(map(letter,t)) for t in w)
def from_fort(w):return tuple(tuple(map(code,t)) for t in w)
def q(w,fw):return {from_fort(m):c for m,c in fw.Q_multi(to_fort(w)).items() if c}
def all_q(expr,fw):
    r={}
    for w,c in expr.items():add_scaled(r,q(w,fw),c)
    return r

def canonical(w):
    sign=1;parts=[]
    for t in w:
        s,p=canonical_trace(t)
        if not s:return 0,()
        sign*=s;parts.append(p)
    s,p=canonical_multi(tuple(parts));return s*sign,p

def solve(target,preds,fw):
    # Each stored row is the Q0 image of the same exact combination of predecessor trace words.
    basis={};used_count=0
    def reduce(g,comb):
        while g:
            k=min(g)
            if k not in basis:return k
            u,v=basis[k];c=g[k];add_scaled(g,u,-c);add_scaled(comb,v,-c)
        return None
    for i,w in enumerate(preds):
        g=polynomial(q(w,fw)); v={i:1};k=reduce(g,v)
        if k is not None:
            s=Fraction(1,g[k]);basis[k]=({m:s*c for m,c in g.items()},{j:s*c for j,c in v.items()})
    # For a valid normal form, also reduce later pivots after the first free coordinate.
    rem={};v={};g=target.copy()
    while g:
        k=min(g)
        if k in basis:
            u,uv=basis[k];c=g[k];add_scaled(g,u,-c);add_scaled(v,uv,c)
        else:rem[k]=g.pop(k)
    pre={preds[i]:c for i,c in v.items()};assert not any(k in basis for k in rem)
    check=polynomial(all_q(pre,fw));add_scaled(check,rem,1);assert check==target
    return len(basis),rem,pre,basis

def run(root,runtime,mode):
    sys.path.insert(0,str(runtime.resolve()));import sp_environment
    from fortuity.full import words as fw
    from fortuity.full.letters import Q_letter, letters_up_to, is_ferm
    from fortuity.full.wordenum import multitraces_with_content
    from g2_full_classical_screen_v1 import contents
    # Each ordered quadratic term pairs with its graded reverse to make a colour supercommutator.
    # Central abelian coordinates then drop from the classical action under block restriction.
    letters=letters_up_to(1,1);pairs=0
    for x in letters:
        D=Counter()
        for c,y,z in Q_letter(x):D[y,z]+=c
        assert all(D[y,z]==-(-1)**(is_ferm(y)*is_ferm(z))*D[z,y] for y,z in D)
        assert is_ferm(x)==bool(odd(code(x)))
        pairs+=len(D)
    base=root/'block_first';manifest=json.loads((base/'extraction.json').read_bytes());assert manifest['parent_rhs_sha256']==RHS_SHA and manifest['common_denominator']==120
    results=[]
    for ch in manifest['channels']:
        raw=gzip.decompress((base/(ch['name']+'.jsonl.gz')).read_bytes());assert hashlib.sha256(raw).hexdigest()==ch['trace_sha256'];expr={words(w):Fraction(c,120) for c,w in map(json.loads,raw.splitlines())};assert len(expr)==ch['formal_terms'];rhs=polynomial(expr);qr=polynomial(all_q(expr,fw));print(ch['name'],len(expr),'matrix monomials',len(rhs),'closure residual',len(qr),flush=True);assert not qr
        mu=tuple(ch['mu']);cs=contents(mu,6);allpreds=[from_fort(w) for cont in cs for w in multitraces_with_content(cont,3)];assert len(set(allpreds))==len(allpreds);preds=sorted({p for w in allpreds for s,p in [canonical(w)] if s});r=dict(name=ch['name'],degree=mu,trace_terms=len(expr),matrix_monomials=len(rhs),closure_residual=len(qr),predecessor_contents=len(cs),predecessor_words=len(allpreds),signed_distinct_predecessor_words=len(preds));print('predecessors',len(allpreds),len(preds),flush=True)
        if not rhs:
            r.update(boundary=True,primitive=[],residual_monomials=0)
        elif mode=='solve':
            rank,res,pre,basis=solve(rhs,preds,fw);r.update(boundary=not res,boundary_rank=rank,primitive=[[str(c),w] for w,c in sorted(pre.items())],residual_monomials=len(res))
            if res:
                k=min(res);ell={k:1}
                for m in sorted(basis,reverse=True):
                    row,_=basis[m];c=sum(ell.get(i,0)*a for i,a in row.items());
                    if c:add_scaled(ell,{m:1},-c)
                def pairing(x):return sum(v*x.get(m,0) for m,v in ell.items())
                assert pairing(rhs)==res[k] and all(pairing(polynomial(q(w,fw)))==0 for w in preds)
                r['witness']={'entries':[[str(m[0]),str(m[1]),str(a)] for m,a in sorted(ell.items())],'pairing':str(res[k])}
        results.append(r)
    report=dict(classical_rule_letters=len(letters),ordered_quadratic_coefficients=pairs,channels=results)
    (base/('reduction.json' if mode=='solve' else 'polynomial_check.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps([{k:v for k,v in r.items() if k not in ('primitive','witness')} for r in results]),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE);p.add_argument('--runtime',type=Path,default=HERE.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync');p.add_argument('--mode',choices=('polynomial','solve'),default='solve');a=p.parse_args();run(a.root,a.runtime,a.mode)
