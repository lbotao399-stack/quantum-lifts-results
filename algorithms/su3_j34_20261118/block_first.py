"""Restrict the *parent* -Q1 O34 to traceless 2x2 matrices and extract abelian coefficients over Q."""
import argparse, gzip, hashlib, json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path

HERE=Path(__file__).resolve().parent
RHS_SHA='9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'
# Integer keys in the program only: three low bits encode theta; upper bits store (z+,z-).
def code(letter):
    a,b,s=letter; assert 0<=a<=1 and 0<=b<=1
    return (a*2+b)*8+sum(1<<i for i in s)
def letter(x):return (x//16,(x//8)%2,tuple(i for i in range(3) if x&(1<<i)))
def odd(x):return (1+(x&7).bit_count())%2
# Names are for coefficient extraction; pass original full-letter triples in the published output.
PHI1=code((0,0,(0,))); CHANNELS=(code((0,1,())), code((0,0,(2,))))
@lru_cache(None)
def parity(t):return sum(map(odd,t))%2
@lru_cache(None)
def canonical_trace(t):
    """Signed cyclic/reversal canonical form for a traceless 2x2 trace."""
    if len(t)<2:return 0,()
    f=sum(map(odd,t)); signs={}
    for seq,start in ((t,1),(t[::-1],(-1)**(len(t)+f*(f-1)//2))):
        for i,x in enumerate(seq):
            rotated=seq[i:]+seq[:i]
            if signs.setdefault(rotated,start)!=start:return 0,()
            if odd(x) and (f-1)%2:start=-start
    best=min(signs);return signs[best],best
@lru_cache(131072)
def canonical_multi(traces):
    fs=[t for t in traces if parity(t)]
    if len(fs)!=len(set(fs)):return 0,()
    sign=(-1)**sum(t>u for i,t in enumerate(fs) for u in fs[i+1:]);return sign,tuple(sorted(traces))
@lru_cache(None)
def trace_options(trace, other):
    """Set a trace letter to its upper traceless block or to its abelian diagonal. Move the extracted odd coordinate to the left."""
    scal=[i for i,x in enumerate(trace) if x==PHI1]; distinguished=[i for i,x in enumerate(trace) if x==other]
    rows=Counter()
    for n in range(min(3,len(scal))+1):
        for picks in combinations(scal,n):
            for j in [None]+distinguished:
                rem=tuple(x for i,x in enumerate(trace) if i not in picks and i!=j)
                sign=(-1)**sum(odd(x) for x in trace[:j]) if j is not None and odd(other) else 1
                if rem: s,w=canonical_trace(rem)
                else: s,w=2+(-2)**len(trace),()
                if s:rows[n,int(j is not None),w]+=s*sign
    return tuple((n,k,w,c) for (n,k,w),c in rows.items() if c)
def extract_word(word,other):
    if sum(t.count(PHI1) for t in word)<3 or not any(other in t for t in word):return ()
    outs=Counter(); opts=tuple(trace_options(t,other) for t in word); ps=tuple(parity(t) for t in word); L=len(word)
    def visit(i,n,k,par,c,factors):
        if i==L:
            if (n,k)==(3,1):
                s,w=canonical_multi(factors)
                if s:outs[w]+=s*c
            return
        for a,b,t,v in opts[i]:
            if n+a>3 or k+b>1 or (i==L-1 and (n+a,k+b)!=(3,1)):continue
            visit(i+1,n+a,k+b,par^ps[i],(-v if b and par and odd(other) else v)*c, factors+((t,) if t else ()))
    visit(0,0,0,0,1,());return outs.items()
def encoded_stream(path):
    with gzip.open(path,'rb') as f:
        for line in f:
            c,w=json.loads(line); frac=Fraction(c)*120; assert frac.denominator==1
            yield line,int(frac),tuple(tuple(map(code,t)) for t in w)
def dump_rows(path,rows):
    raw=b''.join((json.dumps([c,w],separators=(',',':'))+'\n').encode() for w,c in sorted(rows.items()) if c)
    with path.open('wb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:z.write(raw)
    return hashlib.sha256(raw).hexdigest()
def extract(root,out):
    out.mkdir(parents=True,exist_ok=True); h=hashlib.sha256(); rr=[defaultdict(int) for _ in CHANNELS];n=0
    for line,c,w in encoded_stream(root/'negative_Q1.jsonl.gz'):
        h.update(line);n+=1
        for other,r in zip(CHANNELS,rr):
            for k,v in extract_word(w,other):r[k]+=v*c
        if n%100000==0:print('parent rows',n,flush=True)
    assert h.hexdigest()==RHS_SHA and n==517468
    result={'parent_rhs_sha256':RHS_SHA,'parent_rows':n,'abelian_colour_diagonal':[1,1,-2],'common_denominator':120,'channels':[]}
    for tag,other,r in zip(('lambda_minus','phi3'),CHANNELS,rr):
        reduced={w:c for w,c in r.items() if c};assert all(sum(map(len,w))==7 for w in reduced)
        mu=(1,0,4,4,3) if tag=='lambda_minus' else (1,1,4,4,2)
        for w in reduced:
            ls=[letter(x) for t in w for x in t];assert (sum(x[0] for x in ls),sum(x[1] for x in ls),*(sum(i in x[2] for x in ls) for i in range(3)))==mu
        digest=dump_rows(out/(tag+'.jsonl.gz'),reduced);result['channels'].append(dict(name=tag,abelian_letters=[letter(PHI1)]*3+[letter(other)],mu=mu,length=7,formal_terms=len(reduced),trace_sha256=digest))
    (out/'extraction.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE);p.add_argument('--out',type=Path,default=HERE/'block_first');a=p.parse_args();extract(a.root,a.out)
