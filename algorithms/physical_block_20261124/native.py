"""Exact native traces and simultaneous highest-weight equations.

No basis construction by averaging or by an extremal operator is used.
The complete word enumeration and Q letter rule are pinned repository inputs.
All linear equations in this module are solved over the rational numbers.
"""
from pathlib import Path
import sys, math, gzip, json, hashlib
from functools import lru_cache
from itertools import product
from collections import Counter, defaultdict

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RUNTIME=ROOT/'algorithms/su3_second_correction_20261118/runtime/work'
sys.path.insert(0,str(RUNTIME/'cluster-next-v2.nosync/work/so8-bmn/vendor'))
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm
from fortuity.full.wordenum import multitraces_with_content
from flint import fmpz_mat

def encode(w): return [[[a,b,list(s)] for a,b,s in tr] for tr in w]
def decode(w): return tuple(tuple((a,b,tuple(s)) for a,b,s in tr) for tr in w)
def nferm(w): return sum(is_ferm(l) for tr in w for l in tr)
def compact(obj): return (json.dumps(obj,separators=(',',':'))+'\n').encode()
def write(path,obj):
    data=compact(obj)
    path=Path(path)
    path.write_bytes(gzip.compress(data,mtime=0) if path.suffix=='.gz' else data)
def read(path):
    data=Path(path).read_bytes()
    return json.loads(gzip.decompress(data) if str(path).endswith('.gz') else data)
def export_rows(rows): return [[[int(n),encode(w)] for w,n in sorted(r.items())] for r in rows]
def import_rows(rows): return [{decode(w):int(n) for n,w in r} for r in rows]

@lru_cache(None)
def contents(mu,length):
    if min(mu)<0 or length<0 or sum(mu)<length: return ()
    alphabet=[]
    for a,b in product(range(mu[0]+1),range(mu[1]+1)):
        for bits in product((0,1),repeat=3):
            wt=(a,b)+bits
            if 1<=sum(wt)<=sum(mu)-length+1 and all(x<=y for x,y in zip(wt,mu)):
                alphabet.append((wt,(a,b,tuple(i for i,v in enumerate(bits) if v))))
    alphabet.sort(reverse=True)
    answer=[]
    def visit(start,remain,left,word):
        if not left:
            if not any(remain): answer.append(tuple(sorted(word)))
            return
        for j in range(start,len(alphabet)):
            wt,letter=alphabet[j]
            maximum=min([left]+[remain[k]//wt[k] for k in range(5) if wt[k]])
            for count in range(1,maximum+1):
                visit(j+1,tuple(remain[k]-count*wt[k] for k in range(5)),left-count,word+(letter,)*count)
    visit(0,mu,length,())
    return tuple(answer)

def canonical(multi):
    traces=[];total=1
    for word in multi:
        sig,value=fw.canon_word(word)
        reverse_sig,reverse=fw.canon_word(word[::-1])
        odd=sum(is_ferm(l) for l in word)
        reverse_sig*=(-1)**(len(word)+odd*(odd-1)//2)
        if not sig or (value==reverse and sig!=reverse_sig): return 0,None
        if reverse<value: value,sig=reverse,reverse_sig
        total*=sig;traces.append(value)
    sig,value=fw.canon_multi(traces)
    return total*sig,value

def canonical_expr(expr):
    out=Counter()
    for w,n in expr.items():
        sig,v=canonical(w)
        if sig:out[v]+=n*sig
    return {w:n for w,n in out.items() if n}

def shape(w):
    def cyc(seq):return min(v[j:]+v[:j] for v in (seq,seq[::-1]) for j in range(len(v)))
    return tuple(sorted(cyc(tuple((a+b,len(s)) for a,b,s in tr)) for tr in w))

@lru_cache(None)
def raise_word(w,k):
    out=Counter()
    for t,tr in enumerate(w):
        for j,(a,b,s) in enumerate(tr):
            if k==3:
                if not b:continue
                factor=a+1;new=(a+1,b-1,s)
            else:
                src,dst=((1,0) if k==1 else (2,1))
                if src not in s or dst in s:continue
                factor=1;new=(a,b,tuple(sorted(dst if x==src else x for x in s)))
            sig,v=canonical(w[:t]+(tr[:j]+(new,)+tr[j+1:],)+w[t+1:])
            if sig:out[v]+=sig*factor
    return {v:n for v,n in out.items() if n}

def highest(words):
    equations=[{(k,v):n for k in (1,2,3) for v,n in raise_word(w,k).items()} for w in words]
    keys=sorted({k for r in equations for k in r});indices={k:i for i,k in enumerate(keys)}
    a=fmpz_mat(len(keys),len(words))
    for j,row in enumerate(equations):
        for k,n in row.items():a[indices[k],j]=n
    null,rank=a.nullspace()
    for j in range(rank):
        co=[int(null[i,j]) for i in range(len(words))]
        divisor=math.gcd(*co)
        if next(n for n in co if n)<0:divisor=-divisor
        co=[n//divisor for n in co]
        assert (a*fmpz_mat([[n] for n in co])).is_zero()
        yield {w:n for w,n in zip(words,co) if n}

@lru_cache(200000)
def qword(w):return tuple(canonical_expr(fw.Q_multi(w)).items())
def q(expr):
    out=Counter()
    for w,n in expr.items():
        for v,c in qword(w):out[v]+=n*c
    return {v:n for v,n in out.items() if n}

def basis(mu,length):
    blocks=defaultdict(set)
    for content in contents(mu,length):
        for w in multitraces_with_content(content):
            sig,v=canonical(w)
            if sig:blocks[shape(v)].add(v)
    rows=defaultdict(list)
    for sh in sorted(blocks):
        for row in highest(sorted(blocks[sh])):
            degree={nferm(w) for w in row}
            assert len(degree)==1
            rows[degree.pop()].append(row)
    return dict(rows),sum(map(len,blocks.values()))

def physical_sector(mu,length):
    theta=sum(mu[2:])
    return (2*mu[0]+theta-length,2*mu[1]+theta-length)+tuple(mu[0]+mu[1]+2*x for x in mu[2:])
