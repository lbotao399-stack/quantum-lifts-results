"""Exact rational free-trace tests for -Q1 O34, before the SU(3) trace identities.

The index of a letter is its position in sorted(letters_up_to(1,1)). Its order
therefore agrees with the archived graded cyclic canonicalization. No projection
onto a target profile is made. All internal integers store numerators over 120,
the common denominator of the input; exported coefficients are rational.
"""
import argparse, gzip, hashlib, json, sys
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
from itertools import combinations
HERE=Path(__file__).resolve().parent
RT=HERE.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
sys.path.insert(0,str(RT));import sp_environment
from fortuity.full.letters import Q_letter,is_ferm,letters_up_to
from fortuity.full.words import Q_multi
ALPH=sorted(letters_up_to(1,1));ID={l:i for i,l in enumerate(ALPH)};ODD=tuple(map(is_ferm,ALPH));DEGREES=tuple((a,b,*(int(i in s) for i in range(3))) for a,b,s in ALPH)
SPLITS=tuple(tuple((c,bytes((ID[u],ID[v]))) for c,u,v in Q_letter(l)) for l in ALPH)
MERGE={p:i for i,terms in enumerate(SPLITS) for _,p in terms};assert len(MERGE)==sum(map(len,SPLITS))
INPUT='9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'
DEN=120
@lru_cache(200000)
def par(t):return sum(ODD[x] for x in t)%2
@lru_cache(160000)
def canon(t):
    if len(t)<2:return 0,b''
    even=not par(t);rot=t;best=t;sgn=v=1
    for i in range(1,len(t)+1):
        if even and ODD[rot[0]]:v=-v
        rot=rot[1:]+rot[:1]
        if rot==t:
            if v<0:return 0,b''
            break
        if rot<best:best,sgn=rot,v
    return sgn,best

def multi(w):
    odd=[t for t in w if par(t)]
    if len(set(odd))!=len(odd):return 0,()
    return (-1)**sum(x>y for x,y in combinations(odd,2)),tuple(sorted(w))

@lru_cache(15000)
def qtrace(t):
    z=Counter();s=1
    for i,x in enumerate(t):
        for c,v in SPLITS[x]:
            e,r=canon(t[:i]+v+t[i+1:]);
            if e:z[r]+=s*c*e
        if ODD[x]:s=-s
    return tuple((c,r) for r,c in z.items() if c)

def q(w):
    d=Counter();s=1
    for i,t in enumerate(w):
        for c,r in qtrace(t):
            e,v=multi(w[:i]+(r,)+w[i+1:]);
            if e:d[v]+=s*c*e
        if par(t):s=-s
    return {w:c for w,c in d.items() if c}

@lru_cache(30000)
def unsplit(t):
    if len(t)<3:return ()
    out=set()
    for i in range(len(t)):
        rot=t[i:]+t[:i];x=MERGE.get(rot[:2]);
        if x is not None:
            s,v=canon(bytes([x])+rot[2:]);
            if s:out.add(v)
    return tuple(sorted(out))

def parents(w):
    out=set()
    for i,t in enumerate(w):
        for v in unsplit(t):
            s,r=multi(w[:i]+(v,)+w[i+1:]);
            if s:out.add(r)
    return out

@lru_cache(300000)
def trace_degree(t):return tuple(map(sum,zip(*(DEGREES[x] for x in t))))

def key(w):return tuple(sorted(trace_degree(t) for t in w))

def display(w):return [t.hex() for t in w]
def unpack(w):return tuple(bytes.fromhex(t) for t in w)
@lru_cache(4096)
def numerator(c):
    f=Fraction(c)*DEN;assert f.denominator==1;return int(f)
@lru_cache(250000)
def intern(t):return bytes(t)

def sources(path):
    blocks=defaultdict(dict);h=hashlib.sha256();num=0
    with gzip.open(path,'rb') as f:
        for line in f:
            h.update(line);c,v=json.loads(line);w=tuple(intern(tuple(ID[a,b,tuple(s)] for a,b,s in t)) for t in v);b=blocks[key(w)];assert w not in b;b[w]=numerator(c);num+=1
    assert h.hexdigest()==INPUT and num==517468 and len(blocks)==733
    return blocks

def check_word(w):
    dec=lambda w:tuple(tuple(ALPH[x] for x in t) for t in w)
    assert {dec(k):v for k,v in q(w).items()}==dict(Q_multi(dec(w)))
    for v in q(w):assert w in parents(v) and key(w)==key(v)

def closure(rhs):
    # Q0^2=0 makes a single nonzero outgoing coefficient a whole-block no-go.
    # Inverse adjacent-letter merges enumerate every incident input, wrap included.
    tried=set()
    for w in list(rhs)[:4]:
        for v in q(w):
            if v in tried:continue
            tried.add(v);pre=sorted(set(rhs).intersection(parents(v)));arr=[(x,q(x).get(v,0)) for x in pre];n=sum(rhs[x]*t for x,t in arr)
            if n:return {'outgoing':display(v),'coefficient':str(Fraction(n,DEN)),'contributors':[[display(x),str(Fraction(rhs[x],DEN)),t] for x,t in arr if t]},False
    z=Counter()
    for w,c in rhs.items():
        for v,t in q(w).items():z[v]+=t*c
    for v,n in z.items():
        if n:
            arr=[(x,q(x).get(v,0)) for x in sorted(parents(v)) if x in rhs];assert n==sum(rhs[x]*t for x,t in arr)
            return {'outgoing':display(v),'coefficient':str(Fraction(n,DEN)),'contributors':[[display(x),str(Fraction(rhs[x],DEN)),t] for x,t in arr if t]},True
    return None,True

def primitive(rhs,cap):
    # Look in the support-adjacent predecessor span, including every output coordinate
    # of each selected predecessor. Failure here is deliberately unclassified.
    from flint import fmpq_mat
    W=sorted({v for w in rhs for v in parents(w) if par_degree_ok(v)})
    if len(W)>cap:return 'over_limit',None,(len(W),None)
    images=[q(v) for v in W];coords=sorted(set(rhs).union(*(set(i) for i in images)));size=(len(W),len(coords))
    if len(coords)*max(1,len(W))>cap**2:return 'over_limit',None,size
    A=fmpq_mat([[row.get(w,0) for row in images]+[rhs.get(w,0)] for w in coords]);R,r=A.rref();x=[Fraction(0)]*len(W)
    for i in range(r):
        j=next(j for j in range(len(W)+1) if R[i,j]);
        if j==len(W):return 'unsolved_selected',None,size
        x[j]=Fraction(str(R[i,len(W)]))/DEN
    z=Counter()
    for c,row in zip(x,images):
        if c:
            for w,n in row.items():z[w]+=c*n
    assert {w:c for w,c in z.items() if c}=={w:Fraction(c,DEN) for w,c in rhs.items()}
    return 'exact',[[display(v),str(c)] for v,c in zip(W,x) if c],size

def par_degree_ok(w):return sum(map(len,w))==10

def save(path,data):
    path.parent.mkdir(exist_ok=True,parents=True);buf=(json.dumps(data,separators=(',',':'))+'\n').encode();path.write_bytes(gzip.compress(buf,compresslevel=9,mtime=0) if path.suffix=='.gz' else buf)

def run(a):
    blocks=sources(a.rhs);rows=[];counts=Counter();totals=Counter();print('partitioned',len(blocks),flush=True)
    for i,(k,rhs) in enumerate(sorted(blocks.items(),key=lambda item:(len(item[1]),item[0]))):
        check_word(next(iter(rhs))); w,full=closure(rhs);rec={'degrees':k,'rhs_terms':len(rhs),'closure_expanded':full}
        if w:kind='not_formally_closed';rec['witness']=w
        else:
            kind,prim,sz=primitive(rhs,a.cap);rec['selected_shape']=sz
            if prim:rec['primitive']=prim
        rec['classification']=kind;rows.append(rec);counts[kind]+=1;totals[kind]+=len(rhs)
        if len(rows)%100==0:print(len(rows),dict(counts),flush=True)
    assert sum(totals.values())==517468
    save(a.out,dict(input_stream_sha256=INPUT,letter_table=ALPH,source_length=10,target_length=11,degree_order=['z+','z-','theta1','theta2','theta3'],candidate_cap=a.cap,block_counts=dict(counts),term_counts=dict(totals),blocks=rows));print('final',counts,totals,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=HERE/'negative_Q1.jsonl.gz');p.add_argument('--out',type=Path,required=True);p.add_argument('--cap',type=int,default=220);a=p.parse_args();assert __debug__;run(a)
