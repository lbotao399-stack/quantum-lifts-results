"""Re-test archived commuting substitutions on a different rational first correction.

The sign is Q0 O1 + Q1 O0.  Each parent is checked from its graded split and
all matrix commutators are checked before specializing the exterior algebra.
A nonzero value excludes only the displayed set of permitted parent contents.
"""
import argparse, hashlib, json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import parent_frontier as front
import replay_centralizer_sweep as rc

s=front.s
BASE=s.HERE/'joint_70'
CANDIDATES=('000000000001050505080f','0000000000010105050716','0000000000010105050e0f')

def content(w):return bytes(sorted(b''.join(w))).hex()

def residuals(rhs,vector,targets):
    raw=s.read(vector);digest=hashlib.sha256(json.dumps(raw,separators=(',',':')).encode()).hexdigest()
    rows={t:{} for t in targets};h=hashlib.sha256();n=0;nr=Counter()
    for w,c in s.stream_rows(rhs,h):
        n+=1;t=content(w)
        if t in rows:
            assert w not in rows[t];rows[t][w]=-Fraction(c);nr[t]+=1
    assert n==517468 and h.hexdigest()==s.INPUT
    needed={g for t in targets for g in rc.old.split_pairs(t)}
    for tr,c in raw:
        w=s.unpack_hex(tr)
        if content(w) not in needed:continue
        v=Fraction(c)
        for u,k in rc.old.Q_multi(s.expanded(w)).items():
            a=s.packed(u);t=content(a)
            if t in rows:rows[t][a]=rows[t].get(a,0)+k*v
    return digest,nr,{t:{w:c for w,c in d.items() if c} for t,d in rows.items()}

def mm(a,b):
    return tuple(a[i]*b[j]+a[i+1]*b[3+j]+a[i+2]*b[6+j] for i in (0,3,6) for j in (0,1,2))

def sign(a,b):return -1 if sum((b&((1<<i)-1)).bit_count() for i in range(a.bit_length()) if a>>i&1)%2 else 1

def evaluate(t,terms,entry,allowed,control=True):
    wit=entry['witness']; bits=[i for i in range(7) if wit['mask']>>i&1];n=len(bits);full=(1<<n)-1;letters=set(bytes.fromhex(t));odd=rc.ft.ODD
    assert n==sum(odd[i] for i in bytes.fromhex(t)) and allowed<=set(entry['allowed']);edges={rc.old.split_pairs(t)[g][1:3] for g in allowed}
    assert edges<=set(map(tuple,entry['edges']));table={(i,j):tuple(sum(wit['matrices'][f'{i:02x}:{j}'],[])) for i in letters for j in (range(7) if odd[i] else [-1])}
    assert len(table)==len(wit['matrices']) and all(a[0]+a[4]+a[8]==0 for a in table.values());checks=0
    for i,j in sorted(edges):
        for h in (range(7) if odd[i] else [-1]):
            for k in (range(7) if odd[j] else [-1]):
                a,b=table[i,h],table[j,k];assert mm(a,b)==mm(b,a);checks+=1
    # Right multiplication preserves trace order. The sign counts reversals of exterior generators.
    I=(1,0,0,0,1,0,0,0,1);z=(0,)*9
    @lru_cache(maxsize=64000)
    def trace_prefix(w):
        if not w:return {0:I}
        i=w[-1];out={}
        for old,A in trace_prefix(w[:-1]).items():
            for m,B in ([(1<<j,table[i,h]) for j,h in enumerate(bits)] if odd[i] else [(0,table[i,-1])]):
                if old&m:continue
                k=old|m;v=mm(A,B);e=sign(old,m);out[k]=tuple(x+e*y for x,y in zip(out.get(k,z),v))
        return {k:v for k,v in out.items() if any(v)}
    @lru_cache(maxsize=24000)
    def trace(w):return {k:a[0]+a[4]+a[8] for k,a in trace_prefix(w).items() if a[0]+a[4]+a[8]}
    def word(w):
        d={0:1}
        for tr in w:
            out={}
            for a,u in d.items():
                for b,v in trace(tr).items():
                    if not a&b:out[a|b]=out.get(a|b,0)+sign(a,b)*u*v
            d=out
        return d.get(full,0)
    if control:
        M={(i,j):tuple(tuple(A[k:k+3]) for k in (0,3,6)) for i in letters for j in (range(n) if odd[i] else [-1]) for A in [table[i,bits[j] if j>=0 else -1]]};checked=0
        for w in sorted(terms):
            val=word(w)
            if not val:continue
            assert rc.ft.exterior_value(w,M,n)==val;checked+=1
            if checked==2:break
        assert checked or not any(word(w) for w in terms)
    value=sum((c*v for w,c in terms.items() if (v:=word(w))),Fraction());trace.cache_clear();trace_prefix.cache_clear();return value,checks

def run(out):
    z=BASE/'support_cut';c=s.read(z/'preserve_296.census.json');sweep=s.read(BASE/'centralizer_sweep.json.gz');F,B=front.incidence();passed=set(c['exact_zero']);rhs=s.read(z/'rhs_comparison.json')['formal_rhs_words_by_content'];assert len(passed)==313 and not c['unresolved'];raw_dig,nr,rows=residuals(s.HERE/'negative_Q1.jsonl.gz',z/'preserve_296.vector.json.gz',CANDIDATES);assert raw_dig==c['vector_sha'];report={}
    for t in CANDIDATES:
        assert nr[t]==rhs[t]>0 and t in c['mathematically_nonzero'];parents=B[bytes.fromhex(t)];allowed={g.hex() for g in parents if not passed&{k.hex() for k in F[g]}};prev=rc.test(sweep['targets'][t],0);v,k=evaluate(t,rows[t],prev,allowed);dims={g:front.j.E.parity_frame_dims(bytes.fromhex(g))[0] for g in sorted(allowed)};assert v.denominator%1000003;mod=v.numerator*pow(v.denominator,-1,1000003)%1000003
        report[t]=dict(rhs_trace_words=nr[t],residual_trace_words=len(rows[t]),target_dimension=front.j.E.parity_frame_dims(bytes.fromhex(t))[0],allowed_parent_dimensions=dims,split_edges=sorted({parents[bytes.fromhex(g)] for g in allowed}),commutators_checked=k,mask=prev['witness']['mask'],prior_residue_mod_prime=prev['witness']['value'],rational_value=str(v),residue_mod_prime=mod,excluded=bool(v));print(t,sum(dims.values()),len(rows[t]),mod,flush=True)
    result=dict(prime=1000003,source_sha256=raw_dig,rhs_sha256=s.INPUT,sweep_sha256=hashlib.sha256((BASE/'centralizer_sweep.json.gz').read_bytes()).hexdigest(),preserved_contents=len(passed),targets=report,sign='Q0 O1 + Q1 O0',scope='Restricted to parents whose individual splits avoid the 313 successes. Correlated changes of touching parents are not tested.');s.write(out,result);return result
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path);args=ap.parse_args();assert __debug__;run(args.out)
