"""Extract every Q0(-Q1 O0) word in an exact outgoing letter content.

The input to the companion integer colour-polynomial evaluator is plain text:
number of words, then integer coefficient, trace count, and hex-coded traces.
"""
import argparse, hashlib
from collections import Counter
from itertools import combinations
from fractions import Fraction
from math import lcm
from pathlib import Path
import replay_formal_blocks as s
from fortuity.full.words import Q_multi

DEFAULT='000000000000040505050710'
def incident(content):
    found=set()
    for i,j in combinations(range(len(content)),2):
        a,b,S=s.ALPH[content[i]];d,e,T=s.ALPH[content[j]]; merged=a+d,b+e,tuple(sorted(S+T))
        if set(S).isdisjoint(T) and merged in s.LETTER:
            u=bytearray(content);del u[j];del u[i];u.append(s.LETTER[merged]);found.add(bytes(sorted(u)))
    return found

def run(rhs,keys,out):
    parents={c:incident(c) for c in keys}; needed=set.union(*parents.values());total=0;seen=Counter();h=hashlib.sha256();data={k:Counter() for k in keys}
    for w,coef in s.stream_rows(rhs,h):
        total+=1;cw=bytes(sorted(b''.join(w)))
        if cw not in needed: continue
        seen[cw]+=1;coef=Fraction(coef)
        for v,n in Q_multi(s.expanded(w)).items():
            v=s.packed(v);cv=bytes(sorted(b''.join(v)))
            if cv in data:data[cv][v]+=coef*n
    assert total==517468 and h.hexdigest()==s.INPUT
    out.mkdir(exist_ok=True,parents=True);meta={}
    for c in sorted(keys):
        z={w:k for w,k in data[c].items() if k};den=lcm(*(k.denominator for k in z.values()));raw=(str(len(z))+'\n'+''.join(f"{int(v*den)} {len(w)} {' '.join(s.hexes(w))}\n" for w,v in sorted(z.items()))).encode();assert all((v*den).denominator==1 for v in z.values());tag=c.hex();(out/(tag+'.txt')).write_bytes(raw)
        meta[tag]=dict(input_rhs_sha256=s.INPUT,input_rhs_terms=total,incident_rhs_terms={x.hex():seen[x] for x in sorted(parents[c])},formal_terms=len(z),common_denominator=den,numerator_sum_abs=sum(abs(int(v*den)) for v in z.values()),polynomial_input_sha256=hashlib.sha256(raw).hexdigest())
    s.write(out/'meta.json',meta);return meta
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contents',nargs='+',default=[DEFAULT]);p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__;print(run(a.rhs,set(map(bytes.fromhex,a.contents)),a.out),flush=True)
