"""Inventory all 500 contents of the first equation, including homogeneous ones.

The original right side is streamed in canonical words. For the fixed correction
use the upstream Q_multi, not the frame solver. An empty formal residual implies
exact equality before all finite-colour identities. A nonempty formal residual
alone places no restriction on the existence of a primitive.
"""
import argparse, hashlib
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import replay_scalar_layers as r
s,l=r.s,r.lay

def inventory(rhs,vector,emit=None):
    domain=set(l.contents(11,0,l.CHARGES)); assert len(domain)==500
    words=s.read(vector); contrib=defaultdict(dict); rhs_seen=Counter(); image=Counter(); remaining=Counter(); h=hashlib.sha256(); n=0
    for w,c in words:
        for u,co in r.Q_multi(s.expanded(s.unpack_hex(w))).items():
            v=s.packed(u);g=r.content(v);assert g in domain
            l.p.axpy(contrib[g],{v:co},Fraction(c))
    for g,row in contrib.items():image[g]=len(row)
    def deliver(g,w,val):
        if val:
            remaining[g]+=1
            if emit is not None:emit(g,w,val)
    for w,c in s.stream_rows(rhs,h):
        g=r.content(w);assert g in domain;n+=1;rhs_seen[g]+=1
        deliver(g,w,contrib[g].pop(w,0)-Fraction(c))
    assert n==517468 and h.hexdigest()==s.INPUT
    for g,data in contrib.items():
        for w,val in data.items():deliver(g,w,val)
    return dict(scope='Full length-eleven charge sector, fixed first correction. Formal zero is sufficient, not necessary, for physical zero.',rhs_sha256=s.INPUT,source_sha256=hashlib.sha256(s.json.dumps(words,separators=(',',':')).encode()).hexdigest(),source_words=len(words),target_contents=len(domain),rhs_terms=n,nonempty_rhs_contents=len(rhs_seen),nonempty_image_contents=sum(bool(t) for t in image.values()),nonempty_formal_residual_contents=len(remaining),contents={g.hex():[rhs_seen[g],image[g],remaining[g]] for g in sorted(domain)},column_order=['formal_rhs_words','formal_Q0_vector_words','formal_residual_words'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');p.add_argument('--vector',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__;x=inventory(a.rhs,a.vector);s.write(a.out,x);print({k:v for k,v in x.items() if k not in ('contents','scope')})
