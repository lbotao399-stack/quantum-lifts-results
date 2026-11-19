"""Enumerate every graded cyclic multitrace of one specified parent content.

Fix the smallest letter of a trace at its first position. Remove the trace's
multiset from the remaining letters, and append traces in lexicographic order.
Cyclic signs and repeated odd traces use the same classical differential rules.
This includes trace identities at N=3 as linear dependencies rather than imposing
an independence assumption on the parent generators.
"""
from functools import lru_cache
from itertools import product
from collections import Counter
import formal_trace_blocks as f

@lru_cache(maxsize=256)
def all_words(content):
    content=bytes(sorted(content));counts=Counter(content);labels=sorted(counts);traces={}
    def extensions(x,counts):
        if not x or (len(x)>=2 and f.canon(x)==(1,x)):
            t=bytes(sorted(x));traces.setdefault(t,[]).append(x)
        for i in labels:
            if counts[i] and (not x or i>=x[0]):
                counts[i]-=1;extensions(x+bytes([i]),counts);counts[i]+=1
    extensions(b'',counts)
    @lru_cache(None)
    def parts(rem):
        if not rem:return ((),)
        bag=Counter(rem);out=[]
        for mult in product(*(range(1 if i==rem[0] else 0,bag[i]+1) for i in labels)):
            n=sum(mult)
            if n<2 or n==len(rem)-1:continue
            piece=b''.join(bytes([i])*m for i,m in zip(labels,mult));rest=b''.join(bytes([i])*(bag[i]-m) for i,m in zip(labels,mult));tails=parts(rest)
            for tr in traces.get(piece,()):
                for w in tails:
                    if w and (tr>w[0] or (tr==w[0] and f.par(tr))):continue
                    out.append((tr,)+w)
        return tuple(sorted(out))
    words=parts(content);assert len(set(words))==len(words);return words
