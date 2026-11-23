"""Generate universal Q0 parents by merging adjacent native letters."""
from fast_paths import *
from fortuity.full.letters import Q_letter
from collections import defaultdict
import gzip
import json
import pickle
import time


@lru_cache(maxsize=4096)
def merge_letter(left,right):
    if set(left[2]) & set(right[2]):
        return None
    letter=(left[0]+right[0],left[1]+right[1],tuple(sorted(left[2]+right[2])))
    coefficient=sum(c for c,a,b in Q_letter(letter) if (a,b)==(left,right))
    return letter if coefficient else None


def parents(word):
    result=set()
    for k,trace in enumerate(word):
        if len(trace)<3:continue
        for i in range(len(trace)):
            rotated=trace[i:]+trace[:i]
            letter=merge_letter(*rotated[:2])
            if letter is None:continue
            sign,source=pf.m.so_canon(word[:k]+((letter,)+rotated[2:],)+word[k+1:])
            if sign:result.add(source)
    return result


def block(word):
    return tuple(sorted((sum(l[0] for l in tr),sum(l[1] for l in tr),
                         *(sum(i in l[2] for l in tr) for i in range(3))) for tr in word))


def analyse(defect):
    groups=defaultdict(dict)
    for word,c in defect.items():groups[block(word)][word]=c
    records=[]
    for key,part in groups.items():
        source=set().union(*(parents(w) for w in part))
        target={v for w in source for v,c in pf.m.differential(w)}
        records.append(dict(trace_weights=key,defect_words=len(part),parents=len(source),
                            target_words=len(target),missing_words=len(set(part)-target)))
    records.sort(key=lambda x:x['parents'],reverse=True)
    return records


if __name__=='__main__':
    start=time.monotonic()
    folder=BASE_PACKAGE/'symplectic/sp5' if BASE_PACKAGE else HERE.parent/'fortuity_frontier_delivery_v2_20261123/symplectic/sp5'
    data=json.loads(gzip.decompress((folder/'sp5_new_fortuitous_operator.json.gz').read_bytes()))
    from fractions import Fraction
    import math
    polynomial={pf.m.decode(w):Fraction(c) for c,w in data['terms']}
    denominator=math.lcm(*(c.denominator for c in polynomial.values()))
    integral={w:int(c*denominator) for w,c in polynomial.items()}
    defect=pf.m.q(integral)
    records=analyse(defect)
    result=dict(rank=5,operator_terms=len(integral),defect_terms=len(defect),
                trace_weight_blocks=len(records),blocks=records,seconds=time.monotonic()-start)
    (HERE/'inverse_support_sp5.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='blocks'}),flush=True)
    print(json.dumps(records[:10],indent=2),flush=True)
