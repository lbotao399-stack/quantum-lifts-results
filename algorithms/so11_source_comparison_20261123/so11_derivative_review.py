#!/usr/bin/env python3
"""Independent exact Python replay of two chosen SO11 stable derivatives.

Implements graded orthogonal traces, a rational scalar extremal projector,
and the six ordered strict odd substitutions. Does not import Dot's JS or
Codex's trace canonicalizer. This is not a quotient-by-Q0(I) calculation.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import argparse, hashlib, json, time

HERE=Path(__file__).resolve().parent
QUERIES=[('ABC','A+','A-','B+','B-','C+','C-'),('ABC','ABC','ABC','+-','+-','+-')]

def odd(c):return c in '+-'

@lru_cache(None)
def canonical_trace(word):
    if len(word)<2:return 0,''
    count=sum(map(odd,word));variants={}
    for source,initial in [(word,1),(word[::-1],(-1)**(len(word)+count*(count-1)//2))]:
        for j in range(len(word)):
            cut=sum(map(odd,source[:j]));sign=initial*(-1)**(cut*(count-cut))
            name=source[j:]+source[:j]
            variants.setdefault(name,set()).add(sign)
    name=min(variants)
    if len(variants[name])>1:return 0,''
    return next(iter(variants[name])),name

@lru_cache(None)
def canonical(traces):
    words=[];sign=1
    for trace in traces:
        local,name=canonical_trace(trace)
        if not local:return 0,()
        sign*=local;words.append(name)
    oddwords=[w for w in words if sum(map(odd,w))%2]
    if len(set(oddwords))<len(oddwords):return 0,()
    sign*=(-1)**sum(x>y for i,x in enumerate(oddwords) for y in oddwords[i+1:])
    return sign,tuple(sorted(words))

def action(poly,source,target):
    output=defaultdict(Fraction)
    operator_odd=odd(source)!=odd(target)
    for word,coefficient in poly.items():
        prefix=0
        for i,trace in enumerate(word):
            for j,letter in enumerate(trace):
                if letter==source:
                    changed=list(word);changed[i]=trace[:j]+target+trace[j+1:]
                    sign,name=canonical(tuple(changed))
                    if sign:output[name]+=coefficient*sign*(-1)**(prefix if operator_odd else 0)
                prefix+=odd(letter)
    return {w:c for w,c in output.items() if c}

def project(poly):
    for source,target,shift in [('C','B',1),('C','A',2),('B','A',1)]:
        result=defaultdict(Fraction,poly);raised=poly;coefficient=Fraction(1)
        for power in range(1,6):
            raised=action(raised,source,target)
            if not raised:break
            coefficient/=-power*(shift+power)
            lowered=raised
            for _ in range(power):lowered=action(lowered,target,source)
            for w,c in lowered.items():result[w]+=coefficient*c
        poly={w:c for w,c in result.items() if c}
    assert not action(poly,'B','A') and not action(poly,'C','B')
    return poly

def one_bank(directory):
    blocks={tuple(sorted(map(len,q))):defaultdict(Fraction) for q in QUERIES}
    hashes={}
    for name,coefficient in [('cycle.raw',Fraction(1)),('path3.raw',Fraction(-5)),('path3-alt.raw',Fraction(500,3))]:
        path=directory/name;hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        for line in path.read_text().splitlines():
            trace_string,value=line.split()
            if trace_string=='!':continue
            traces=tuple(trace_string.rstrip('.').split('.'))
            shape=tuple(sorted(map(len,traces)))
            if shape not in blocks:continue
            sign,word=canonical(traces)
            if sign:blocks[shape][word]+=coefficient*sign*int(value)
    result=[]
    for query in QUERIES:
        shape=tuple(sorted(map(len,query)))
        poly={w:c for w,c in blocks[shape].items() if c}
        poly=project(poly)
        for scalar in 'ABC':
            for ghost in '-+':poly=action(poly,scalar,ghost)
        sign,word=canonical(query)
        value=sign*poly.get(word,0)
        result.append({'query':list(query),'coefficient':str(value)})
    return {'input_sha256':hashes,'queries':result}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archived',type=Path,required=True)
    parser.add_argument('--fresh',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=HERE/'so11_derivative_review.json')
    args=parser.parse_args()
    start=time.monotonic()
    archived=one_bank(args.archived)
    fresh=one_bank(args.fresh)
    assert archived['queries']==fresh['queries']
    assert [x['coefficient'] for x in fresh['queries']]==['7372800','0']
    result={'passed':True,'independent_implementation':'Python Fraction with its own graded cyclic/reversal and multitrace canonicalization; imports neither JS engine nor repository trace implementation',
       'dot_source_derivative':'P^6 Pi G, because Q0 h_B(G)=G. No extra 1/2 factor.',
       'archived':archived,'fresh_enumeration':fresh,'seconds':time.monotonic()-start,
       'scope':'Two ordered stable-derivative coefficients only. No classical or quantum cohomology independence or full Q0(I) annihilation claimed.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
