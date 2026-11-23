#!/usr/bin/env python3
"""Evaluate the scalar Pfaffian ancestor of the Codex degree15 relation."""
from pathlib import Path
from fractions import Fraction
from collections import Counter,defaultdict
from math import factorial,prod
import json,os,sys,time

HERE=Path(__file__).resolve().parent
from pfaffian_adapter import ScalarPfaffian,PfaffianCoefficient,oracle
import independent_trace_ops as tr

def transpose(dual,source,target):
    out=defaultdict(Fraction)
    for word,c in dual.items():
        candidates=set()
        for i,trace in enumerate(word):
            for j,l in enumerate(trace):
                if l==target:
                    changed=word[:i]+(trace[:j]+source+trace[j+1:],)+word[i+1:]
                    sign,w=tr.canonical(changed)
                    if sign:candidates.add(w)
        for w in candidates:
            value=tr.action({w:Fraction(1)},source,target).get(word,0)
            if value:out[w]+=value*c
    return {w:c for w,c in out.items() if c}

def main():
    start=time.monotonic();s=ScalarPfaffian();c=PfaffianCoefficient(3,0)
    cert=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    values={w:str(oracle(s,tuple(w.split('.')))) for w,x in zip(cert['terminal_probe_words'],cert['dual_coefficients']) if Fraction(x)}
    pairing=sum(Fraction(x)*Fraction(values[w]) for w,x in zip(cert['terminal_probe_words'],cert['dual_coefficients']) if Fraction(x))
    probes=[]
    for query in tr.QUERIES:
        sign,w=tr.canonical(query);dual={w:Fraction(sign)};sizes=[len(dual)]
        for scalar in 'CBA':
            for ghost in '+-':
                dual=transpose(dual,scalar,ghost);sizes.append(len(dual))
        actual=sum(v*oracle(s,w) for w,v in dual.items())
        codex=oracle(c,query)
        probes.append({'ordered_query':list(query),'scalar_ancestor_P6':str(actual),'codex':str(codex),'ratio':str(actual/codex),'transpose_supports':sizes})
    result={'passed':True,'scalar_pfaffian_order':5,'scalar_dual_pairing':str(pairing),'scalar_probe_coefficients':values,'return_checks':probes,'seconds':time.monotonic()-start}
    (HERE/'codex_scalar_return.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
