"""Extend the top correction with the full word list of one isolated parent content.

Prioritize new target words already reached by the previous parent list. Then
try complete trace shapes and a fixed-seed sample. No conclusion about the
full differential is drawn from a modular gap in a selected subset.
"""
import argparse,random
from pathlib import Path
from collections import Counter
from fractions import Fraction
import content_parents
import scalar_filtration as l
s,f,p,t=l.s,l.f,l.p,l.tri
SOURCE=bytes.fromhex('00000000010205050710')

def construct(work,h,rhsfile,extra,mode,seed):
    hi={k for k in l.contents(11,0,l.CHARGES) if k.count(0)>=6};keep=hi|set(t.CONTENTS);rhs=p.rhs_rows(rhsfile,{t.MID});high=l.load_vector(h)
    for w,c in high.items():p.axpy(rhs,{u:a for u,a in f.q(w).items() if p.content(u)==t.MID},-c)
    ws=content_parents.all_words(SOURCE);pool={}
    for w in ws:
        r={u:k for u,k in f.q(w).items() if p.content(u) in keep}
        if r:assert all(p.content(u)==t.MID for u in r);pool[w]=r
    sel={w for u in rhs for w in f.parents(u) if p.content(w)==SOURCE and w in pool};original=len(sel);hit=set(rhs);steps=[]
    while True:
        more={w for w in pool if w not in sel and not hit.isdisjoint(pool[w])};prior=len(hit);sel.update(more);hit.update(u for w in sel for u in pool[w]);steps.append([len(sel),len(hit)])
        if not more and len(hit)==prior:break
    more=sorted(pool.keys()-sel);rng=random.Random(seed)
    if mode=='sparse':
        rng.shuffle(more);more.sort(key=lambda w:-len(w))
    else:rng.shuffle(more)
    sel.update(more[:extra]);im={w:pool[w] for w in sorted(sel)};z=l.recipe(work,im,rhs,{t.MID},'Complete enumeration of one isolated parent content; selected extra parents and formal incidence closure.');s.write(work/'selection.json',dict(content=SOURCE.hex(),enumerated=len(ws),nonzero_on_projection=len(pool),old_incident=original,closure=steps,extra=extra,ordering=mode,seed=seed,selected=len(sel),target_trace_words=len(z['word'])));return z
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--work',type=Path,required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--top',type=Path,default=s.HERE/'scalar_layers/top_vector.json.gz');a.add_argument('--extra',type=int,default=0);a.add_argument('--seed',type=int,default=340170);a.add_argument('--ordering',choices=['random','sparse'],default='random');b=a.parse_args();assert __debug__ and b.extra>=0;construct(b.work,b.top,b.rhs,b.extra,b.ordering,b.seed)
