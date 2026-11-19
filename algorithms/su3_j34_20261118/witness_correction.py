"""Use identities on the *right-hand side*, then solve complete incoming incidence
on the union of the supports of 733 outgoing coefficient functionals.

All coordinates and all coefficients are rational; outside this union the same
trace expression is not asserted to solve the physical correction equation.
"""
import argparse, hashlib, os
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import chain
from pathlib import Path
import replay_formal_blocks as s
import finite_trace_identity as f
from fortuity.full.words import Q_multi
from projected_finite_correction import axpy, rref

@lru_cache(20000)
def image(w):return {s.packed(v):n for v,n in Q_multi(s.expanded(w)).items() if n}

def incoming(w):return s.parents(s.expanded(w))

def run(rhs, replay):
    data=s.read(replay)['coefficients']; outs=[s.unpack_hex(r['outgoing']) for r in data];weights=[];which={}
    for j,u in enumerate(outs):
        z={w:k for w in sorted(incoming(u)) if (k:=image(w).get(u,0))}
        for w in z:
            assert w not in which and s.block(w)==s.block(u)
            which[w]=j
        weights.append(z)
    blocks=[{} for _ in outs]; h=hashlib.sha256();n=0
    for w,c in s.stream_rows(rhs,h):
        n+=1
        if w in which:blocks[which[w]][w]=Fraction(c)
    assert h.hexdigest()==s.INPUT and n==517468
    report=[]; solution={}; full_R={}; no_pair=[];partial=[]
    for j,(wt,b) in enumerate(zip(weights,blocks)):
        pairing=sum((wt[w]*c for w,c in b.items()),Fraction());assert pairing==Fraction(data[j]['coefficient'])!=0
        parents=sorted(set().union(*(incoming(w) for w in wt)))
        basis=[{w:k for w,k in image(v).items() if w in wt} for v in parents]
        for v,col in zip(parents,basis):
            assert all(which.get(w,j)==j for w in image(v));assert not sum(wt[w]*k for w,k in col.items())
        # An identity hitting this witness and having no coordinates in another
        # witness support. Counting multiplicity and signs on *all* RHS parents.
        chosen=None
        for w in sorted(b):
            for pieces in f.choices(w):
                R=f.identity(pieces)
                if any(which[t]!=j for t in R if t in which):continue
                H={t:v for t,v in R.items() if t in wt}; m=sum(wt[t]*v for t,v in H.items())
                if not m:continue
                chosen=(pieces,R,H,m); break
            if chosen:break
        if not chosen: no_pair.append(j); continue
        pieces,R,H,m=chosen; lam=pairing/m; rem=b.copy();axpy(rem,H,-lam)
        assert sum(wt[w]*c for w,c in rem.items())==0
        raw,x0,_=rref(basis,b);adj,x,_=rref(basis,rem)
        assert raw['rank']==adj['rank'] and raw['augmented_rank']==raw['rank']+1
        resid=rem.copy()
        for i,c in x.items():axpy(resid,basis[i],-c)
        if resid:partial.append(j)
        else:
            for i,c in x.items():
                assert parents[i] not in solution;solution[parents[i]]=c
        axpy(full_R,R,lam)
        report.append(dict(witness=j,coordinates=len(wt),incident_predecessors=len(parents),rank=raw['rank'],formal_augmented=raw['augmented_rank'],identity_augmented=adj['augmented_rank'],pieces=s.hexes(pieces),Q0R_pairing=m,rhs_pairing=str(pairing),multiplier=str(lam),primitive_terms=len(x) if not resid else None))
    print('coordinates',len(which),'rank',sum(z['rank'] for z in report),'identities',len(report),'misses',no_pair,'nonzero residual blocks',partial,flush=True)
    # Replay the actual joint solution, not just the 733 scalar pairings.
    left=Counter(full_R)
    for w,c in solution.items():axpy(left,image(w),c)
    filed={k:v for b in blocks for k,v in b.items()}; residual={w:filed.get(w,0)-left.get(w,0) for w in which if filed.get(w,0)!=left.get(w,0)}
    for w in residual:assert which[w] in set(no_pair+partial)
    outside=None
    if not residual:
        h=hashlib.sha256();d=left.copy();count=0;first=None
        for w,k in s.stream_rows(rhs,h):
            y=Fraction(k)-d.pop(w,0)
            if y:
                count+=1
                if first is None:first=[s.hexes(w),str(k),str(y)]
        count+=sum(bool(v) for v in d.values());assert h.hexdigest()==s.INPUT and count and first
        outside=dict(formal_residual_after_this_fixed_decomposition=count,first_rhs_word_and_residual=first,extra_off_source_words=sum(bool(v) for v in d.values()))
    return dict(rhs_sha256=s.INPUT,coordinates=len(which),rhs_nonzero_coordinates=len(filed),witness_blocks=len(outs),separated_functionals=len(report),unpaired=no_pair,unsolved_blocks=partial,
                full_formal_projected_rank=sum(z['rank'] for z in report) if not no_pair else None,
                identity_adjusted_nonzero_projected_residuals=len(residual),fixed_solution_terms=len(solution),combined_identity_terms=len(full_R),
                identity_sha256=hashlib.sha256(sjson(full_R)).hexdigest(),solution=[[s.hexes(w),str(v)] for w,v in sorted(solution.items())],blocks=report,off_projection=outside)

def sjson(z):
    import json
    return json.dumps([[s.hexes(w),str(v)] for w,v in sorted(z.items())],separators=(',',':')).encode()
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');p.add_argument('--replay',type=Path,default=s.HERE/'finite_colour/replayed_733.json.gz');p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__;d=run(a.rhs,a.replay);s.write(a.out,d);print({k:v for k,v in d.items() if k not in ('blocks','solution')},flush=True)
