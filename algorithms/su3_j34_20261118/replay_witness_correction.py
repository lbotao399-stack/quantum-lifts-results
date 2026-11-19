"""Replay filed coefficients without solving a system. Use the separate packed Q0
implementation for the forward/reverse words, and the graded four-slot identity.
"""
import json, hashlib
from fractions import Fraction
from pathlib import Path
import replay_formal_blocks as s
import formal_trace_blocks as q
import finite_trace_identity as ft
from projected_finite_correction import axpy, content

D = s.HERE/'finite_colour'

def main():
    w = s.read(D/'witness_correction.json.gz'); u=s.read(D/'replayed_733.json.gz')['coefficients']
    identified={}; ids={}; projection=set(); expected={}; total=0
    assert not w['unsolved_blocks'] and not w['unpaired'] and len(u)==len(w['blocks'])==733
    for j,(record,target) in enumerate(zip(w['blocks'],u)):
        assert record['witness']==j
        t=s.unpack_hex(target['outgoing']); weights={p:n for p in q.parents(t) if (n:=q.q(p).get(t,0))}
        assert len(weights)==record['coordinates'] and projection.isdisjoint(weights); projection.update(weights)
        R=ft.identity(s.unpack_hex(record['pieces']));m=sum(v*weights.get(p,0) for p,v in R.items())
        assert m==record['Q0R_pairing']!=0 and Fraction(record['rhs_pairing'])==Fraction(target['coefficient'])==m*Fraction(record['multiplier'])
        axpy(ids,R,Fraction(record['multiplier']))
    serial=json.dumps([[s.hexes(p),str(v)] for p,v in sorted(ids.items())],separators=(',',':')).encode()
    assert len(ids)==w['combined_identity_terms'] and hashlib.sha256(serial).hexdigest()==w['identity_sha256']
    def apply(raw):
        result={}
        for word,v in raw:axpy(result,q.q(s.unpack_hex(word)),Fraction(v))
        return result
    lhs=apply(w['solution']);axpy(lhs,ids,1)
    four=s.read(D/'four_contents_first.json.gz'); targets=set(map(bytes.fromhex,four['contents'])); lhs4=apply(four['solution_or_partial']); assert len(w['solution'])==w['fixed_solution_terms'] and len(lhs4)>0
    h=hashlib.sha256();unrestricted=outside4=inside4=tested=inside=off=0; first=None
    for v,c in s.stream_rows(s.HERE/'negative_Q1.jsonl.gz',h):
        total+=1; k=Fraction(c);delta=k-lhs.pop(v,0);delta4=k-lhs4.pop(v,0)
        if v in projection:tested+=1;inside+=bool(delta)
        if delta and first is None:first=[s.hexes(v),str(c),str(delta)]
        unrestricted+=bool(delta);outside4+=bool(delta4)
        if content(v) in targets:inside4+=bool(delta4);off+=1
    assert h.hexdigest()==w['rhs_sha256']==s.INPUT and total==517468
    inside+=sum(v in projection for v in lhs);inside4+=sum(content(v) in targets for v in lhs4);outside4+=len(lhs4);unrestricted+=len(lhs)
    assert inside==w['identity_adjusted_nonzero_projected_residuals']==0 and len(projection)==w['coordinates'] and tested==w['rhs_nonzero_coordinates']
    assert four['projected_terms']==off and inside4==0
    assert (unrestricted,first,len(lhs))==(w['off_projection']['formal_residual_after_this_fixed_decomposition'],w['off_projection']['first_rhs_word_and_residual'],w['off_projection']['extra_off_source_words'])
    return dict(rhs_terms=total,rhs_sha256=h.hexdigest(),projected_coordinates=len(projection),rhs_nonzero_in_projection=tested,projected_residuals=inside,unrestricted_formal_residuals=unrestricted,four_contents_rhs_terms=off,four_contents_projected_residuals=inside4,four_contents_unrestricted_formal_residuals=outside4)

if __name__=='__main__':
    assert __debug__; print(json.dumps(main(),sort_keys=True,indent=2))
