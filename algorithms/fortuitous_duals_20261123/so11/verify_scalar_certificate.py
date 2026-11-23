#!/usr/bin/env python3
"""Independent standard-library rational elimination of the compact certificate."""
from fractions import Fraction
from pathlib import Path
import json,time

HERE=Path(__file__).resolve().parent

def rank(rows):
    a=[[Fraction(v) for v in row] for row in rows];r=0
    for c in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][c]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for k in range(r+1,len(a)):
            v=a[k][c]
            if v:a[k]=[x-v*y for x,y in zip(a[k],a[r])]
        r+=1
        if r==len(a):break
    return r

def main():
    t=time.monotonic();r=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    rs=rank(r['source_basis_probe_matrix']);rt=rank(r['terminal_basis_probe_matrix']);rq=rank(r['Q0_relation_matrix_on_probes'])
    assert (rs,rt,rq)==(21,17,16)
    dual=list(map(Fraction,r['dual_coefficients']))
    assert sum(bool(v) for v in dual)==6
    pairs=[sum(Fraction(v)*c for v,c in zip(row,dual)) for row in r['Q0_relation_matrix_on_probes']]
    assert not any(pairs)
    bpairs=[sum(Fraction(v)*c for v,c in zip(row,dual)) for row in r['terminal_basis_probe_matrix']]
    assert bpairs==[0]*16+[1]
    g=r['dot_seed_pairings'];assert Fraction(g['cycle'])==1 and Fraction(g['path3'])==Fraction(-1,2) and Fraction(g['path3-alt'])==Fraction(-1,50)
    score=Fraction(g['cycle'])-5*Fraction(g['path3'])+Fraction(500,3)*Fraction(g['path3-alt']);assert score==Fraction(1,6)
    out={'passed':True,'implementation':'independent Python standard-library Fraction elimination',
         'source_basis_rank':rs,'terminal_basis_rank':rt,'Q0_relation_rank':rq,
         'relation_quotient_dimension':rt-rq,'dual_support':6,'zero_Q0_pairings':len(pairs),
         'Dot_scalar_Gram_combination_pairing':str(score),'seconds':time.monotonic()-t,
         'scope':'Checks every entry of the compact rational matrices, all21annihilation identities, and the6coordinate dual. Character completeness and membership as genuine Gram identities are separate structural inputs reconstructed by the full driver.'}
    (HERE/'independent_scalar_certificate_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
