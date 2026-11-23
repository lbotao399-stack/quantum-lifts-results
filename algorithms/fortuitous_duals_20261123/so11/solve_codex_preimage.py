#!/usr/bin/env python3
"""Exact scalar Pfaffian ancestor coordinates and a constructive Q0(I) identity.

The full return identity is a separate, independently audited input. This
program certifies the scalar identity only and never infers finite exactness.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,time
from flint import fmpq_mat
from codex_scalar_return import ScalarPfaffian,oracle

HERE=Path(__file__).resolve().parent

def main():
    start=time.monotonic()
    cert=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    words=cert['terminal_probe_words'];pf=ScalarPfaffian()
    codex=[oracle(pf,tuple(w.split('.'))) for w in words]
    dot=[]
    raw={}
    for label in ['cycle','path3','path3-alt']:
        r=json.loads((HERE/'terminal_factor5'/f'{label}.json').read_text())
        raw[label]={w:Fraction(c)/int(r['denominator']) for w,c in r['projected']}
    dot=[raw['cycle'].get(w,0)-5*raw['path3'].get(w,0)+Fraction(500,3)*raw['path3-alt'].get(w,0) for w in words]
    dual=list(map(Fraction,cert['dual_coefficients']))
    assert sum(x*y for x,y in zip(dual,codex))==Fraction(-189,128)
    assert sum(x*y for x,y in zip(dual,dot))==Fraction(1,6)
    wanted=[9*d+Fraction(64,63)*c for d,c in zip(dot,codex)]
    b=fmpq_mat(cert['terminal_basis_probe_matrix'])
    x=b.transpose().solve(fmpq_mat([[str(v)] for v in wanted]))
    assert x[16,0]==0
    assert b.transpose()*x==fmpq_mat([[str(v)] for v in wanted])
    terms=[{'source_Gram_index':index,'coefficient':str(x[i,0])} for i,index in enumerate(cert['terminal_basis_Q0_source_indices']) if x[i,0]]
    result={
       'passed':True,
       'terminal_probe_words':words,
       'scalar_Pfaffian_ancestor_coefficients':list(map(str,codex)),
       'Dot_scalar_combination_coefficients':list(map(str,dot)),
       'combined_coefficients':list(map(str,wanted)),
       'combination':{'Dot_scalar_Gram_combination':'9','Codex_scalar_Pfaffian_ancestor':'64/63'},
       'dual_pairings':{'Dot_scalar_Gram_combination':'1/6','Codex_scalar_Pfaffian_ancestor':'-189/128','combination':'0'},
       'predecessor_linear_combination':terms,
       'cycle_coordinate':'0',
       'scalar_full_identity_reason':'The scalar Pfaffian ancestor is a degree15 coefficient of the degree11 Pfaffian in SO(11), transformed by det(d/dX)^5, so it is a genuine finite-colour relation of scalar flavour weight(5,5,5). The complete17coordinate frame is injective on this entire relation space. All17 equality checks therefore prove the full scalar Q0(I) identity.',
       'post_return_claim_status':'The separately verified global identity P6 G_C=25200 R_C in ../so11_independent_audit transports this actual scalar predecessor to the entire stable classical derivative of the aligned closed combination. Its classical fortuitous obstruction vanishes. Finite cohomology exactness and a lift for the complete quantum differential are not claimed.',
       'input_certificate_sha256':hashlib.sha256((HERE/'complete_scalar_relation_quotient.json').read_bytes()).hexdigest(),
       'seconds':time.monotonic()-start}
    (HERE/'codex_scalar_preimage.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
