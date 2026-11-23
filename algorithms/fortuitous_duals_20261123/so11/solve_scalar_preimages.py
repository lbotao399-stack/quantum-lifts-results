#!/usr/bin/env python3
"""Explicit relation-ideal predecessors for differences of the Dot Gram seeds."""
from pathlib import Path
from fractions import Fraction
import json
from flint import fmpq_mat

HERE=Path(__file__).resolve().parent

def main():
    cert=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    b=fmpq_mat(cert['terminal_basis_probe_matrix']);words=cert['terminal_probe_words']
    raw={}
    for label in ['cycle','path3','path3-alt']:
        r=json.loads((HERE/'terminal_factor5'/f'{label}.json').read_text());den=int(r['denominator'])
        raw[label]={w:Fraction(c)/den for w,c in r['projected']}
    definitions={
      'path3_plus_half_cycle':{'path3':'1','cycle':'1/2'},
      'alt_plus_one_fiftieth_cycle':{'path3-alt':'1','cycle':'1/50'},
      'Dot_combination_minus_one_sixth_cycle':{'cycle':'5/6','path3':'-5','path3-alt':'500/3'}}
    records=[]
    for name,combination in definitions.items():
        wanted=[sum(Fraction(c)*raw[label].get(w,0) for label,c in combination.items()) for w in words]
        x=b.transpose().solve(fmpq_mat([[str(v)] for v in wanted]));assert x[16,0]==0
        terms=[{'source_Gram_index':index,'coefficient':str(x[i,0])} for i,index in enumerate(cert['terminal_basis_Q0_source_indices']) if x[i,0]]
        source={r['source_Gram_index']:json.loads((HERE/'source_factor5'/f"projected_{r['source_Gram_index']}.json").read_text()) for r in terms}
        lhs={}
        for term in terms:
            r=source[term['source_Gram_index']]
            for w,c in r['Q0_projected']:lhs[w]=lhs.get(w,Fraction())+Fraction(term['coefficient'])*Fraction(c)/int(r['denominator'])
        rhs={}
        for label,c in combination.items():
            for w,v in raw[label].items():rhs[w]=rhs.get(w,Fraction())+Fraction(c)*v
        lhs={w:v for w,v in lhs.items() if v};rhs={w:v for w,v in rhs.items() if v};assert lhs==rhs
        records.append({'name':name,'terminal_linear_combination':combination,'predecessor_linear_combination':terms,
                        'recorded_high_factor_coordinates_checked':len(set(lhs)|set(rhs)),
                        'full_identity_reason':'Both sides lie in the complete17dimensional scalar relation space; the certified17coordinate restriction is injective.'})
    out={'passed':True,'terminal_probe_words':words,
         'Dot_seed_probe_coefficients':{label:[str(poly.get(w,0)) for w in words] for label,poly in raw.items()},
         'relations':records,'scope':'Exact full scalar relation-ideal preimages. These identities transport under the declared six odd returns. No identification with the separately specified Codex source is used.'}
    (HERE/'scalar_preimages.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
