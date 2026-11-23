#!/usr/bin/env python3
"""Standard-library verification of every compact scalar preimage coordinate."""
from fractions import Fraction
from pathlib import Path
import json,time

HERE=Path(__file__).resolve().parent

def main():
    start=time.monotonic()
    cert=json.loads((HERE/'complete_scalar_relation_quotient.json').read_text())
    seeds=json.loads((HERE/'scalar_preimages.json').read_text())
    extra=json.loads((HERE/'codex_scalar_preimage.json').read_text())
    assert seeds['terminal_probe_words']==extra['terminal_probe_words']==cert['terminal_probe_words']
    words=cert['terminal_probe_words'];width=len(words)
    images={i:list(map(Fraction,row)) for i,row in zip(cert['source_relation_basis_indices'],cert['Q0_relation_matrix_on_probes'])}
    raw={k:list(map(Fraction,row)) for k,row in seeds['Dot_seed_probe_coefficients'].items()}
    dual=list(map(Fraction,cert['dual_coefficients']))
    dot=[raw['cycle'][i]-5*raw['path3'][i]+Fraction(500,3)*raw['path3-alt'][i] for i in range(width)]
    assert list(map(str,dot))==extra['Dot_scalar_combination_coefficients']
    seedpairs={label:sum(x*y for x,y in zip(dual,row)) for label,row in raw.items()}
    assert seedpairs=={'cycle':Fraction(1),'path3':Fraction(-1,2),'path3-alt':Fraction(-1,50)}
    def combine_predecessors(terms):
        return [sum(Fraction(t['coefficient'])*images[t['source_Gram_index']][i] for t in terms) for i in range(width)]
    count=0
    for record in seeds['relations']:
        rhs=[sum(Fraction(c)*raw[k][i] for k,c in record['terminal_linear_combination'].items()) for i in range(width)]
        lhs=combine_predecessors(record['predecessor_linear_combination'])
        assert lhs==rhs,record['name'];count+=width
    codex=list(map(Fraction,extra['scalar_Pfaffian_ancestor_coefficients']))
    assert sum(x*y for x,y in zip(dual,codex))==Fraction(-189,128)
    wanted=[9*d+Fraction(64,63)*c for d,c in zip(dot,codex)]
    assert list(map(str,wanted))==extra['combined_coefficients']
    assert combine_predecessors(extra['predecessor_linear_combination'])==wanted;count+=width
    assert sum(x*y for x,y in zip(dual,wanted))==0
    out={'passed':True,'implementation':'Python standard-library Fraction; matrix-times-witness checks, no solver used',
         'Dot_preimage_identities':len(seeds['relations']),'Codex_Dot_scalar_preimage_identities':1,
         'exact_coordinate_equalities':count,'Dot_pairing':'1/6','Codex_scalar_pairing':'-189/128','combined_pairing':'0',
         'scope':'Full scalar identities follow from certified injectivity of the17-coordinate relation frame. Post-return claims additionally require the global return identity audit.',
         'seconds':time.monotonic()-start}
    (HERE/'independent_scalar_preimages_replay.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
