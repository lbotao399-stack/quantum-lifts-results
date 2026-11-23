#!/usr/bin/env python3
"""Character-complete rational relation quotient in the scalar (5,5,5) sector."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,sys,time
from flint import fmpq_mat,nmod_mat

HERE=Path(__file__).resolve().parent
PRIME=1000003

def pivots(rows):
    a=nmod_mat(rows,PRIME);r,rank=a.rref()
    return [next(j for j in range(r.ncols()) if r[i,j]) for i in range(rank)]

def minor(expressions):
    words=sorted(set().union(*(e for e in expressions)))
    rows=[[int(e.get(w,0))%PRIME for w in words] for e in expressions]
    cols=pivots(rows);assert len(cols)==len(expressions)
    selected=[words[j] for j in cols]
    square=[[e.get(w,0) for w in selected] for e in expressions]
    det=int(nmod_mat([[int(v)%PRIME for v in row] for row in square],PRIME).det());assert det
    return selected,square,det

def rational_rows(exprs,words,den):
    return [[str(Fraction(e.get(w,0),den)) for w in words] for e in exprs]

def main():
    start=time.monotonic();folder=HERE/'source_factor5'
    history=json.loads((folder/'rank.json').read_text())
    selected=[r['index'] for r in history['rows'] if r['independent']]
    assert len(selected)==21
    source=[];images=[];receipts=[]
    for i in selected:
        path=folder/f'projected_{i}.json';r=json.loads(path.read_text())
        den=int(r['denominator'])
        source.append({w:int(c) for w,c in r['projected']})
        images.append({w:int(c) for w,c in r['Q0_projected']})
        receipts.append({'index':i,'projected_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                         'source_spec':history['rows'][i]['source']})
    qindices=[r['index'] for r in history['rows'] if r['image_independent']]
    assert len(qindices)==16
    target=[images[selected.index(i)] for i in qindices]
    dot={}
    for name in ['cycle','path3','path3-alt']:
        r=json.loads((HERE/'terminal_factor5'/f'{name}.json').read_text())
        assert int(r['denominator'])==den
        dot[name]={w:int(c) for w,c in r['projected']}
    target.append(dot['cycle'])
    swords,smat,sdet=minor(source)
    twords,tmat,tdet=minor(target)
    smat=rational_rows(source,swords,den);tmat=rational_rows(target,twords,den)
    image_rows=rational_rows(images,twords,den)
    b=fmpq_mat(tmat);a=fmpq_mat(image_rows)
    assert b.rank()==17 and a.rank()==16
    dual=b.solve(fmpq_mat([[0]]*16+[[1]]))
    assert all(v==0 for row in (a*dual).tolist() for v in row)
    dotscores={name:str((fmpq_mat(rational_rows([expr],twords,den))*dual)[0,0]) for name,expr in dot.items()}
    combined=Fraction(dotscores['cycle'])-5*Fraction(dotscores['path3'])+Fraction(500,3)*Fraction(dotscores['path3-alt'])
    # Verify the all-coefficient recorded projected Q0 maps, not only minors.
    allwords=sorted(set().union(*target))
    for e in images:
        coords=b.transpose().solve(fmpq_mat([[str(Fraction(e.get(w,0),den))] for w in twords]))
        assert coords[16,0]==0
        for word in allwords:
            assert sum(coords[k,0]*target[k].get(word,0) for k in range(17))==e.get(word,0)
        assert not(set(e)-set(allwords))
    record={'passed':True,'field_for_basis_minors':PRIME,'scalar_highest_weight':[5,5,5],
       'Taylor_diagram':[3,3,3,3,3],'source_length':14,'terminal_length':15,
       'character_source_relation_dimension':21,'character_terminal_relation_dimension':17,
       'source_relation_basis_indices':selected,'source_basis_probe_words':swords,
       'source_basis_probe_matrix':smat,'source_minor_mod_prime':sdet,
       'terminal_basis_Q0_source_indices':qindices,'terminal_basis_last':'projected Dot cycle Gram identity',
       'terminal_probe_words':twords,'terminal_basis_probe_matrix':tmat,'terminal_minor_mod_prime':tdet,
       'Q0_relation_matrix_on_probes':image_rows,'rational_Q0_rank':16,'relation_quotient_dimension':1,
       'dual_coefficients':[str(x[0]) for x in dual.tolist()],
       'dual_pairing_on_terminal_basis':['0']*16+['1'],
       'dual_pairing_on_all21_Q0_relations':['0']*21,'dot_seed_pairings':dotscores,
       'dot_combination_pairing':str(combined),'all_recorded_high_factor_Q0_coordinates_checked':len(images)*len(allwords),
       'proof_scope':'Complete rational scalar-(5,5,5) BMN relation quotient, using exact orthogonal character dimensions and Gram identities. Trace factors>=5 are used as an injective coordinate frame, certified by17independent actual terminal relations. All21predecessor relations form a complete basis. No physical full-singlet conclusion after six odd returns without separate Taylor-component membership proof.',
       'input_receipts':receipts,'seconds':time.monotonic()-start}
    (HERE/'complete_scalar_relation_quotient.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ['input_receipts','source_basis_probe_matrix','terminal_basis_probe_matrix','Q0_relation_matrix_on_probes','source_basis_probe_words','terminal_probe_words']},indent=2))
if __name__=='__main__':main()
