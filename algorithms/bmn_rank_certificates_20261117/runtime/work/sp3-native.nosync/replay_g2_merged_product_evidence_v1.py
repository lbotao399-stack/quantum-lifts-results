"""Replay merged G2 minors, rational probes, and all generic product values.

Source minor determinants are verified but only the saved source probes are
independently reconstructed. Every archived known-product value is regenerated.
"""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import argparse
import json
import pickle
import time

import g2_bmn_complex_v1 as base
import numpy as np
from flint import nmod_mat
from hw_products import Raising
from sp_graviton_singlets import product_values


def replay(evidence):
    started=time.monotonic()
    hashes=json.loads((evidence/'SHA256SUMS.json').read_text())
    for name,digest in hashes.items():assert base.sha(evidence/name)==digest,name
    manifest=json.loads((evidence/'frame_manifest.json').read_text())
    weights,prime=tuple(manifest['weights']),manifest['prime']
    base.setup(weights,tuple(range(1,10)),prime)
    assert json.loads(json.dumps(base.CONTEXT['dimensions']))==manifest['dimensions']
    assert json.loads(json.dumps(base.CONTEXT['frame_metadata']))==manifest['frames']
    records=[]
    for spin in range(1,10):
        payload=pickle.loads((evidence/('s%d.pkl'%spin)).read_bytes())
        assert payload['spin']==spin and payload['prime']==prime
        minor=payload['minor'];rank=minor.shape[0]
        assert minor.shape==(rank,rank) and len(payload['recipes'])==rank
        determinant=int(nmod_mat(minor.tolist(),prime).det())
        assert determinant==payload['determinant'] and determinant
        images=[]
        for word in payload['probe_recipes']:
            expression=base.hw_extremal.project({word:1},(spin,)+weights)
            assert not base.hw_extremal.action(expression,0,1)
            assert not base.hw_extremal.action(expression,1,2)
            image=base.native.sp_census.base.action_expr(expression)
            assert not base.native.sp_census.base.action_expr(image)
            images.append({word:Fraction(value).numerator*pow(Fraction(value).denominator,-1,prime)%prime
                for word,value in image.items() if value})
        coordinates=nonzero=0
        for probe in payload['probes']:
            original=base.CONTEXT['frames'][spin-1][probe['degree']]
            for letter,(odd,matrix) in probe['letters'].items():
                assert original.L[letter][0]==odd
                assert np.array_equal(original.L[letter][1][:probe['points']],matrix)
            frame=base.native.Frame(probe['degree'],probe['points'],1,prime,k=probe['k'])
            frame.L=probe['letters']
            expressions=[{word:value for word,value in image.items()
                if sum(3<=letter<=5 for trace in word[1] for letter in trace)==probe['degree']}
                for image in images]
            actual=base.original_evaluate(expressions,frame)
            assert np.array_equal(actual,probe['values']),(spin,probe['degree'])
            coordinates+=actual.size;nonzero+=int(np.count_nonzero(actual))
        row=dict(spin=spin,rank_lower_bound=rank,minor_determinant=determinant,
            rational_projected_recipes=len(images),probe_coordinates=coordinates,nonzero_probe_coordinates=nonzero)
        records.append(row);print(json.dumps(dict(stage='source_rank_replay',**row)),flush=True)

    ranks={row['spin']:row['rank_lower_bound'] for row in records};ranks.update({0:0,10:0})
    dimensions={spin:sum(parts.values()) for spin,parts in base.CONTEXT['dimensions'].items()}
    bounds=[dict(spin=spin,dimension=dimensions[spin],outgoing_rank_lower_bound=ranks[spin],
        incoming_rank_lower_bound=ranks[spin+1],cohomology_upper_bound=dimensions[spin]-ranks[spin]-ranks[spin+1])
        for spin in range(10)]
    assert [row['cohomology_upper_bound'] for row in bounds]==[0,0,0,0,0,4,0,0,0,0]
    summary=json.loads((evidence/'result.json').read_text())
    for row,original in zip(bounds,summary['bounds']):
        assert all(original[key]==value for key,value in row.items())

    closure=set();raising_rows=0;products=[]
    for spin in (1,3,5,7):
        family=pickle.loads((evidence/('family_s%d.pkl'%spin)).read_bytes())
        payload=pickle.loads((evidence/('products_s%d.pkl'%spin)).read_bytes())
        incoming=pickle.loads((evidence/('s%d.pkl'%(spin+1))).read_bytes())
        assert payload['incoming_recipes']==incoming['recipes']
        kernel=family['kernel']
        for generator in family['generators']:
            fingerprint=repr(sorted(generator['expr'].items()))
            if fingerprint not in closure:
                assert not base.native.sp_census.base.action_expr(generator['expr'])
                closure.add(fingerprint)
        raising=Raising(family['generators'],prime)
        for row_index in range(kernel.shape[0]):
            row=kernel[row_index]
            for root in (1,2):
                residual=defaultdict(int)
                for index,value in zip(row.indices,row.data):
                    for word,coefficient in raising.product(family['products'][int(index)],root).items():
                        residual[word]+=int(value)*coefficient
                assert not any(value%prime for value in residual.values()),(spin,row_index,root)
            raising_rows+=1
        needed=sorted(set(map(int,kernel.indices)))
        assert int(max(np.diff(kernel.indptr),default=0))*(prime-1)**2<2**63
        blocks=[]
        for degree,frame in base.CONTEXT['frames'][spin].items():
            raw=product_values(family['generators'],[family['products'][index] for index in needed],frame,degree)
            blocks.append(np.asarray(kernel[:,needed]@raw)%prime)
        values=np.concatenate(blocks,axis=1)
        assert np.array_equal(values,payload['values']),(spin,'product_values')
        minor=payload['minor']
        determinant=int(nmod_mat(minor.tolist(),prime).det())
        assert determinant==payload['determinant'] and determinant
        incoming_rank=ranks[spin+1]
        selected_products=[index for index,row in enumerate(payload['rows']) if row>=incoming_rank]
        for index in selected_products:
            assert np.array_equal(minor[index],values[payload['rows'][index]-incoming_rank,payload['columns']])
        # Check every overlapping production incoming entry between certificates.
        selected_incoming=[index for index,row in enumerate(payload['rows']) if row<incoming_rank]
        column_map={column:index for index,column in enumerate(incoming['columns'])}
        common=[index for index,column in enumerate(payload['columns']) if column in column_map]
        original_rows=[payload['rows'][index] for index in selected_incoming]
        original_columns=[column_map[payload['columns'][index]] for index in common]
        assert np.array_equal(minor[np.ix_(selected_incoming,common)],
                              incoming['minor'][np.ix_(original_rows,original_columns)])
        incoming_upper=dimensions[spin+1]-ranks.get(spin+2,0)
        lower=max(0,len(minor)-incoming_upper)
        assert incoming_upper==incoming_rank
        assert lower==bounds[spin]['cohomology_upper_bound']
        result=payload['result']
        assert result['combined_rank_lower_bound']==len(minor)
        assert result['known_product_cohomology_lower_bound']==lower
        record=dict(spin=spin,formal_highest_weight_rows=len(values),minor_determinant=determinant,
            product_coordinates_reconstructed=int(values.size),incoming_rank=incoming_rank,
            combined_rank_lower_bound=len(minor),known_product_cohomology_lower_bound=lower,
            cohomology_upper_bound=bounds[spin]['cohomology_upper_bound'],
            matched_incoming_minor_coordinates=len(selected_incoming)*len(common),
            known_products_exhaust_cohomology=True)
        products.append(record);print(json.dumps(dict(stage='generic_product_replay',**record)),flush=True)

    result=dict(status='passed_merged_rank_minors_rational_probes_and_all_generic_product_values',
        group='G2',weights=weights,prime=prime,payload_hashes=len(hashes),
        source_ranks=records,bounds=bounds,products=products,
        exact_closed_generators=len(closure),modular_highest_weight_rows=raising_rows,
        rational_projected_recipes=sum(row['rational_projected_recipes'] for row in records),
        independent_source_probe_coordinates=sum(row['probe_coordinates'] for row in records),
        product_coordinates_reconstructed=sum(row['product_coordinates_reconstructed'] for row in products),
        nonzero_determinants=len(records)+len(products),
        all_source_minor_entries_independently_reconstructed=False,
        all_product_values_reconstructed=True,
        cohomology_dimensions=[row['cohomology_upper_bound'] for row in bounds],
        no_new_fortuitous_family_in_this_column=True,
        d1_increment=0,new_generator_increment=0,d2_increment=0,d3_increment=0,np_increment=0,
        driver_sha256=base.sha(__file__),seconds=time.monotonic()-started,
        scope='All saved determinants, exact character dimensions, 27 rational projected source recipes and their saved dense-matrix Q0 probes; exact generator closure, modular raising equations and every generic product value. Source minors remain production data beyond the independently reconstructed probes. No quantum differential was computed.')
    (evidence/'replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('source_ranks','bounds','products')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--evidence',type=Path,required=True)
    replay(parser.parse_args().evidence)
