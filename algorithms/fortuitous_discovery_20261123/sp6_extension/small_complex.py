"""Exact small-complex obstruction at Sp(6), without a large stable census.

A common-kernel defect is only a candidate until a universal Q0 preimage
has been independently constructed and checked.
"""
from fast_paths import *
import argparse
from fractions import Fraction
import hashlib
import json
import math
import pickle
import time
import abelian_obstruction as ab
from flint import nmod_mat,fmpq_mat,fmpz_mat


def dump(path,obj):
    path.write_bytes(pickle.dumps(obj,protocol=4))


def main():
    cache_native_canonicalization()
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,default=6)
    args=ap.parse_args();rank=args.rank
    folder=HERE/f'frames_sp{rank}';output=HERE/f'small_sp{rank}';output.mkdir(exist_ok=True)
    start=time.monotonic();lengths=(2*rank+2,2*rank+3,2*rank+4)
    frames={l:pickle.loads((folder/f'relations_L{l}.pkl').read_bytes())[0] for l in lengths}
    assert list(map(len,frames.values()))==[67,188,183]
    def qrows(length):
        cache=folder/f'Q{length}';cache.mkdir(exist_ok=True);result=[]
        for i,row in enumerate(frames[length]):
            path=cache/f'row_{i}.pkl';digest=hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest()
            if path.exists():
                old,image=pickle.loads(path.read_bytes());assert old==digest
            else:
                image=pf.m.q(row);dump(path,(digest,image))
            result.append(image)
            if i%30==0:print('Q',rank,length,i,'seconds',time.monotonic()-start,flush=True)
        return result
    qbefore=qrows(lengths[0]);qsource=qrows(lengths[1]);qtarget=qrows(lengths[2])
    before_cert=ab.rank_certificate(qbefore,1000003)
    source_cert=ab.rank_certificate(qsource,1000003)
    assert before_cert['rank']==67 and source_cert['rank']==121
    # This is an exact identity of the native differential, checked on every
    # admissible letter in these Taylor and theta bounds.
    from fortuity.full.letters import letters_up_to,Q_letter
    for letter in letters_up_to(3,3):
        twice=Counter()
        for c,left,right in Q_letter(letter):
            for d,a,b in Q_letter(left):twice[a,b,right]+=c*d
            for d,a,b in Q_letter(right):twice[left,a,b]+=c*d*(-1)**pf.is_ferm(left)
        assert not any(twice.values()),('free-associative Q^2',letter)
    # Check the inclusion of 67 independent first images in the 188 source
    # identities over Q, so the source differential has rank at most 121.
    source_frame=frames[lengths[1]]
    frame_cert=ab.rank_certificate(source_frame,1000003)
    inclusion_coordinates=ab.rational_rank_upper(source_frame+qbefore,frame_cert)[188:]
    relations=frames[lengths[2]]
    cartan=[ab.restrict(row) for row in relations]
    combined=[{**{('Q0',w):c for w,c in q.items()},**{('Cartan',w):c for w,c in a.items()}}
              for q,a in zip(qtarget,cartan)]
    qcert=ab.rank_certificate(qtarget,1000003)
    acert=ab.rank_certificate(cartan,1000003)
    ccert=ab.rank_certificate(combined,1000003)
    print('ranks',qcert['rank'],acert['rank'],ccert['rank'],'seconds',time.monotonic()-start,flush=True)
    coordinates=ab.rational_rank_upper(combined,ccert)
    qcoordinates=ab.rational_rank_upper(qtarget,qcert)
    acoordinates=ab.rational_rank_upper(cartan,acert)
    selected=ccert['rows'];columns=ccert['columns']
    determinant=int(fmpz_mat([[combined[i].get(w,0) for w in columns] for i in selected]).det())
    relation_cert=ab.rank_certificate(relations,1000003)
    assert relation_cert['rank']==183
    prime=2147483647;probes=relation_cert['columns']
    matrix=nmod_mat([[row.get(w,0)%prime for row in relations] for w in probes],prime)
    rhs=nmod_mat([[row.get(w,0)%prime for row in qsource] for w in probes],prime)
    boundary_coordinates=matrix.solve(rhs).transpose()
    rr,boundary_rank=boundary_coordinates.rref();assert boundary_rank==121
    pivots=[];j=0
    for i in range(boundary_rank):
        while not rr[i,j]:j+=1
        pivots.append(j);j+=1
    candidates=[]
    for i in range(183):
        if i in selected:continue
        coeff=[Fraction(0) for _ in range(183)];coeff[i]=1
        for j,c in zip(selected,coordinates[i]):coeff[j]-=Fraction(c)
        denominator=math.lcm(*(c.denominator for c in coeff))
        ints=[int(c*denominator) for c in coeff]
        divisor=math.gcd(*ints);ints=[c//divisor for c in ints]
        mod=[c%prime for c in ints]
        for row,pivot in enumerate(pivots):
            c=mod[pivot]
            if c:mod=[(v-c*int(rr[row,j]))%prime for j,v in enumerate(mod)]
        if any(mod):
            estimate=sum(len(row) for row,c in zip(relations,ints) if c)
            candidates.append((estimate,i,ints))
    candidates.sort();assert candidates
    polynomial=Counter();chosen=candidates[0]
    for row,c in zip(relations,chosen[2]):
        if c:
            for w,d in row.items():polynomial[w]+=c*d
    polynomial=clean(polynomial)
    assert not pf.m.q(polynomial) and not ab.restrict(polynomial)
    dump(output/'candidate_defect.pkl',polynomial)
    dump(output/'candidate_relation_coefficients.pkl',chosen[2])
    # Non-boundary certificate using ONLY the connected, single-trace part.
    single=[{w:c for w,c in row.items() if len(w)==1} for row in qsource]
    single.append({w:c for w,c in polynomial.items() if len(w)==1})
    dualcert=ab.rank_certificate(single,1000003)
    assert dualcert['rank']==122 and 188 in dualcert['rows']
    dualmatrix=fmpq_mat([[single[i].get(w,0) for w in dualcert['columns']] for i in dualcert['rows']])
    dualrhs=fmpq_mat([[int(i==188)] for i in dualcert['rows']])
    dual=dualmatrix.solve(dualrhs)
    values=[Fraction(str(c)) for c in dual.entries()]
    for i,row in enumerate(single):
        assert sum(c*row.get(w,0) for w,c in zip(dualcert['columns'],values))==int(i==188)
    result=dict(group=f'Sp({rank})',native_weight=[3,3,2*rank-1,0,0],
        relation_dimensions=[67,188,183],incoming_relation_rank=67,source_relation_Q_rank=121,
        target_Q_rank=qcert['rank'],stable_abelianization_rank=acert['rank'],
        combined_rank=ccert['rank'],combined_exact_minor_determinant=str(determinant),
        combined_minor_rows=selected,
        combined_minor_columns=[(kind,pf.m.encode(w)) for kind,w in columns],
        common_kernel_dimension=183-ccert['rank'],
        common_kernel_mod_boundaries_dimension=183-ccert['rank']-121,
        all_183_combined_rows_reconstructed_over_Q=True,
        free_associative_Q2_letter_checks=len(letters_up_to(3,3)),
        all_67_preceding_images_in_source_reconstructed_over_Q=True,
        candidate_terms=len(polynomial),candidate_relation_index=chosen[1],
        candidate_relation_coefficients=chosen[2],candidate_exact_Q0_closed=True,
        candidate_stable_abelianization_zero=True,candidate_exact_nonboundary=True,
        single_trace_boundary_rank=121,single_trace_augmented_rank=122,
        dual_support=[(str(c),pf.m.encode(w)) for w,c in zip(dualcert['columns'],values) if c],
        dual_field='Q',dual_pairing='1',universal_preimage_NOT_YET_proved=True,
        fortuitous_class_NOT_YET_claimed=True,seconds=time.monotonic()-start)
    (output/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'rational_coordinates.json').write_text(json.dumps(dict(
        combined_basis=selected,combined_coordinates=coordinates,
        Q0_basis=qcert['rows'],Q0_coordinates=qcoordinates,
        abelianization_basis=acert['rows'],abelianization_coordinates=acoordinates,
        preceding_image_source_coordinates=inclusion_coordinates),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('dual_support','candidate_relation_coefficients')}),flush=True)


if __name__=='__main__':main()
