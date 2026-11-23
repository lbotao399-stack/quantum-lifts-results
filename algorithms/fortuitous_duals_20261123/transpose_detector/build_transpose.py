"""Transpose the certified40-term relation detector into sparse source coordinates.

This is an exact coefficient functional on the fixed highest-weight sector.
No Q0 image or harmonic representative is constructed when scoring a new
candidate. The source functional is obtained by enumerating all possible
adjacent inverse splittings of the40 supported target trace words.
"""
from pathlib import Path
from fractions import Fraction
import sys,pickle,json,hashlib,time,gzip,statistics,gc
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'positive_sector'))
from probe_stable import m

def combine(a,b):
    if set(a[2]).intersection(b[2]):return None
    return a[0]+b[0],a[1]+b[1],tuple(sorted(a[2]+b[2]))

def score(functional,polynomial):
    return sum(c*polynomial.get(w,0) for w,c in functional.items())

def main():
    started=time.monotonic()
    positive=HERE.parent/'positive_sector'
    _,centre,incoming,outgoing=pickle.loads((positive/'stable_rows.pkl').read_bytes())
    native_words=sorted({w for row in outgoing for w in row})
    saved=json.loads((HERE.parent/'enumeration/sp3_intrinsic_dual.json').read_text())
    target={m.decode(row['ordered_native_traces']):row['integer_coefficient']
            for row in saved['detector']}
    assert all(native_words[row['index']]==m.decode(row['ordered_native_traces'])
               for row in saved['detector'])
    assert len(target)==40
    candidates=set()
    for word in target:
        for ti,trace in enumerate(word):
            for i in range(len(trace)):
                rotated=trace[i:]+trace[:i]
                new_letter=combine(rotated[0],rotated[1])
                if new_letter is None:continue
                predecessor=word[:ti]+((new_letter,)+rotated[2:],)+word[ti+1:]
                sign,canonical=m.so_canon(predecessor)
                if sign:candidates.add(canonical)
    source={w:score(target,m.q({w:1})) for w in sorted(candidates)}
    source={w:c for w,c in source.items() if c}
    # Exhaustive action on the complete formal highest-weight source basis.
    all_scores=[score(source,row) for row in centre]
    reference=[score(target,row) for row in outgoing]
    assert all_scores==reference
    assert all(score(source,row)==0 for row in incoming)
    omega_data=json.loads(gzip.decompress((positive/'omega_rational.json.gz').read_bytes()))
    omega={m.decode(w):Fraction(c) for c,w in omega_data['terms']}
    defect=m.q(omega)
    assert score(source,omega)==score(target,defect)==1
    relation_data=json.loads(gzip.decompress((positive/'rational_relation_basis.json.gz').read_bytes()))
    for row in relation_data['basis_coefficients']:
        assert sum(Fraction(c)*all_scores[i] for i,c in row)==0
    # Recheck the exact arithmetic proof gate rather than treating an
    # unvalidated rational reconstruction as an identity.
    gate=json.loads((positive/'rational_relations_result.json').read_text())
    assert gate['passed'] and gate['crt_modulus']>2*gate['integer_evaluation_absolute_bound']
    assert gate['rational_relation_dimension']==181 and gate['full_finite_frame_minor_determinant']
    times=[]
    # Paired warm-cache timings: both ways score the exact same rational
    # input, while the full-image method benefits from all cached letter Qs.
    assert score(target,m.q(omega))==score(source,omega)==1
    for trial in range(7):
        values={}
        for label in (['full_image','transpose'] if trial%2==0 else ['transpose','full_image']):
            gc.collect();t=time.perf_counter()
            answer=score(target,m.q(omega)) if label=='full_image' else score(source,omega)
            elapsed=time.perf_counter()-t;assert answer==1
            values[label]=elapsed
        values['ratio']=values['full_image']/values['transpose'];times.append(values)
    payload={'weight':[3,3,5,0,0],'length':9,'source_terms':[[m.encode(w),str(c)] for w,c in source.items()],
             'target_terms':[[m.encode(w),str(c)] for w,c in target.items()]}
    data=(json.dumps(payload,separators=(',',':'))+'\n').encode()
    (HERE/'omega_source_dual.json').write_bytes(data)
    result={'passed':True,'field':'Q','inverse_splitting_candidates':len(candidates),
            'target_support':len(target),'source_support':len(source),
            'source_terms_present_in_omega':sum(w in omega for w in source),
            'input_omega_terms':len(omega),'full_Q0_omega_terms':len(defect),
            'formal_source_basis_equations':len(centre),'formal_boundary_equations':len(incoming),
            'rational_identity_equations':len(relation_data['basis_coefficients']),
            'omega_pairing':'1','coefficients_are_integers':all(isinstance(x,int) for x in source.values()),
            'max_abs_source_coefficient':max(map(abs,source.values())),
            'paired_warm_cache_timings_seconds':times,
            'explicit_warmup_both_methods':True,
            'median_pair_ratio':statistics.median(t['ratio'] for t in times),
            'source_dual_sha256':hashlib.sha256(data).hexdigest(),
            'rational_relation_basis_sha256':gate['rational_basis_sha256'],
            'seconds':time.monotonic()-started,
            'scope':'The functional descends and detects the full fortuitous quotient only in the complete fixed highest-weight sector certified above. The inverse-splitting transpose is exact on all raw cochains of this weight, but finite-colour descent outside the selected highest-weight sector has not been certified. Timing measures only repeated candidate scoring after the fixed detector is built, not end-to-end census speed.'}
    (HERE/'transpose_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
