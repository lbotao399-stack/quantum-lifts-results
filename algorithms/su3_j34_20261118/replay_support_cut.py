"""Rebuild a support-cut correction, check the three-colour identities, and classify all 500 contents.

Does not use the producer's basis projection, frame matrices, solver, or pivots.
The two homogeneous equations are checked separately from the original 70.
"""
import argparse
import hashlib

from collections import Counter
from fractions import Fraction
from pathlib import Path
import full_content_probes as probe
import parent_frontier
from su_trace_adapter_v1 import parity_project

r=probe.cen.r
s=r.s
HERE=Path(__file__).resolve().parent
BASE=HERE/'joint_70'
VARIANTS=('preserve_296','all_nonincident')

def ch(w):
    return bytes(sorted(b''.join(w))).hex()

def checked_vector(path):
    doc=s.read(path)
    if isinstance(doc,dict):
        raw=doc['words']; den=int(doc['denominator']); assert den>0
        values={s.unpack_hex(traces):Fraction(int(c),den) for traces,c in raw}
    else:
        raw=doc; values={s.unpack_hex(traces):Fraction(c) for traces,c in raw}
    assert len(values)==len(raw) and all(values.values())
    return values

def rebuild(variant,files,receipt):
    assert variant in VARIANTS and receipt['job']=='5665560'
    for key,p in [('vector_shared.json.gz',BASE/'vector_shared.json.gz'),('census_2724.json',BASE/'census_2724.json'),('negative_Q1.jsonl.gz',HERE/'negative_Q1.jsonl.gz'),('a8_5664529.json',BASE/'a8_5664529.json')]:
        assert hashlib.sha256(p.read_bytes()).hexdigest()==receipt['inputs'][key]
    for name in ('setup.json',variant+'.certificate.json.gz',variant+'.vector.json.gz',variant+'.audit.json'):
        payload=(files/name).read_bytes(); assert [len(payload),hashlib.sha256(payload).hexdigest()]==receipt['files'][name]
    setup=s.read(files/'setup.json');cert=s.read(files/(variant+'.certificate.json.gz')); original=checked_vector(BASE/'vector_shared.json.gz')
    old=set(s.read(BASE/'census_2724.json')['exact_zero']); targets=set(setup['targets']); fail=set(setup['homogeneous_rhs_word_counts']); assert len(old)==308 and len(targets)==12 and targets<=old and len(fail)==2 and not fail&old
    F,B=parent_frontier.incidence(); fwd={g.hex():{t.hex() for t in ts} for g,ts in F.items()}; parents=lambda ts:{g.hex() for t in ts for g in B[bytes.fromhex(t)]}
    cut=parents(fail); reachable=parents(targets); restricted={g for g in reachable if not ((fwd[g]&old)-targets)}-cut
    assert sorted(cut)==setup['forbidden']==cert['forbidden'] and sorted(restricted)==setup['additional_reset_groups']
    assert len(cut)==19 and len(restricted)==11 and len(reachable-cut)==32
    allowed=restricted if variant==VARIANTS[0] else reachable-cut
    assert sorted(allowed)==cert['allowed'] and cert['kind']=='exact_solution'
    assert dict(Counter(ch(w) for w in original if ch(w) in cut))==setup['removed_words_from_2724']
    x={w:c for w,c in original.items() if ch(w) not in cut|restricted}; assert len(x)==setup['starting_words']==682
    seen=set()
    for g,word,c in cert['frame_solution']:
        w=s.unpack_hex(word); value=Fraction(c);assert ch(w)==g and g in allowed and w not in seen and value;seen.add(w)
        for u,m in parity_project({s.expanded(w):1},1).items():
            key=s.packed(u); x[key]=x.get(key,0)+value*m
            if not x[key]:del x[key]
    assert x==checked_vector(files/(variant+'.vector.json.gz')) and len(x)==cert['vector_words']
    groups={ch(w) for w in set(x)|set(original) if x.get(w,0)!=original.get(w,0)}; affected=set.union(*(fwd[g] for g in groups)); preserved=old-affected
    assert all(ch(w) not in cut for w in x) and not any(setup['homogeneous_rhs_word_counts'].values())
    if variant==VARIANTS[0]:assert len(preserved)>=296 and old-targets<=preserved
    return dict(words=len(x),frame_terms=len(seen),changed_parents=sorted(groups),previous_zeros_untouched_by_support=sorted(preserved)),targets,fail,old

def run(variant,files,receipt_file,columns,exterior,work,out,jobs):
    work.mkdir(exist_ok=True,parents=True);out.mkdir(exist_ok=True,parents=True)
    receipt=s.read(receipt_file);rec=receipt.get('receipt',receipt); result,targets,fail,old=rebuild(variant,files,rec)
    print(variant,'reconstructed',result['words'],result['frame_terms'],flush=True)
    vec=files/(variant+'.vector.json.gz'); rhs=HERE/'negative_Q1.jsonl.gz'; gold=s.read(files/(variant+'.audit.json'))
    selected=r.domain('joint'); must=selected if variant==VARIANTS[0] else {bytes.fromhex(t) for t in targets}; extra={bytes.fromhex(t) for t in fail}
    raw=r.replay(rhs,vec,'joint',columns,work,jobs,extra | (selected-must),must); lookup={a['content']:a for a in raw['contents']+raw['holdouts']}
    assert len(lookup)==72 and gold['ordered_source_sha256']==raw['ordered_source_sha256'];assert {a['content'] for a in gold['contents']}==targets and {a['content'] for a in gold['holdouts']}==fail
    assert all(lookup[a['content']]==a and not a['nonzero_monomials'] for a in gold['contents']+gold['holdouts'])
    assert all(not lookup[t]['rhs_terms'] and not lookup[t]['nonzero_parent_word_incidences'] for t in fail)
    s.write(out/(variant+'.replay72.json'),raw)
    # The census API expects 70 distinct polynomial checks. Keep both other exact identities as independent evidence.
    audit=dict(raw,holdouts=[a for a in raw['holdouts'] if a['content'] not in fail]); auditfile=work/(variant+'.audit70.json');s.write(auditfile,audit)
    print(variant,'exact polynomial checks',len(must),len(raw['holdouts']),flush=True)
    census=probe.run(rhs,vec,auditfile,exterior,out/(variant+'.census.json'),work,[817234,419071]); zeros=set(census['exact_zero']); bad=set(census['mathematically_nonzero']); missing=set(census['unresolved'])
    # The support check can additionally transport exact identities from the original fully classified vector.
    independent_zeros=set(result['previous_zeros_untouched_by_support']) | fail;assert not independent_zeros&bad
    newly_proved=independent_zeros-zeros; assert newly_proved<=missing
    if newly_proved:
        zeros|=newly_proved;missing-=newly_proved;census.update(exact_zero=sorted(zeros),unresolved=sorted(missing),extra_support_or_homogeneous_zero=sorted(newly_proved));s.write(out/(variant+'.census.json'),census)
    assert len(zeros|bad|missing)==500 and len(zeros)+len(bad)+len(missing)==500
    result.update(selected_passes=sum(lookup[t.hex()]['nonzero_monomials']==0 for t in selected),zero=len(zeros),nonzero=len(bad),unresolved=sorted(missing),gained=sorted(zeros-old),lost_to_nonzero=sorted(old&bad),previously_zero_now_unresolved=sorted(old&missing),source_sha256=census['vector_sha']);s.write(out/(variant+'.comparison.json'),result)
    print(variant,'verified',result['selected_passes'],len(zeros),len(bad),len(missing),'gains/losses',len(zeros-old),len(old&bad),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('variant',choices=VARIANTS); ap.add_argument('--files',type=Path,default=BASE/'support_cut');ap.add_argument('--receipt-file',type=Path,default=BASE/'support_cut'/'a8_5665560.json');ap.add_argument('--columns',type=Path,required=True);ap.add_argument('--exterior',type=Path,required=True);ap.add_argument('--work',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--jobs',type=int,default=2);assert __debug__;run(**vars(ap.parse_args()))
