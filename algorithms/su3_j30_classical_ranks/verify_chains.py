"""Bind the four R-007 reductions to their actual matrices and coordinate frames.

Offline, standard library only. This checks original stage records, not a new
elimination of the production arrays. See replay_su3_j30_packed_residuals.py
for the independent multiplication and su3_j30_frame_audit for the minors.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from math import prod
from pathlib import Path
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / 'data'
SCIENCE = ROOT / 'algorithms/su3_independent_prime_20261118'
RUNTIME = SCIENCE / 'runtime/work/sp3-native.nosync'
P = 1000033
h = lambda b: sha256(b).hexdigest()

def pinned(p, digest):
    raw = p.read_bytes()
    assert h(raw) == digest, p
    return json.loads(raw) if p.suffix == '.json' else p

def stage_pivots(rec, seen, n):
    piv = rec['pivots']
    assert rec['rank'] == len(piv) == len(set(piv))
    assert all(0 <= i < n and i not in seen for i in piv)
    seen.update(piv)

def run():
    assert __debug__, 'Enable assertions'
    prod_manifest = pinned(SCIENCE/'SOURCE_MANIFEST.json', '3cd48eea41430f6cd08449f651f729db8125417ce04b1d6b8dbd6f4b77867149')
    plans = {}
    for tag, filename in [('first','su3_correction_profiles_v1'),('second_Q1','su3_second_quantum_profiles_v1'),('second_Q0','su3_second_correction_profiles_v1')]:
        f = RUNTIME/(filename+'.json'); plans[tag] = pinned(f,prod_manifest[f.relative_to(SCIENCE).as_posix()])
    assert plans['second_Q1']['targets'] == plans['second_Q0']['targets']
    source = pinned(DATA/'N4_SU3_j30_source_order_20261119.json','e7c39b5545afd3c57a2104a6e451d24a2908127aa3c3136681620983564a26a4')
    profiles = {}; offset = 0
    for eq,kind in [(1,'first'),(2,'second_Q0')]:
        for i,s in enumerate(plans[kind]['sources']):
            row = source['profiles'][len(profiles)]; profiles[eq,i] = row
            assert (row['equation'],row['index'],row['start'],row['stop'],row['dimension'],row['profile']) == (eq,i,offset,offset+s['dimension'],s['dimension'],s['profile'])
            path = ('first_basis/profile%02d' if eq==1 else 'second_basis/shard0/sources%02d') % i + '/independent_seeds.pkl'
            assert h((SCIENCE/path).read_bytes()) == prod_manifest[path] == row['seeds_sha256']
            offset=row['stop']
        assert offset == (28311 if eq==1 else 77639)
    assert len(profiles)==len(source['profiles'])==74
    independent=pinned(DATA/'N4_SU3_j30_independent_residuals_20261118.json','ebac8b013044e79a3b7bf246555b2c68d3207f2798a580a8f0059b4c0845beb9')
    inc=pinned(DATA/'N4_SU3_j30_incoming_rank_20261119.json','1d174399cd6b0cb8eb606c410368e548cdf0ce2dc54957d92464ae86cde7b686')
    joint=pinned(DATA/'N4_SU3_j30_joint_rank_20261119.json','8d9f8f13fc8b5e640d2cb48a269d4ad65a997667adbd6267d4c15b9b3789e677')
    jp=pinned(HERE/'joint_plan.json',joint['plan']); assert jp['source_dimension']==63983
    ar=pinned(DATA/'N4_SU3_j30_original_stages_p1000033.zip','d46e9a0091c1b4fd423921e4ae22c71fa04e61b652eec7f2ce3b37aa119a87a0')
    fr=pinned(DATA/'N4_SU3_j30_target_frames_p1000033.zip','648bcb49a40980fb98291b940717af4d999428a0fca1424469ff0afc2cafae50')
    rr=pinned(DATA/'N4_SU3_j30_solution_receipts_20261118.zip','b4b4e1df18bfa3666d5c0f9b09443d06e0087a41a78fcb742078f8e378a9557e')
    assert inc['prime']==joint['prime']==independent['prime']==P and inc['matrix_receipt_sha256']==h((DATA/'N4_SU3_j30_independent_residuals_20261118.json').read_bytes())
    replay={(r['equation'],r['target']):r for r in independent['targets']}; inp={}; visited=set(); piv=set(); previous=[]; ranks=[]; counts=Counter(); target_sizes=Counter(); widths=Counter(); digest_links=0; entries=0
    with ZipFile(ar) as z, ZipFile(fr) as frames, ZipFile(rr) as reduced:
        for eq,dirname in [(1,'first'),(2,'second')]:
            tags=['first'] if eq==1 else ['second_Q0','second_Q1']; plan=plans[tags[0]]; n=28311 if eq==1 else 77639
            stages=json.loads(reduced.read(dirname+'_solve_compact/stages.json')); result=json.loads(reduced.read(dirname+'_solve_compact/result.json')); assert len(stages)==len(plan['targets']) and [s['folder'][-2:] for s in stages]==[f'{i:02d}' for i in result['processed_targets']]
            for step,short in enumerate(stages):
                t=int(short['folder'][-2:]); target=plan['targets'][t]; assert target['index']==t
                path=f'{dirname}/stage{step:02d}_target{t:02d}.json'; raw=z.read(path); visited.add(path); rec=json.loads(raw); m=rec['manifest']; r=replay.pop((eq,t)); assert m['previous_stages']==previous; previous.append(h(raw)); assert previous[-1]==short['result_sha256']
                assert all(rec[k]==short[k] for k in ('pivots','rank','rows_sha256','columns_sha256','consistent')) and rec['consistent'] is True
                stage_pivots(rec,piv,n)
                assert (m['prime'],m['target'],m['variables'],m['rhs_sha256'],m['columns'])==(P,t,n,r['rhs_sha256'],r['coordinates']) and r['joint_residual_nonzero']==0
                if eq==1:assert r['standalone_first_residual_nonzero']==0 and m['plan_sha256']==prod_manifest['runtime/work/sp3-native.nosync/su3_correction_profiles_v1.json']
                else:
                    assert (m['first_variables'],m['second_variables'])==(28311,49328)
                    assert all(m[k]==prod_manifest['runtime/work/sp3-native.nosync/'+fname+'.json'] for k,fname in [('first_plan_sha256','su3_correction_profiles_v1'),('second_plan_sha256','su3_second_correction_profiles_v1'),('quantum_plan_sha256','su3_second_quantum_profiles_v1')])
                base=f'e{eq}_t{t:02d}/'; fmraw=frames.read(base+'manifest.json');fm=json.loads(fmraw);f=json.loads(frames.read(base+'result.json')); w=f['witness'];dim=target['dimension'];width=r['coordinates'];
                assert (fm['target'],fm['equation'],fm['prime'],fm['frame'],fm['rhs_sha256'])==(target,eq,P,m['rhs_frame'],r['rhs_sha256']) and fm['frame']['columns']==width
                assert (f['manifest_sha256'],f['equation'],f['target'],f['prime'],f['width'],f['dimension'],f['injective'],w['rank'])==(h(fmraw),eq,t,P,width,dim,True,dim)
                assert len(set(map(tuple,w['rows'])))==len(w['rows'])==len(set(w['columns']))==len(w['columns'])==dim and all(0<=i<width for i in w['columns'])
                assert sum(k for k,d in w['panel_orders_and_determinants'])==dim and all(k>0 and 0<d<P for k,d in w['panel_orders_and_determinants'])
                det=1
                for _,d in w['panel_orders_and_determinants']:det=det*d%P
                assert det==w['determinant_in_selected_order']!=0
                widths[eq]+=width;target_sizes[eq]+=dim
                if eq==2:
                    js=jp['sources'][t]; assert js['recipe_sha256']==fm['recipe_sha256'] and all(js[k]==target[k] for k in target) and f['seed_rank']==dim
                matrices=r['matrices']; blocks=m['blocks'];assert len(matrices)==len(blocks) and set(x[0] for x in matrices)<=set(tags)
                assert len(blocks)==sum(len([s for s in plans[tag]['sources'] if t in s['targets']]) for tag in tags)
                if eq==2:assert [b['kind'] for b in blocks]==sorted([b['kind'] for b in blocks],reverse=True) # original Q1 then Q0
                for tag in tags:
                    expected=[s['index'] for s in plans[tag]['sources'] if t in s['targets']]; rows=[x for x in matrices if x[0]==tag]; bs=[b for b in blocks if eq==1 or b['kind']==tag[-2:]]
                    assert [x[1] for x in rows]==[b['source'] for b in bs]==expected
                    for x,b in zip(rows,bs,strict=True):
                        i=b['source'];s=profiles[2 if tag=='second_Q0' else 1,i];assert x==[tag,i,t,s['dimension'],width,b['sha256']] and (b['start'],b['stop'])==(0,s['dimension'])
                        if eq==2:assert (b['global_start'],b['source_basis_sha256'])==(s['start'],s['seeds_sha256'])
                        counts[tag]+=1;entries+=s['dimension']*width
                    if tag=='second_Q0':inp[t]=[b['sha256'] for b in bs]
                extras={b['source']:b for b in f['packed_Q0_sources']}
                for i,b in extras.items():assert eq==1 and [x[-1] for x in matrices if x[1]==i]==[b['sha256']]
                for label in w['rows']:
                    if label[0]=='seed':assert len(label)==2 and 0<=label[1]<fm['recipe_count']
                    else:assert eq==1 and label[0]=='image' and len(label)==3 and 0<=label[2]<extras[label[1]]['rows_used_through']<=profiles[1,label[1]]['dimension']
            assert previous==result['stage_certificates'] and len(piv)==result['equation_rank'];ranks.append(len(piv))
        assert not replay and len(frames.namelist())==100
        seen=set(); stamps=[]
        for s in inc['stages']:
            t=s['target'];path=f'incoming/stage{t:02d}.json';raw=z.read(path);visited.add(path);r=json.loads(raw);g=r['guard'];assert g['previous']==stamps;stamps.append(h(raw));assert (stamps[-1],g['target'],g['prime'],g['input_digests'],len(g['input_digests']),r['rank'])==(s['stage_sha256'],t,P,inp.pop(t),s['matrices'],s['rank_increment']); assert (g['plan_sha256'],g['driver_sha256'])==(inc['input_plan_sha256'],inc['driver_sha256']);digest_links+=s['matrices'];stage_pivots(r,seen,49328)
        assert not inp and len(seen)==inc['incoming_rank_over_prime'];incoming=len(seen)
        seen=set();stamps=[];outgoing=0
        for step,s in enumerate(joint['stages']):
            key=s['stage'];path='joint/'+key+'.json';raw=z.read(path);visited.add(path);r=json.loads(raw);g=r['guard'];assert g['previous']==stamps;stamps.append(h(raw));assert (stamps[-1],g['key'],g['width'],r['rank'],g['inputs'],g['plan'])==(s['stage_sha256'],key,s['coordinates'],s['rank_increment'],joint['physical_inputs'],joint['plan']);assert (step<19)==key.startswith('q');assert g['driver']==h((HERE/'reduce_joint.py').read_bytes()) and g['worker']==h((HERE/'joint.py').read_bytes())
            stage_pivots(r,seen,63983)
            if step<19:outgoing+=r['rank']
        assert len(visited)==len(z.namelist())==138 and set(z.namelist())==visited
    joint_rank=len(seen); hidden=63983-incoming-joint_rank; incoming_d1=ranks[1]-ranks[0]-incoming
    assert entries==independent['matrix_entries']==1739328032 and sum(counts.values())==378 and counts['second_Q0']==digest_links==221
    assert ranks==[18003,50147] and incoming==32144 and incoming_d1==0 and outgoing==joint['outgoing_measured_rank']==31230 and joint_rank==joint['joint_measured_rank']==31838 and hidden==1
    return dict(group='SU(3)',source='R-007',level=30,prime=P,original_stage_records=len(visited),matched_packed_matrices=dict(counts),matrix_entries_in_independent_product=entries,matched_right_hand_sides=50,certified_frame_dimensions_by_equation=list(target_sizes.values()),evaluation_coordinates_by_equation=list(widths.values()),source_lists_checked=len(profiles),incoming_Q0_hashes_reconciled=digest_links,first_and_joint_equation_ranks=ranks,incoming_Q0_rank=incoming,incoming_d1_rank_on_recorded_sources_mod_prime=incoming_d1,outgoing_Q0_lower_bound=outgoing,outgoing_plus_Cartan_lower_bound=joint_rank,rational_classical_cohomology_upper_bound=63983-incoming-outgoing,rational_Cartan_invisible_cohomology_upper_bound=hidden,rational_d2_decided=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);a=parser.parse_args();result=json.dumps(run(),indent=2)+'\n';print(result,end='');
    if a.output:a.output.write_text(result)
