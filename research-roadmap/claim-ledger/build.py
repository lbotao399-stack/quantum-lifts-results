#!/usr/bin/env python3
"""Generate canonical claims, the explicit frontier inputs and an audit queue.

Ingestion, semantic review, mathematical validation and gate acceptance are
four different states.  This generator never labels an unreviewed source as
semantically complete merely because its bytes were archived.
"""
import argparse
import gzip
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from canonical import BASELINE, REPOSITORY, claims, frontiers, incoming, native, outgoing, reference

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIMARY_PATHS = ['README.md', 'STATUS.md', 'RESEARCH_DISCUSSION.md',
                 'supplement/supplement.tex', 'nonperturbative-lift/claims.json']
REQUIRED_RELEASES = {'b05','physical','f4'}
RELEASE_ROOTS = {'b05':'algorithms/su3_generator_20261124',
                 'physical':'algorithms/physical_block_20261124',
                 'f4':'algorithms/f4_complete_source_20261124'}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def json_bytes(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + '\n').encode()


def git_show(path, commit=BASELINE):
    return subprocess.check_output(['git', '-C', str(ROOT), 'show', commit + ':' + path])


def write_json(path, value):
    path.write_bytes(json_bytes(value))


def gzip_rows(path, rows):
    with path.open('wb') as stream:
        with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as archive:
            for row in rows:
                archive.write(json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode() + b'\n')


def native_possible(row):
    weight=row['native_weight']; length=row['letter_length']
    if min(weight) < 0 or length < 0:
        return False
    if length == 0:
        return not any(weight)
    return max(weight[:3]) <= length <= sum(weight)


def frontier_manifest(records):
    for record in records:
        base=record['native_input']
        adjacent={}
        maps=[]

        def include(row, role):
            if not native_possible(row):
                return
            key=','.join(map(str,[*row['native_weight'],row['letter_length']]))
            adjacent.setdefault(key,{'native_input':row,'roles':[]})
            if role not in adjacent[key]['roles']:
                adjacent[key]['roles'].append(role)
            for shift,label in [(-1,'classical_predecessor'),(1,'classical_successor')]:
                neighbour=native(row['native_weight'][:3],row['native_weight'][3:],row['letter_length']+shift)
                if native_possible(neighbour):
                    key2=','.join(map(str,[*neighbour['native_weight'],neighbour['letter_length']]))
                    adjacent.setdefault(key2,{'native_input':neighbour,'roles':[]})
                    full_role=role+'/'+label
                    if full_role not in adjacent[key2]['roles']:
                        adjacent[key2]['roles'].append(full_role)

        include(base,'compulsory_target_or_source')
        for direction,maximum in [('incoming',min(base['native_weight'][3:])),
                                   ('outgoing',min(base['native_weight'][:3]))]:
            for order in range(1,maximum+1):
                first=incoming(base,order) if direction=='incoming' else base
                last=base if direction=='incoming' else outgoing(base,order)
                if not native_possible(first) or not native_possible(last):
                    continue
                corrections=[]
                for k in range(order):
                    row=native([v-k for v in first['native_weight'][:3]],
                               [v+k for v in first['native_weight'][3:]],
                               first['letter_length']-k)
                    if not native_possible(row):
                        raise ValueError('Admissible map has an inadmissible correction slot: '+record['id'])
                    corrections.append(row)
                    include(row,direction+'/order%d/correction%d'%(order,k))
                include(last,direction+'/order%d/final_target'%order)
                maps.append({'direction':direction,'order':order,'source':first,'target':last,
                             'all_correction_cochain_slots':corrections,
                             'homogeneous_freedom':'Every classically closed cochain in each listed correction slot remains allowed.',
                             'status':'kinematically_admitted; scientific exclusions require the referenced proof'})
        record['kinematically_admitted_maps']=maps
        record['adjacent_and_correction_slots']=list(adjacent.values())
    return {'schema_version':1,'scope_version':'1.0','source_commit':BASELINE,
            'scope':'Every finite named compulsory native input in SCOPE.md, including those outside the core groups or level cutoff.',
            'not_a_cohomology_count':True,'kinematic_only':True,
            'native_coordinate_order':['theta1','theta2','theta3','n_plus','n_minus'],
            'fields':'All native charges and all doubled physical charges are integers.',
            'homogeneous_correction_rule':'No chosen earlier primitive removes a complete correction cochain slot.',
            'finite_input_count':len(records),
            'finite_adjacent_and_correction_slot_count_with_frontier_multiplicity':sum(len(r['adjacent_and_correction_slots']) for r in records),
            'finite_map_count':sum(len(r['kinematically_admitted_maps']) for r in records),
            'records':records,
            'unbounded_families':[
                {'id':'SU3-F1','claim_id':'su3-f1-formal-tower','parameter':'every integer n>=3',
                 'native_weight':['n+4',4,4,0,0],'letter_length':'n+7','physical_level':'24+2*n',
                 'finite_seeds':['SU3-R007','SU3-F1-second-seed']},
                {'id':'SU3-B1','claim_id':'su3-b1-formal-tower','parameter':'every integer n>=3',
                 'native_weight':['n+5',5,4,0,0],'letter_length':'n+8','physical_level':'28+2*n',
                 'finite_seed_native_inputs':[native([8,5,4],[0,0],11),native([9,5,4],[0,0],12)]},
                {'id':'orthogonal-first-family','claim_id':'orthogonal-first-family',
                 'native_input_reference':reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/first_family.md','d320363f2eb54898e6ff554943b0d6d66705c120'),
                 'parameters':'integers p>=4, N>=7; source closed if N<p+5; target nonzero if p<N-2',
                 'source_native_weight':[1,1,1,'p-1','p-1'],'source_letter_length':'2*p',
                 'target_native_weight':[0,0,0,'p','p'],'target_letter_length':'2*p',
                 'order':1,'coefficient':'3*2^(p-5)*(2*N-3*p-4)/(p-4)!',
                 'nonzero_uniform_subfamilies':[{'N':'p+3','p':'p>=4'},{'N':'p+4','p':'p>=5'}],
                 'exceptional_single_source_zero':{'N':8,'p':4}}]
            }


def located_native_inputs(value, location, group):
    """Index stated integer tuples; do not infer group assignments from prose."""
    if isinstance(value,dict):
        local_group=value.get('group',group)
        weight=value.get('native_weight')
        length=value.get('letter_length')
        if isinstance(weight,list) and len(weight)==5 and all(isinstance(v,int) for v in weight) and isinstance(length,int):
            yield location,value,local_group
        for key,item in value.items():
            yield from located_native_inputs(item,location+'/'+key,local_group)
    elif isinstance(value,list):
        for index,item in enumerate(value):
            yield from located_native_inputs(item,location+'/'+str(index),group)


def advertised_input_manifest(records):
    """Supplement the fixed scope with every explicitly stated finite tuple.

    One tuple may occur in different groups or in a subgroup control. Its
    exact occurrence retains the containing claim's group scope and JSON
    address; no Cartesian group expansion is inferred from a mixed packet.
    Quantified charge domains are indexed verbatim, including uniform families
    and finite aggregate censuses, separately from finite named point inputs.
    """
    inputs={}
    for claim in records:
        for key in ['charges','source','target']:
            for location,row,group in located_native_inputs(claim[key],'/'+key,claim['group']):
                weight=row['native_weight'];length=row['letter_length']
                point=','.join(map(str,[*weight,length]))
                item=inputs.setdefault(point,{'id':'advertised-native-'+point.replace(',','-'),
                    'native_input':native(weight[:3],weight[3:],length),
                    'claim_ids':[],'input_references':[],'uses':[]})
                if claim['id'] not in item['claim_ids']:
                    item['claim_ids'].append(claim['id'])
                item['uses'].append({'claim_id':claim['id'],'json_pointer':location,
                    'stated_group_scope':group,'native_input_with_symmetry_labels':row,
                    'supersession_state':claim['supersession_state'],
                    'scope_rule':'The exact group, subgroup, parity and representation scope is this native occurrence inside its linked claim. A mixed claim group list is not a Cartesian multiplication.'})
                for ref in claim['certificate']:
                    if ref not in item['input_references']:
                        item['input_references'].append(ref)
    points=frontier_manifest([inputs[k] for k in sorted(inputs,key=lambda s:tuple(map(int,s.split(','))))])
    points['scope']='Every explicit finite native tuple in the complete canonical ledger, including historical, withdrawn and additional named frontier inputs outside the fixed core. Each scientific use retains its exact claim and group scope.'
    points['not_a_group_or_result_count']=True
    points['finite_native_occurrence_count']=sum(len(r['uses']) for r in points['records'])
    points['quantified_charge_domains']=[{'claim_id':c['id'],'group':c['group'],
        'charges':c['charges'],'source':c['source'],'target':c['target'],
        'scope':c['scope'],'certificate':c['certificate'],
        'classification':'Exact full statement domain; an aggregate census or uniform family is not a new finite point count.'}
        for c in records]
    points['quantified_domain_count']=len(records)
    return points


def deduplication_manifest(records,audit):
    """Make repeated evidence and non-additive mathematical consequences explicit."""
    by_id={r['id']:r for r in records}
    usage=defaultdict(list)
    seen=set()
    for doc in audit['documents']:
        if doc['sha256'] in seen:
            continue
        seen.add(doc['sha256'])
        review=doc.get('semantic_review') or {}
        for assertion in review.get('covered_assertions',[]):
            for identifier in assertion['claim_ids']:
                usage[identifier].append({'source_text_sha256':doc['sha256'],
                    'location':assertion['location'],'review_packet':review.get('review_packet',doc.get('review_packet'))})
    rows=[]
    for claim in records:
        definition=claim.get('deduplication',{})
        role=definition.get('role')
        if not role:
            role=('historical_or_withdrawn_statement' if claim['supersession_state']!='current'
                  else 'canonical_proposition')
        rows.append({'id':claim['id'],'role':role,'proposition_key':claim['proposition']['key'],
            'canonical_semantic_links':definition.get('canonical_semantic_links',[]),
            'supersedes':claim['supersedes'],'supersession_state':claim['supersession_state'],
            'assertion_uses':usage[claim['id']],
            'scientific_count_increment':None,
            'counting_rule':'This row indexes a proposition or attributed evidence packet. It never increments a target, generator, multiplet, map-rank or physical-lift total.'})
    relations=[
        {'type':'same_target_independent_source_witnesses',
         'members':['so11-power6-d3','so11-power6-d3-codex'],
         'rule':'Two normalized source chains hit the same SO11 sixth-power target; the source comparison is not a second target discovery.'},
        {'type':'product_consequence',
         'members':['so13-power7-d3','so13-power8-d3-product'],
         'rule':'The eighth-power chain is the matched product consequence of the seventh-power chain, not an additional independent generator family.'},
        {'type':'count_then_representative',
         'members':['su3-seven-letter-generator-count','su3-seven-letter-generator-explicit'],
         'rule':'The explicit representative and complete dual realize the already counted one-dimensional quotient.'},
        {'type':'distinct_withdrawn_witnesses_same_target',
         'members':['f4-power6-original-witness-withdrawn','f4-power6-eight-albert-withdrawn'],
         'rule':'The original exact witness and the later eight-Albert sign error are separate historical withdrawals; neither is an active positive second differential.'},
        {'type':'finite_order_result_and_stronger_formal_claim',
         'members':['f4-power6-full-second-zero-report','f4-power6-formal-survival-report'],
         'rule':'The full second-map zero and the separately attributed permanence theorem retain different evidence/acceptance scopes.'},
        {'type':'historical_count_different_scope',
         'members':['historical-n4-campaign-totals-20261118','all-loop-historical-census'],
         'rule':'Campaign pairings, manuscript charge-sector counts and current target results have different counting units and dates; they are not summed.'}]
    for relation in relations:
        if any(identifier not in by_id for identifier in relation['members']):
            raise ValueError('Unknown explicit semantic relation')
    return {'schema_version':1,'record_count':len(rows),
        'scientific_result_count':None,'discovery_count':None,
        'rule':'Text hashes deduplicate identical prose. Canonical identifiers deduplicate assertions. Supporting packets retain every detailed intermediate claim without pretending their row count is a scientific total.',
        'records':rows,'explicit_nonadditive_relations':relations}


def import_nonperturbative():
    raw=git_show('nonperturbative-lift/claims.json')
    original=json.loads(raw)
    records=[]
    excluded=[]
    for item in original['claims']:
        if not item['id'].startswith('n4-'):
            excluded.append({'id':item['id'],'theory':item['theory'],'reason':'Not an advertised undeformed-N4 claim in this source claim index.'})
            continue
        paths=[item['derivation'],*item.get('evidence',[])]
        refs=[]
        for path in paths:
            path,_,anchor=path.partition('#')
            refs.append(reference('nonperturbative-lift/'+path,anchor=anchor or None))
        identifier='np-index-'+item['id']
        state='current'
        proof={'proved':'published_proof','conditional':'conditional_physical',
               'literature-input':'literature_input','open':'open'}[item['status']]
        group={'n4-scalar-detector':'SU(N)',
               'n4-konishi':['SO(7)','Sp(3)=USp(6)','SU(N)'],
               'n4-fourier-positivity':'Every group satisfying the explicitly stated physical hypotheses',
               'n4-nonscalar-image':'All undeformed-N4 groups in this collection',
               'n4-total-descendants':'SU(N)'}[item['id']]
        records.append({'id':identifier,'theory':'undeformed four-dimensional N=4 SYM',
                        'group':group,'title':item['title'],
                        'proposition':{'key':'np-index/'+item['id'],'value':item['statement']},
                        'charges':{'quantified_domain':item['sector'],'hypotheses':item['hypotheses']},
                        'source':{'quantified_domain':item['sector']},
                        'target':{'quantified_domain':item['sector']},
                        'order':item['coefficient_ring'],'direction':'as stated in quantified theorem',
                        'field':item['coefficient_ring'],
                        'normalization':{'native_hypotheses':item['hypotheses'],
                                         'convention_reference':refs[0]},
                        'scope':{'description':item['sector'],'completeness':'Exactly the stated quantified theorem and hypotheses; no extension to all nonscalar classes.',
                                 'physical_identification':'As explicitly stated in the source claim, including its conditions.'},
                        'certificate':refs,'supersession_state':state,'proof_status':proof,
                        'supersedes':[],'limitations':[item['limits']],
                        'unresolved_metadata':[],
                        'advertisements':['nonperturbative-index:'+item['id']],
                        'original_source_record':item})
    return records, {'path':'nonperturbative-lift/claims.json','commit':BASELINE,
                     'sha256':digest(raw),'source_claim_count':len(original['claims']),
                     'in_scope_claim_ids':[r['id'] for r in records],'excluded_records':excluded,
                     'semantic_import_complete_for_this_index':True,
                     'limit':'This index is only one of the public claim sources; this does not certify semantic coverage of its linked documents.'}


def source_units(source):
    """Archive exact nonempty paragraphs with line spans, not inferred claims."""
    text=source['text'];lines=text.splitlines(keepends=True)
    section='preamble';start=None;buffer=[]
    for index,line in enumerate(lines,1):
        heading=re.match(r'^(#{1,6})\s+(.+)',line)
        tex=re.search(r'\\(?:sub)*section\{(.+)\}(?:\\label\{([^}]+)\})?',line)
        if heading:
            section=heading.group(2).strip()
        elif tex:
            section=tex.group(2) or tex.group(1)
        if line.strip():
            if start is None:
                start=index
            buffer.append(line)
        elif buffer:
            body=''.join(buffer)
            yield {'source':source['id'],'line_start':start,'line_end':index-1,
                   'section':section,'text':body,'sha256':digest(body.encode())}
            start=None;buffer=[]
    if buffer:
        body=''.join(buffer)
        yield {'source':source['id'],'line_start':start,'line_end':len(lines),
               'section':section,'text':body,'sha256':digest(body.encode())}


def review_packets():
    packets=[]
    for path in sorted((HERE/'reviews').glob('*.json')):
        data=json.loads(path.read_text())
        if data.get('schema_version')!=1:
            raise ValueError('Unsupported review packet schema: '+str(path))
        packets.append({'path':str(path.relative_to(HERE)),'sha256':digest(path.read_bytes()),'data':data})
    return packets


def release_packets():
    packets=[]
    for path in sorted((HERE/'releases').glob('*.json')):
        data=json.loads(path.read_text())
        if data.get('schema_version')!=1 or not data.get('release_id'):
            raise ValueError('Unsupported release review schema: '+str(path))
        packets.append({'path':str(path.relative_to(HERE)),
            'sha256':digest(path.read_bytes()),'data':data})
    return packets


def apply_release_claims(records,packets):
    for packet in packets:
        records.extend(packet['data'].get('additional_claims',[]))
    by_id={r['id']:r for r in records}
    for packet in packets:
        for update in packet['data'].get('claim_updates',[]):
            claim=by_id[update['id']]
            allowed={'id','append_certificate','proof_status','append_limitations','acceptance_note'}
            if set(update)-allowed:
                raise ValueError('Unsupported release claim update fields: '+str(set(update)-allowed))
            for ref in update.get('append_certificate',[]):
                if ref not in claim['certificate']:
                    claim['certificate'].append(ref)
            claim['limitations'].extend(update.get('append_limitations',[]))
            if 'proof_status' in update:
                claim['proof_status']=update['proof_status']
            if 'acceptance_note' in update:
                claim.setdefault('release_acceptance_notes',[]).append(update['acceptance_note'])


def release_audit(packets):
    documents=[];sources={};ids=[]
    for packet in packets:
        data=packet['data'];release_id=data['release_id'];ids.append(release_id)
        if not re.fullmatch('[0-9a-f]{40}',data.get('commit','')):
            raise ValueError('Release has no immutable publication pin: '+release_id)
        for review in data['documents']:
            commit=review.get('commit',data['commit']);path=review['path']
            if not re.fullmatch('[0-9a-f]{40}',commit) or path.startswith('/') or '..' in Path(path).parts:
                raise ValueError('Invalid release document locator: '+release_id+':'+path)
            if any(s in ('/'+path).lower() for s in ['/isolated/','highly_innovative','/quarantine/','/private/']):
                raise ValueError('Forbidden release source path')
            raw=git_show(path,commit)
            if digest(raw)!=review['sha256']:
                raise ValueError('Release document hash mismatch: '+commit+':'+path)
            source_id=commit+':'+path
            source={'id':source_id,'kind':'git-document','commit':commit,'path':path,
                'sha256':digest(raw),'text':raw.decode(),'release_id':release_id}
            sources[source_id]=source
            documents.append({'release_id':release_id,'source':source_id,
                'commit':commit,'path':path,'sha256':digest(raw),
                'review_packet':packet['path'],'review_packet_sha256':packet['sha256'],
                'semantic_review':review,
                'semantic_coverage':'complete' if review.get('review_complete') and not review.get('metadata_gaps') else 'incomplete'})
    raw=json_bytes(sources)
    (HERE/'release-sources.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    missing=sorted(REQUIRED_RELEASES-set(ids))
    return {'schema_version':1,'required_release_ids':sorted(REQUIRED_RELEASES),
        'release_ids':ids,'missing_release_ids':missing,
        'scope':'Finite post-baseline mathematical release addendum: accepted B05, F02/F03 and D03 evidence. Generated ledger, navigation and score reports do not recursively expand their own review queue.',
        'release_document_count':len(documents),
        'release_distinct_text_count':len({d['sha256'] for d in documents}),
        'semantic_gate':'complete' if not missing and len(set(ids))==len(ids) and all(d['semantic_coverage']=='complete' for d in documents) else 'incomplete',
        'review_packets':[{k:p[k] for k in ['path','sha256']} for p in packets],
        'sources_sha256':digest((HERE/'release-sources.json.gz').read_bytes()),
        'documents':documents}


def source_audit(curated,np_import,packets):
    path=HERE/'inventory'/'public_sources.json.gz'
    inventory_path=HERE/'inventory'/'manifest.json'
    if path.exists() and inventory_path.exists():
        sources=json.loads(gzip.decompress(path.read_bytes()))
        inventory=json.loads(inventory_path.read_text())
        ingestion='completed_at_inventory_pin'
    else:
        sources={}
        for path0 in PRIMARY_PATHS[:-1]:
            raw=git_show(path0);identifier=BASELINE+':'+path0
            sources[identifier]={'id':identifier,'kind':'git-document','commit':BASELINE,
                                 'path':path0,'text':raw.decode(),'sha256':digest(raw)}
        snapshot=json.loads((HERE/'issue-1-snapshot.json').read_text())
        for item in [snapshot['issue'],*snapshot['comments']]:
            identifier='issue-1:'+str(item.get('id','issue-body'))
            sources[identifier]={'id':identifier,'kind':'issue-snapshot','text':item['body'],
                                 'sha256':digest(item['body'].encode())}
        inventory={'document_count':len(sources),'missing_reference_count':None}
        ingestion='incomplete_inventory_not_yet_written'
    reviewed={}
    for packet in packets:
        for review in packet['data'].get('documents',[]):
            if review['sha256'] in reviewed:
                raise ValueError('Two reviews claim the same document text: '+review['sha256'])
            reviewed[review['sha256']]={**review,'review_packet':packet['path'],
                                       'review_packet_sha256':packet['sha256']}
    units=[];documents=[];groups={};unit_count_by_hash={}
    for identifier,source in sorted(sources.items()):
        groups.setdefault(source['sha256'],[]).append(identifier)
    # Identical document bytes can share a textual review.  Their immutable
    # occurrences remain distinct: relative evidence links may point at
    # different underlying code at different commits.
    canonical_sources={key:next((i for i in identifiers if i.startswith(BASELINE+':')),identifiers[0])
                       for key,identifiers in groups.items()}
    for identifier,source in sorted(sources.items()):
        source.setdefault('id',identifier)
        linked=[]
        for claim in curated:
            if any(r.get('kind')=='git' and r.get('commit')==source.get('commit') and r.get('path')==source.get('path') for r in claim['certificate']):
                linked.append(claim['id'])
        count=0
        review=reviewed.get(source['sha256'])
        complete=bool(review and review.get('review_complete') and not review.get('metadata_gaps'))
        for unit in source_units(source):
            count+=1
            if canonical_sources[source['sha256']]!=identifier:
                continue
            unit['id']='span-'+digest((identifier+':'+str(unit['line_start'])+':'+unit['sha256']).encode())[:24]
            unit['unit']='Exact source paragraph; it may contain zero, one or several claims.'
            unit['semantic_disposition']=('reviewed_document_context' if review['classification'] in ('in_scope','mixed') else 'excluded_scope') if complete else 'unreviewed'
            unit['claim_ids']=review.get('claim_ids',[]) if complete else []
            if complete:
                unit['document_review_packet']=review['review_packet']
            unit['identical_text_source_occurrences']=groups[source['sha256']]
            unit['duplicate_rule']='Textual equality only; immutable evidence and relative-link contexts remain separate.'
            units.append(unit)
        unit_count_by_hash[source['sha256']]=count
        documents.append({'source':identifier,'kind':source['kind'],'sha256':source['sha256'],
                          'source_span_count':count,'claims_with_this_certificate':linked,
                          'textual_review_source':canonical_sources[source['sha256']],
                          'semantic_coverage':'complete' if complete else 'not_certified',
                          'semantic_review':review,
                          'remaining_action':'Review every in-scope assertion, connect it to canonical claims, and explicitly classify every other span.'})
    gzip_rows(HERE/'advertisement-spans.jsonl.gz',units)
    return {'schema_version':1,'baseline':BASELINE,'ingestion_state':ingestion,
            'ingested_source_count':len(documents),'distinct_source_content_count':len(groups),
            'source_span_count':len(units),
            'source_span_occurrence_count':sum(d['source_span_count'] for d in documents),
            'definition':'Source spans are an audit queue, not a claim count and not evidence of complete semantic coverage.',
            'semantic_gate':'complete' if all(d['semantic_coverage']=='complete' for d in documents) else 'unresolved',
            'semantically_complete_document_count':sum(d['semantic_coverage']=='complete' for d in documents),
            'semantically_complete_distinct_contents':len({d['sha256'] for d in documents if d['semantic_coverage']=='complete'}),
            'unreviewed_span_count':sum(u['semantic_disposition']=='unreviewed' for u in units),'documents':documents,
            'review_packets':[{k:p[k] for k in ('path','sha256')} for p in packets],
            'primary_sources':PRIMARY_PATHS,'nonperturbative_index':np_import,
            'inventory_missing_reference_count':inventory.get('missing_reference_count'),
            'advertisement_spans_sha256':digest((HERE/'advertisement-spans.jsonl.gz').read_bytes())}


def summary(ledger,audit,frontier,core,advertised=None,dedup=None,releases=None):
    counts=Counter(r['proof_status'] for r in ledger['claims'])
    unresolved=[r for r in ledger['claims'] if r['unresolved_metadata']]
    complete=ledger['gate_status']=='complete'
    lines=['# Canonical public claim ledger','',
           ('**A07 semantic/status coverage is complete.** Historical attribution, mathematical acceptance, selected-family scope and physical hypotheses remain explicit.' if complete else '**A07 remains incomplete.** This artifact separates normalized claim records, archived source spans, mathematical evidence and final gate acceptance.'),'',
           'The immutable public baseline is `'+BASELINE+'`; Issue 1 includes the newer sixth-power proof at `c778cd7c3ceb72ad8082d791b2aee5f0df082164`.','',
           '## Exact state','',
           '| Item | Count or state |','|---|---:|',
           '| Canonical semantic records | %d |'%len(ledger['claims']),
           '| Records with explicitly unresolved metadata | %d |'%len(unresolved),
           '| Baseline public document / issue occurrences | %d |'%audit['ingested_source_count'],
           '| Baseline distinct document contents | %d |'%audit['distinct_source_content_count'],
           '| Post-baseline mathematical release documents | %d |'%(releases['release_document_count'] if releases else 0),
           '| Post-baseline release semantic gate | %s |'%(releases['semantic_gate'] if releases else 'not_generated'),
           '| Distinct source paragraphs awaiting semantic disposition | %d |'%audit['unreviewed_span_count'],
           '| Documents certified semantically exhaustive | %d |'%audit['semantically_complete_document_count'],
           '| Distinct contents certified semantically exhaustive | %d |'%audit['semantically_complete_distinct_contents'],
           '| Named finite compulsory native inputs | %d |'%frontier['finite_input_count'],
           '| Named-frontier native adjacent / correction slots, counting each frontier separately | %d |'%frontier['finite_adjacent_and_correction_slot_count_with_frontier_multiplicity'],
           '| Named-frontier kinematically admitted incoming / outgoing maps | %d |'%frontier['finite_map_count'],
           '| Additional index: distinct explicitly stated native tuples | %d |'%(advertised['finite_input_count'] if advertised else 0),
           '| Additional index: native tuple uses, with exact claim/group scope | %d |'%(advertised['finite_native_occurrence_count'] if advertised else 0),'',
           'Document ingestion is not semantic claim coverage. The paragraph archive contains the exact public statements and their document-level semantic review decisions; it does not assert that each paragraph is a claim. No discovery or lift percentage is reported.','',
           '## Current scientific distinctions','',
           '- The SO11 sixth-power, SO12 seventh-power and SO13 seventh-power third-order witnesses are separate targets. The SO13 eighth-power hit is a product consequence. The two normalized SO11 source constructions do not create a second target.',
           '- The eight-Albert F4 witness is withdrawn. Its two corrected modular zeros are recorded with their limited scope. The root accepted the full incoming sixth-power second-map zero after auditing the published complete-source theorem and an independent closed-factor proof. The published full 3017 minor is an explicit prerequisite, not a newly reconstructed minor. The separate all-page permanence assertion retains its own proof status.',
           '- F4 seventh-power first-page survival and eighth-power first-order nonvanishing retain their native Ward ranks and minors. A selected 518-recipe seventh-power zero does not decide the complete incoming source.',
           '- R-007, the F1 and B1 families, and the nine-letter SU3 target have formal results. Their physical instanton matrix elements remain separate. The explicit seven-letter SU3 representative, complete dual and finite-rank fortuity certificate are accepted and pinned at `eac7db42a6e27f4900826fe7b42cc75365698dc3`.',
           '- Historical modular censuses retain their fields and dated scope. Unavailable row-specific rational receipts remain stated evidence limitations; no fresh rational theorem is inferred. No new SU2, Sp1 or SO3 calculation is scheduled.',
           '- Physical partial-profile zeros, conditional protected-product arguments and long-multiplet controls are not positive microscopic BPS instanton lifts.','',
           '## Finite coverage manifest','']
    if core:
        lines += ['The core has **%d native charge/length slots per group** and **%d complete physical charge blocks per group** over **%d groups**. The explicit Cartesian expansion therefore contains **%d native group slots** and **%d physical group blocks**.'%(core['native_slots_per_group'],core['physical_charge_blocks_per_group'],len(core['groups']),core['expanded_native_group_slots'],core['expanded_physical_group_blocks']),
                  '', 'These are kinematic slots. Gauge-invariant dimensions, classical cohomology dimensions and resolved discovery/lift denominators have not been computed by this enumeration. The maximum kinematically allowed incoming and outgoing orders in the core are both %d. The independent alphabet enumeration checks %d slots through level %d.'%(core['maximum_kinematically_allowed_incoming_order'],core['independent_alphabet_check_slots'],core['independent_alphabet_check_level']), '']
    lines += ['`native-frontiers.json` supplies the required native charges outside the core and every admissible source/correction slot, with complete homogeneous freedom retained. Family parameters are stored separately from finite inputs.',
              '`advertised-native-inputs.json` additionally indexes every explicit finite native tuple in the canonical claims, its exact source address and group/representation context, all admissible neighbours and correction slots, and every full quantified charge domain. Its tuple count is not multiplied by a mixed packet\'s group list.',
              '`deduplication.json` records repeated assertion uses, supporting evidence packets, historical versions and explicit non-additive relations. The record count is not a count of scientific discoveries.','',
              '## Proof-state counts','', '| Proof state | Records |','|---|---:|']
    lines += ['| '+state+' | '+str(count)+' |' for state,count in sorted(counts.items())]
    lines += ['','## Active semantic records','',
              '| Identifier | Group | Proof state | Exact recorded outcome |',
              '|---|---|---|---|']
    for record in ledger['claims']:
        if record['supersession_state']!='current':
            continue
        group=record['group'] if isinstance(record['group'],str) else ', '.join(record['group'])
        value=json.dumps(record['proposition']['value'],ensure_ascii=False,sort_keys=True)
        lines.append('| `'+record['id']+'` | '+group.replace('|','\\|')+' | '+record['proof_status']+' | '+value.replace('|','\\|').replace('\n',' ')+' |')
    lines += ['','## Audit boundary','',
              'The baseline finite source queue contains the pinned public graph and frozen Issue 1 snapshot. All textual variants keep their immutable occurrence and relative-reference context. `release-audit.json` and `release-sources.json.gz` authenticate the separate finite B05, F02/F03 and D03 publication addendum. Later publications require a dated update; they do not silently replace historical statements. Generated ledger/navigation reports are excluded from recursive self-review.',
              'A complete semantic/status audit does not re-prove every historical result or complete the classical, perturbative, physical or instanton censuses. Open, conditional, modular-only and unavailable-replay states remain visible.','',
              '## Explicit unresolved metadata','']
    if not unresolved:
        lines.append('None in the fixed semantic queue. Remaining mathematical and physical limitations are attached to the individual claims.')
    for record in unresolved:
        lines.append('- `'+record['id']+'`: '+' '.join(record['unresolved_metadata']))
    lines += ['','## Reproduce and check','',
              '```sh','python3 research-roadmap/claim-ledger/build.py',
              'python3 research-roadmap/claim-ledger/check.py --negative-controls',
              'python3 research-roadmap/claim-ledger/check.py --require-complete','```','',
              'The checker validates the ledger structure, references, native charges, complete source coverage, explicit coverage manifests, semantic deduplication and deliberate corruption controls. The strict command passes only when every fixed A07 completeness condition is recorded without metadata gaps. No roadmap score is edited by these scripts.','']
    return '\n'.join(lines)


def build():
    records=claims();np_records,np_import=import_nonperturbative();records.extend(np_records)
    packets=review_packets()
    for packet in packets:
        records.extend(packet['data'].get('additional_claims',[]))
    release_inputs=release_packets()
    apply_release_claims(records,release_inputs)
    releases=release_audit(release_inputs)
    frontier=frontier_manifest(frontiers())
    ledger={'schema_version':1,'repository':REPOSITORY,'baseline':BASELINE,
            'scope':'Every advertised undeformed-N4 claim in the fixed reachable public evidence graph, including all historical text versions and discussion-linked immutable commits outside main. Scientific acceptance remains distinct from semantic coverage.',
            'gate_status':'incomplete','counting_rule':'Targets, source witnesses, map ranks, multiplets, generators, products, formal existence and physical matrix elements remain separate.',
            'claims':records}
    audit=source_audit(records,np_import,packets)
    ledger['gate_status']='complete' if (audit['semantic_gate']=='complete'
        and audit['ingestion_state']=='completed_at_inventory_pin'
        and audit['inventory_missing_reference_count']==0
        and releases['semantic_gate']=='complete'
        and all(not r['unresolved_metadata'] for r in records)) else 'incomplete'
    core_path=HERE/'coverage'/'coverage_manifest.json'
    core=json.loads(core_path.read_text()) if core_path.exists() else None
    advertised=advertised_input_manifest(records)
    dedup=deduplication_manifest(records,{'documents':audit['documents']+releases['documents']})
    write_json(HERE/'claims.json',ledger)
    write_json(HERE/'native-frontiers.json',frontier)
    write_json(HERE/'source-audit.json',audit)
    write_json(HERE/'release-audit.json',releases)
    write_json(HERE/'advertised-native-inputs.json',advertised)
    write_json(HERE/'deduplication.json',dedup)
    coverage={'schema_version':1,'scope_version':'1.0','status':'kinematic_manifest_complete; semantic_A07_gate_'+ledger['gate_status'],
              'core_manifest':'coverage/coverage_manifest.json','native_frontiers':'native-frontiers.json',
              'core_manifest_sha256':digest(core_path.read_bytes()) if core else None,
              'native_frontiers_sha256':digest((HERE/'native-frontiers.json').read_bytes()),
              'advertised_native_inputs':'advertised-native-inputs.json',
              'advertised_native_inputs_sha256':digest((HERE/'advertised-native-inputs.json').read_bytes()),
              'deduplication':'deduplication.json',
              'deduplication_sha256':digest((HERE/'deduplication.json').read_bytes()),
              'release_audit':'release-audit.json',
              'release_audit_sha256':digest((HERE/'release-audit.json').read_bytes()),
              'source_of_scope':'research-roadmap/scope.json',
              'exact_classical_cohomology_dimensions_complete':False,
              'discovery_coverage_percentage':None,'lift_coverage_percentage':None}
    write_json(HERE/'coverage-blocks.json',coverage)
    (HERE/'CURRENT_SUMMARY.md').write_text(summary(ledger,audit,frontier,core,advertised,dedup,releases))
    print(json.dumps({'canonical_claims':len(records),'frontier_inputs':len(frontier['records']),
                      'source_spans':audit['source_span_count'],'ingestion':audit['ingestion_state'],
                      'gate_status':ledger['gate_status']},sort_keys=True))


if __name__=='__main__':
    build()
