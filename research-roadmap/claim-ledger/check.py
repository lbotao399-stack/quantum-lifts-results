#!/usr/bin/env python3
"""Check semantic ledger consistency and reject synthetic corruptions.

Exit 0 means the stored partial ledger is internally valid.  Only
--require-complete additionally demands the fixed A07 semantic coverage gate.
No source-document count is used as a scientific coverage denominator.
"""
import argparse
import copy
import gzip
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from build import (HERE, ROOT, digest, native_possible, summary,
                   located_native_inputs, deduplication_manifest, REQUIRED_RELEASES, RELEASE_ROOTS)
from canonical import native

REQUIRED=['id','theory','group','title','proposition','charges','source','target',
          'order','direction','field','normalization','scope','certificate',
          'supersession_state','proof_status','supersedes','limitations',
          'unresolved_metadata','advertisements']
STATES={'current','historical','withdrawn','superseded'}
BLOCKED=('/isolated/','highly_innovative','/quarantine/','/private/')


def refs(value):
    if isinstance(value,dict):
        if value.get('kind') in ('git','workspace','external'):
            yield value
        for item in value.values():
            yield from refs(item)
    elif isinstance(value,list):
        for item in value:
            yield from refs(item)


def natives(value):
    if isinstance(value,dict):
        if (isinstance(value.get('native_weight'),list)
                and len(value['native_weight'])==5
                and all(isinstance(v,int) for v in value['native_weight'])
                and isinstance(value.get('letter_length'),int)):
            yield value
        for item in value.values():
            yield from natives(item)
    elif isinstance(value,list):
        for item in value:
            yield from natives(item)


@lru_cache(maxsize=None)
def git_object(commit,path):
    result=subprocess.run(['git','-C',str(ROOT),'cat-file','-t',commit+':'+path],
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    return result.returncode,result.stdout.decode().strip()


@lru_cache(maxsize=None)
def git_text(commit,path):
    result=subprocess.run(['git','-C',str(ROOT),'show',commit+':'+path],
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    return result.stdout.decode(errors='replace') if not result.returncode else ''


def anchor_exists(text,anchor):
    if re.fullmatch(r'L\d+(?:-L\d+)?',anchor):
        numbers=[int(x) for x in re.findall(r'\d+',anchor)]
        return min(numbers)>=1 and max(numbers)<=len(text.splitlines())
    if re.search(r'(?:id|name)=["\']'+re.escape(anchor)+r'["\']',text):
        return True
    if '\\label{'+anchor+'}' in text:
        return True
    slugs=[]
    for line in text.splitlines():
        match=re.match(r'^#{1,6}\s+(.+)',line)
        if match:
            label=re.sub(r'<[^>]+>','',match.group(1)).lower().strip()
            label=re.sub(r'[^\w\- ]','',label)
            slugs.append(label.replace(' ','-'))
    return anchor in slugs


def validate_reference(ref):
    errors=[]
    path=ref.get('path','')
    if any(x in ('/'+path).lower() for x in BLOCKED):
        return ['quarantined reference is forbidden: '+path]
    if ref['kind']=='git':
        commit=ref.get('commit','')
        if not re.fullmatch('[0-9a-f]{40}',commit):
            errors.append('non-immutable native reference: '+str(ref))
        elif not path or path.startswith('/') or '..' in Path(path).parts:
            errors.append('invalid native reference path: '+path)
        else:
            code,kind=git_object(commit,path)
            if code or kind not in ('blob','tree'):
                errors.append('missing native reference: '+commit+':'+path)
            elif ref.get('anchor') and kind=='blob' and path.endswith(('.md','.tex')):
                if not anchor_exists(git_text(commit,path),ref['anchor']):
                    errors.append('missing native reference anchor: '+commit+':'+path+'#'+ref['anchor'])
    elif ref['kind']=='workspace':
        if not path or path.startswith('/') or '..' in Path(path).parts or not (ROOT/path).is_file():
            errors.append('missing workspace reference: '+path)
    elif ref['kind']=='external' and not ref.get('url'):
        errors.append('external reference has no URL')
    return errors


def validate_claims(ledger,check_references=True):
    errors=[];ids={};active=defaultdict(list)
    if ledger.get('schema_version')!=1:
        errors.append('unsupported claim schema')
    for claim in ledger.get('claims',[]):
        identifier=claim.get('id','<missing>')
        if identifier in ids:
            errors.append('duplicate claim id: '+identifier)
        ids[identifier]=claim
        for field in REQUIRED:
            if field not in claim:
                errors.append(identifier+': missing required field '+field)
        if any(field not in claim for field in REQUIRED):
            continue
        for field in ['theory','group','title','proposition','charges','source','target','order','field','normalization','scope','certificate']:
            if claim[field] in (None,'',[],{}):
                errors.append(identifier+': empty required field '+field)
        if claim['supersession_state'] not in STATES:
            errors.append(identifier+': invalid supersession state')
        if claim['supersession_state']=='current':
            prop=claim['proposition']
            if 'key' not in prop or 'value' not in prop:
                errors.append(identifier+': proposition needs key and value')
            else:
                active[(prop['key'],json.dumps(claim['group'],sort_keys=True),
                        json.dumps(claim['field'],sort_keys=True))].append(claim)
        if check_references:
            seen=set()
            for ref in refs(claim):
                key=json.dumps(ref,sort_keys=True)
                if key not in seen:
                    seen.add(key)
                    errors.extend(identifier+': '+err for err in validate_reference(ref))
        for row in natives(claim):
            if not native_possible(row):
                errors.append(identifier+': invalid native charge/length slot '+str(row))
                continue
            rebuilt=native(row['native_weight'][:3],row['native_weight'][3:],row['letter_length'])
            for key in ['twice_energy','twice_R','twice_J_left','twice_J_right','physical_level']:
                if key in row and row[key]!=rebuilt[key]:
                    errors.append(identifier+': inconsistent native-to-physical charge '+key)
    for key,records in active.items():
        values={json.dumps(r['proposition']['value'],sort_keys=True) for r in records}
        if len(values)>1:
            errors.append('conflicting active claims for '+key[0]+': '+', '.join(r['id'] for r in records))
    for identifier,claim in ids.items():
        for old in claim.get('supersedes',[]):
            if old not in ids:
                errors.append(identifier+': unknown superseded claim '+old)
            elif ids[old].get('supersession_state')=='current':
                errors.append(identifier+': superseded record is still current: '+old)
    visiting=set();done=set()
    def visit(identifier):
        if identifier in visiting:
            errors.append('cyclic supersession at '+identifier);return
        if identifier in done or identifier not in ids:
            return
        visiting.add(identifier)
        for old in ids[identifier].get('supersedes',[]):
            visit(old)
        visiting.remove(identifier);done.add(identifier)
    for identifier in ids:
        visit(identifier)
    return errors


def validate_frontiers(frontier,ledger):
    errors=[];ids={c['id'] for c in ledger['claims']}
    if frontier['finite_input_count']!=len(frontier['records']):
        errors.append('frontier input count mismatch')
    for record in frontier['records']:
        for claim_id in record['claim_ids']:
            if claim_id not in ids:
                errors.append('unknown frontier claim id '+claim_id)
        for row in natives(record):
            if not native_possible(row):
                errors.append('inadmissible frontier slot '+record['id'])
        for item in record['kinematically_admitted_maps']:
            order=item['order'];source=item['source'];target=item['target']
            if target['native_weight']!=[v-order for v in source['native_weight'][:3]]+[v+order for v in source['native_weight'][3:]]:
                errors.append('wrong frontier map charges '+record['id'])
            if target['letter_length']!=source['letter_length']+1-order:
                errors.append('wrong frontier map length '+record['id'])
            corrections=item['all_correction_cochain_slots']
            if len(corrections)!=order:
                errors.append('missing correction freedom '+record['id'])
            for k,row in enumerate(corrections):
                if row['native_weight']!=[v-k for v in source['native_weight'][:3]]+[v+k for v in source['native_weight'][3:]] or row['letter_length']!=source['letter_length']-k:
                    errors.append('wrong correction slot '+record['id'])
    return errors


def validate_coverage(core,coverage):
    errors=[]
    if core is None:
        return ['missing core coverage manifest']
    for filename,expected in core['sha256'].items():
        path=HERE/'coverage'/filename
        if not path.is_file() or digest(path.read_bytes())!=expected:
            errors.append('coverage payload digest mismatch: '+filename)
    if core['expanded_native_group_slots']!=len(core['groups'])*core['native_slots_per_group']:
        errors.append('native Cartesian expansion count mismatch')
    if core['expanded_physical_group_blocks']!=len(core['groups'])*core['physical_charge_blocks_per_group']:
        errors.append('physical Cartesian expansion count mismatch')
    if sum(core['native_slots_by_level'].values())!=core['native_slots_per_group']:
        errors.append('level counts do not sum to native slots')
    for filename,key in [('native_slots.jsonl.gz','native_slots_per_group'),('physical_blocks.jsonl.gz','physical_charge_blocks_per_group')]:
        with gzip.open(HERE/'coverage'/filename,'rt') as stream:
            count=sum(1 for _ in stream)
        if count!=core[key]:
            errors.append('explicit coverage row count mismatch: '+filename)
    for filename,key in [('coverage/coverage_manifest.json','core_manifest_sha256'),('native-frontiers.json','native_frontiers_sha256'),
                         ('advertised-native-inputs.json','advertised_native_inputs_sha256'),('deduplication.json','deduplication_sha256'),
                         ('release-audit.json','release_audit_sha256')]:
        if digest((HERE/filename).read_bytes())!=coverage[key]:
            errors.append('coverage index digest mismatch: '+filename)
    if coverage.get('discovery_coverage_percentage') is not None or coverage.get('lift_coverage_percentage') is not None:
        errors.append('kinematic counts cannot imply a discovery or lift coverage percentage')
    return errors


def validate_advertised_inputs(advertised,ledger):
    errors=validate_frontiers(advertised,ledger)
    expected={}
    for claim in ledger['claims']:
        for key in ['charges','source','target']:
            for location,row,group in located_native_inputs(claim[key],'/'+key,claim['group']):
                expected[(claim['id'],location)]=(row,group)
    actual={}
    for record in advertised['records']:
        for use in record['uses']:
            key=(use['claim_id'],use['json_pointer'])
            if key in actual:
                errors.append('duplicate advertised native occurrence '+str(key))
            actual[key]=(use['native_input_with_symmetry_labels'],use['stated_group_scope'])
            if use['native_input_with_symmetry_labels']['native_weight']!=record['native_input']['native_weight'] or use['native_input_with_symmetry_labels']['letter_length']!=record['native_input']['letter_length']:
                errors.append('advertised occurrence assigned to wrong native tuple '+str(key))
    if actual!=expected:
        errors.append('advertised native occurrence coverage mismatch')
    if advertised['finite_native_occurrence_count']!=len(actual):
        errors.append('advertised native occurrence count mismatch')
    domain={d['claim_id']:d for d in advertised['quantified_charge_domains']}
    if set(domain)!={c['id'] for c in ledger['claims']}:
        errors.append('advertised quantified charge domain coverage mismatch')
    for claim in ledger['claims']:
        if claim['id'] in domain and any(domain[claim['id']][key]!=claim[key] for key in ['group','charges','source','target','scope','certificate']):
            errors.append('advertised quantified domain differs from canonical claim '+claim['id'])
    return errors


def validate_deduplication(dedup,ledger,audit,releases):
    errors=[]
    if dedup!=deduplication_manifest(ledger['claims'],{'documents':audit['documents']+releases['documents']}):
        errors.append('semantic deduplication manifest is stale or incomplete')
    identifiers={c['id'] for c in ledger['claims']}
    for record in dedup['records']:
        if record['scientific_count_increment'] is not None:
            errors.append('semantic record count is being used as a scientific increment')
        for target in record['canonical_semantic_links']:
            if target not in identifiers:
                errors.append('unknown canonical deduplication target '+target)
    return errors


def validate_releases(releases,ledger):
    errors=[];packets={};ids={c['id'] for c in ledger['claims']}
    for packet in releases['review_packets']:
        path=HERE/packet['path']
        if not path.is_file() or digest(path.read_bytes())!=packet['sha256']:
            errors.append('release review packet digest mismatch '+packet['path']);continue
        data=json.loads(path.read_text())
        release_id=data['release_id']
        if release_id in packets:
            errors.append('duplicate release identifier '+release_id)
        packets[release_id]=data
        if release_id not in RELEASE_ROOTS:
            errors.append('unscoped post-baseline release '+release_id);continue
        result=subprocess.run(['git','-C',str(ROOT),'ls-tree','-r','--name-only',data['commit'],'--',RELEASE_ROOTS[release_id]],
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        expected={name for name in result.stdout.decode().splitlines() if name.endswith(('.md','.tex'))}
        actual={d['path'] for d in data['documents']}
        if result.returncode or not expected or actual!=expected:
            errors.append('release mathematical text inventory incomplete '+release_id)
    archive=HERE/'release-sources.json.gz'
    if not archive.is_file() or digest(archive.read_bytes())!=releases['sources_sha256']:
        return errors+['release source archive digest mismatch']
    sources=json.loads(gzip.decompress(archive.read_bytes()))
    if set(packets)!=set(releases['release_ids']):
        errors.append('declared release identifiers differ from actual review packets')
    expected_occurrences={(data['release_id'],doc.get('commit',data['commit'])+':'+doc['path'])
        for data in packets.values() for doc in data['documents']}
    actual_occurrences={(doc['release_id'],doc['source']) for doc in releases['documents']}
    if expected_occurrences!=actual_occurrences or len(actual_occurrences)!=len(releases['documents']):
        errors.append('release audit occurrence list differs from complete review packets')
    if set(sources)!={d['source'] for d in releases['documents']}:
        errors.append('release source occurrence coverage mismatch')
    if releases['release_document_count']!=len(releases['documents']):
        errors.append('release document count mismatch')
    if releases['release_distinct_text_count']!=len({d['sha256'] for d in releases['documents']}):
        errors.append('release distinct-text count mismatch')
    for doc in releases['documents']:
        review=doc['semantic_review'];source=sources.get(doc['source'])
        packet=packets.get(doc['release_id'])
        if not packet or review not in packet['documents']:
            errors.append('release review differs from pinned packet '+doc['source'])
        if not source or digest(source['text'].encode())!=doc['sha256'] or digest(git_text(doc['commit'],doc['path']).encode())!=doc['sha256']:
            errors.append('release document bytes differ from immutable source '+doc['source'])
        if review['classification'] not in ['in_scope','mixed','out_of_scope','navigation']:
            errors.append('unknown release semantic classification '+doc['source'])
        if review.get('review_complete') and review['classification'] in ['in_scope','mixed'] and (not review.get('claim_ids') or not review.get('covered_assertions')):
            errors.append('release mathematical text lacks semantic assertions '+doc['source'])
        if review['classification'] in ['out_of_scope','navigation'] and not review.get('scope_exclusions'):
            errors.append('release excluded text lacks scope reason '+doc['source'])
        for identifier in review.get('claim_ids',[]):
            if identifier not in ids:
                errors.append('release refers to absent canonical claim '+identifier)
        for assertion in review.get('covered_assertions',[]):
            if not assertion.get('location') or not assertion.get('description'):
                errors.append('release assertion lacks exact location '+doc['source'])
            for identifier in assertion.get('claim_ids',[]):
                if identifier not in review.get('claim_ids',[]):
                    errors.append('release assertion absent from document claim list '+identifier)
    return errors


def validate_source_audit(audit,ledger):
    errors=[]
    if digest((HERE/'advertisement-spans.jsonl.gz').read_bytes())!=audit['advertisement_spans_sha256']:
        errors.append('advertisement archive digest mismatch')
    with gzip.open(HERE/'advertisement-spans.jsonl.gz','rt') as stream:
        rows=[json.loads(line) for line in stream]
    if len(rows)!=audit['source_span_count']:
        errors.append('source span count mismatch')
    if sum(row['semantic_disposition']=='unreviewed' for row in rows)!=audit['unreviewed_span_count']:
        errors.append('unreviewed span count mismatch')
    ids=set()
    for row in rows:
        if row['id'] in ids:
            errors.append('duplicate source span id '+row['id'])
        ids.add(row['id'])
        if digest(row['text'].encode())!=row['sha256']:
            errors.append('source span text digest mismatch '+row['id'])
    claim_ids={r['id'] for r in ledger['claims']}
    complete_reviews={}
    for packet in audit.get('review_packets',[]):
        path=HERE/packet['path']
        if not path.is_file() or digest(path.read_bytes())!=packet['sha256']:
            errors.append('semantic review packet digest mismatch '+packet['path'])
        else:
            for review in json.loads(path.read_text())['documents']:
                if review['sha256'] in complete_reviews:
                    errors.append('duplicate semantic review text '+review['sha256'])
                complete_reviews[review['sha256']]=review
    actual_complete=0;complete_hashes=set()
    for document in audit['documents']:
        review=document.get('semantic_review')
        if review and review.get('review_complete') and not review.get('metadata_gaps'):
            actual_complete+=1;complete_hashes.add(document['sha256'])
            if document['semantic_coverage']!='complete':
                errors.append('source review state mismatch '+document['source'])
        elif document['semantic_coverage']=='complete':
            errors.append('unreviewed document incorrectly marked complete '+document['source'])
        if review:
            stripped={k:v for k,v in review.items() if k not in ['review_packet','review_packet_sha256']}
            if stripped!=complete_reviews.get(document['sha256']):
                errors.append('embedded source review differs from signed review packet '+document['source'])
    if actual_complete!=audit['semantically_complete_document_count'] or len(complete_hashes)!=audit['semantically_complete_distinct_contents']:
        errors.append('semantic coverage counts do not match actual document reviews')
    checked=set()
    for document in audit['documents']:
        if document['sha256'] in checked:
            continue
        checked.add(document['sha256'])
        review=document.get('semantic_review')
        if not review:
            continue
        if review['classification'] not in ['in_scope','mixed','out_of_scope','navigation']:
            errors.append('unknown semantic document classification '+document['source'])
        if review.get('review_complete'):
            if review['classification'] in ['in_scope','mixed'] and (not review.get('claim_ids') or not review.get('covered_assertions')):
                errors.append('complete in-scope review lacks semantic assertions '+document['source'])
            if review['classification'] in ['out_of_scope','navigation'] and not review.get('scope_exclusions'):
                errors.append('complete excluded review lacks an explicit reason '+document['source'])
        for identifier in review.get('claim_ids',[]):
            if identifier not in claim_ids:
                errors.append('review refers to absent canonical claim '+identifier)
        for assertion in review.get('covered_assertions',[]):
            if not assertion.get('location') or not assertion.get('description'):
                errors.append('review assertion lacks a literal source location '+document['source'])
            for identifier in assertion.get('claim_ids',[]):
                if identifier not in review.get('claim_ids',[]):
                    errors.append('assertion mapping omitted from document claim list '+identifier)
    inventory_path=HERE/'inventory'/'manifest.json'
    if inventory_path.exists():
        inventory=json.loads(inventory_path.read_text())
        archive=HERE/'inventory'/'public_sources.json.gz'
        if not archive.is_file() or digest(archive.read_bytes())!=inventory['sources_sha256']:
            errors.append('public source inventory digest mismatch')
        elif audit['ingestion_state']=='completed_at_inventory_pin':
            sources=json.loads(gzip.decompress(archive.read_bytes()))
            for key,source in sources.items():
                if digest(source['text'].encode())!=source['sha256']:
                    errors.append('source inventory text digest mismatch '+key)
            if len(sources)!=audit['ingested_source_count']:
                errors.append('source audit is stale relative to completed inventory')
    return errors


def gate_blockers(ledger,audit,frontier,releases=None):
    blockers=[]
    unresolved=[r['id'] for r in ledger['claims'] if r['unresolved_metadata']]
    if unresolved:
        blockers.append({'kind':'unresolved_claim_metadata','count':len(unresolved),'claim_ids':unresolved})
    if audit['unreviewed_span_count']:
        blockers.append({'kind':'semantic_source_review_incomplete','unreviewed_source_spans':audit['unreviewed_span_count'],
                         'complete_documents':audit['semantically_complete_document_count']})
    incomplete={d['sha256'] for d in audit['documents'] if not d.get('semantic_review') or
        not d['semantic_review'].get('review_complete') or d['semantic_review'].get('metadata_gaps')}
    if incomplete:
        blockers.append({'kind':'semantic_document_review_incomplete','distinct_text_count':len(incomplete)})
    if audit['ingestion_state']!='completed_at_inventory_pin':
        blockers.append({'kind':'public_graph_ingestion_incomplete','state':audit['ingestion_state']})
    if audit['inventory_missing_reference_count'] is None or audit['inventory_missing_reference_count']:
        blockers.append({'kind':'missing_public_graph_references','count':audit['inventory_missing_reference_count']})
    for item in frontier['unbounded_families']:
        if item.get('unresolved'):
            blockers.append({'kind':'unresolved_uniform_family_input','id':item['id'],'reason':item['unresolved']})
    if releases is not None:
        missing=sorted(REQUIRED_RELEASES-set(releases['release_ids']))
        if missing:
            blockers.append({'kind':'post_baseline_release_review_missing','release_ids':missing})
        incomplete=[d['source'] for d in releases['documents'] if not d['semantic_review'].get('review_complete') or d['semantic_review'].get('metadata_gaps')]
        if incomplete:
            blockers.append({'kind':'post_baseline_release_review_incomplete','sources':incomplete})
    return blockers


def negative_controls(ledger,frontier,audit,advertised,releases):
    results=[]
    baseline=copy.deepcopy(ledger)
    chosen=next(c for c in baseline['claims'] if c['id']=='sp3-omega-exact-fortuity')
    changed=copy.deepcopy(chosen);changed['id']='synthetic-conflicting-active-claim'
    changed['proposition']['value']['cohomology']=2
    baseline['claims'].append(changed)
    errors=validate_claims(baseline,check_references=False)
    results.append({'name':'contradictory_active_outcome','rejected':any('conflicting active claims' in e for e in errors)})
    baseline=copy.deepcopy(ledger)
    baseline['claims'][0]['certificate']=[{'kind':'git','commit':ledger['baseline'],'path':'research-roadmap/claim-ledger/DOES_NOT_EXIST_negative_control.dat'}]
    errors=validate_claims(baseline)
    results.append({'name':'missing_immutable_certificate','rejected':any('missing native reference' in e for e in errors)})
    baseline=copy.deepcopy(ledger)
    baseline['claims'][0]['certificate']=[{'kind':'git','commit':ledger['baseline'],'path':'README.md','anchor':'missing-negative-control-anchor-1234'}]
    errors=validate_claims(baseline)
    results.append({'name':'missing_native_certificate_anchor','rejected':any('missing native reference anchor' in e for e in errors)})
    baseline=copy.deepcopy(ledger);baseline['claims'].append(copy.deepcopy(baseline['claims'][0]))
    errors=validate_claims(baseline,check_references=False)
    results.append({'name':'duplicate_claim_identifier','rejected':any('duplicate claim id' in e for e in errors)})
    baseline=copy.deepcopy(ledger)
    chosen=next(c for c in baseline['claims'] if c['id']=='sp3-omega-exact-fortuity')
    chosen['charges']['twice_energy']+=1
    errors=validate_claims(baseline,check_references=False)
    results.append({'name':'corrupted_physical_charge_translation','rejected':any('inconsistent native-to-physical' in e for e in errors)})
    changed=copy.deepcopy(frontier)
    item=next(m for r in changed['records'] for m in r['kinematically_admitted_maps'] if m['order']>1)
    item['all_correction_cochain_slots'].pop()
    errors=validate_frontiers(changed,ledger)
    results.append({'name':'omitted_homogeneous_correction_slot','rejected':any('missing correction freedom' in e for e in errors)})
    changed=copy.deepcopy(audit)
    changed['documents'][0]['semantic_review']=None
    blockers=gate_blockers(ledger,changed,frontier)
    results.append({'name':'omitted_semantic_document_review','rejected':any(b['kind']=='semantic_document_review_incomplete' for b in blockers)})
    changed=copy.deepcopy(advertised)
    next(r for r in changed['records'] if r['uses'])['uses'].pop()
    errors=validate_advertised_inputs(changed,ledger)
    results.append({'name':'omitted_advertised_native_occurrence','rejected':any('occurrence coverage mismatch' in e for e in errors)})
    changed=copy.deepcopy(releases)
    changed['release_ids']=[r for r in changed['release_ids'] if r!='b05']
    blockers=gate_blockers(ledger,audit,frontier,changed)
    results.append({'name':'omitted_post_baseline_release','rejected':any(b['kind']=='post_baseline_release_review_missing' for b in blockers)})
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--require-complete',action='store_true')
    parser.add_argument('--negative-controls',action='store_true')
    parser.add_argument('--receipt',type=Path)
    args=parser.parse_args()
    ledger=json.loads((HERE/'claims.json').read_text())
    audit=json.loads((HERE/'source-audit.json').read_text())
    frontier=json.loads((HERE/'native-frontiers.json').read_text())
    coverage=json.loads((HERE/'coverage-blocks.json').read_text())
    advertised=json.loads((HERE/'advertised-native-inputs.json').read_text())
    dedup=json.loads((HERE/'deduplication.json').read_text())
    releases=json.loads((HERE/'release-audit.json').read_text())
    corepath=HERE/'coverage'/'coverage_manifest.json'
    core=json.loads(corepath.read_text()) if corepath.exists() else None
    errors=(validate_claims(ledger)+validate_frontiers(frontier,ledger)+validate_coverage(core,coverage)+validate_source_audit(audit,ledger)
        +validate_advertised_inputs(advertised,ledger)+validate_deduplication(dedup,ledger,audit,releases)+validate_releases(releases,ledger))
    if (HERE/'CURRENT_SUMMARY.md').read_text()!=summary(ledger,audit,frontier,core,advertised,dedup,releases):
        errors.append('generated current summary is stale')
    negatives=negative_controls(ledger,frontier,audit,advertised,releases) if args.negative_controls else []
    errors.extend('negative control did not reject corruption: '+r['name'] for r in negatives if not r['rejected'])
    blockers=gate_blockers(ledger,audit,frontier,releases)
    report={'schema_version':1,'command':'python3 research-roadmap/claim-ledger/check.py'+(' --negative-controls' if args.negative_controls else '')+(' --require-complete' if args.require_complete else ''),
            'ledger_consistent':not errors,'A07_complete':not errors and not blockers,
            'canonical_claim_count':len(ledger['claims']),'errors':errors,
            'gate_blockers':blockers,'negative_controls':negatives,
            'artifact_sha256':{name:digest((HERE/name).read_bytes()) for name in ['claims.json','source-audit.json','native-frontiers.json','advertised-native-inputs.json','deduplication.json','release-audit.json','release-sources.json.gz','coverage-blocks.json','CURRENT_SUMMARY.md','advertisement-spans.jsonl.gz']}}
    if args.receipt:
        args.receipt.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))
    raise SystemExit(1 if errors else 2 if args.require_complete and blockers else 0)


if __name__=='__main__':
    main()
