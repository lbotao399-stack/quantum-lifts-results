"""Independently audit the support nesting and replay the strongest witnesses over Q.

The unrestricted test may have no witness. This does not imply that a primitive exists.
The integer check expands exterior coefficients as signed permutations, with no C++ or modular polynomial evaluator.
"""
import argparse, hashlib
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import commuting_parent_test as old
s,ft=old.s,old.ft

def test(row,level):
    r=row['tests'][str(level)]
    return row['tests'][r['same_as']] if 'same_as' in r else r

def run(rhs,shared,census,sweep,receipt):
    job=s.read(receipt); z=s.read(sweep); raw=sweep.read_bytes();assert len(raw)==job['size'] and hashlib.sha256(raw).hexdigest()==job['sha256'] and job['summary']==z['summary'] and job['job']=='5665220'
    cc=s.read(census); p=z['prime'];d=z['exterior_generators'];passed=set(cc['exact_zero']);failed=set(cc['mathematically_nonzero']);assert not cc['unresolved'] and set(z['targets'])==failed and len(failed)==192 and len(passed)==308 and p==1000003 and d==7
    splits={t:old.split_pairs(t) for t in failed|passed}; forward=defaultdict(set)
    for t,ps in splits.items():
        for g in ps:forward[g].add(t)
    costs={g:len(v&passed) for g,v in forward.items()};assert len(costs)==581
    for g,targets in forward.items():
        v=bytes.fromhex(g); assert targets=={bytes(sorted(v[:i]+v[i+1:]+pair)).hex() for i,k in enumerate(v) for _,pair in old.f.SPLITS[k]}
    records={}; best=-1
    for level in z['levels']:
        count=zeros=unreached=0
        for t,x in z['targets'].items():
            r=test(x,level); gs=sorted(g for g in splits[t] if costs[g]<=level);assert r['allowed']==gs and x['parent_success_counts']=={g:costs[g] for g in splits[t]}
            edges=sorted({splits[t][g][1:3] for g in gs}); assert edges==list(map(tuple,r['edges']))
            unreached+=not gs
            if r['witness'] is None:zeros+=1;continue
            count+=1
            if level>best:best=level;records={}
            if level==best:records[t]=r
        assert dict(witnesses=count,unreachable=unreached,tried_without_witness=zeros,eligible_parents=sum(c<=level for c in costs.values()))==z['summary'][str(level)]
    assert best==2 and len(records)==3
    hv,rows=old.residuals(rhs,shared,records);assert z['source_sha']==cc['vector_sha']==hv and z['rhs_sha']==cc['rhs_sha']==s.INPUT
    out={}
    for t,r in sorted(records.items()):
        w=r['witness'];ids=bytes.fromhex(t);bits=[k for k in range(d) if w['mask']>>k&1]; n=len(bits);assert n==sum(ft.ODD[k] for k in ids) and w['mask']<(1<<d)
        parts={(k,h):w['matrices'][f'{k:02x}:{h}'] for k in set(ids) for h in (range(d) if ft.ODD[k] else [-1])};assert len(parts)==len(w['matrices']) and all(sum(A[i][i] for i in range(3))==0 for A in parts.values()); checks=0
        for i,j in r['edges']:
            for (k,h),A in parts.items():
                if k==i:
                    for (l,h2),B in parts.items():
                        if l==j:assert ft.mm(A,B)==ft.mm(B,A);checks+=1
        M={(k,h):parts[k,bits[h] if ft.ODD[k] else -1] for k in set(ids) for h in (range(n) if ft.ODD[k] else [-1])}; value=sum((v*ft.exterior_value(a,M,n) for a,v in rows[t].items()),Fraction());assert 0<w['value']==value.numerator*pow(value.denominator,-1,p)%p and len(rows[t])==z['targets'][t]['formal_words']
        out[t]=dict(allowed=r['allowed'],edges=r['edges'],edge_matrix_checks=checks,words=len(rows[t]),mask=w['mask'],rational=str(value),mod_prime=w['value'],next_parents={g:sorted(forward[g]&passed) for g in splits[t] if costs[g]==3})
    return dict(sha256=job['sha256'],rhs=s.INPUT,source=hv,best_excluding_level=best,all_parent_witnesses=z['summary']['17']['witnesses'],replayed=out,sign=z['sign'],scope='Any repair must use some parent touching at least three previously passing contents. Cancellations on those contents are allowed; no unrestricted obstruction follows.')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--directory',type=Path,default=s.HERE/'joint_70');a.add_argument('--out',type=Path,required=True);v=a.parse_args();q=v.directory;assert __debug__;r=run(v.rhs,q/'vector_shared.json.gz',q/'census_2724.json',q/'centralizer_sweep.json.gz',q/'centralizer_receipt.json');s.write(v.out,r);print([(k,x['rational']) for k,x in r['replayed'].items()])
