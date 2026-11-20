"""Replay rational first corrections from frame numerators and classify every changed content.

The base is reconstructed from the three previously published word vectors.
A support exclusion inherits a base result. For a changed content a nonzero
exterior specialization proves failure; a zero is counted only after exact
formal cancellation or expansion of the invariant polynomial at three colours.
"""
import argparse, hashlib
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import replay_inhomogeneous_repairs as r
s=r.s;HERE=r.HERE
SINGLES=('0000000001010106060715','0000000001010106060c16','0000000001010106060e14')
PAIR='pair/all'; BASE_DIGEST='c3e5f3a0405b643b605a750d95811863c748c0a3442cffb82945d3e14f58c9fd'

def base315():
    p=HERE/'joint_70/inhomogeneous';cc=r.BASE/'preserve_296.census.json';old=s.read(cc);d=s.read(p/'independent.json');assert r.sha(cc)==d['base_census_sha256'] and old['vector_sha']==d['base_vector_sha'] and not old['unresolved']
    x=r.add_rows(r.load_words(r.BASE/'preserve_296.vector.json.gz'),*(r.load_words(p/t/'delta.json.gz') for t in r.TARGETS));v=d['variants']['combined'];assert r.vector_sha(x)==BASE_DIGEST==v['ordered_vector_sha256'] and not v['lost'] and not v['unresolved']
    ok=set(old['exact_zero'])|set(v['gained']);bad=set(old['mathematically_nonzero'])-set(v['gained']);assert len(ok)==315 and len(bad)==185 and not ok&bad
    return x,ok,bad

def unpack(p):
    z=s.read(p);D=int(z['denominator']);terms=z['entries'];assert D>0 and gcd(D,*[int(t[-1]) for t in terms])==1;seen=set();delta={}
    for g,t,n in terms:
        w=s.unpack_hex(t);val=Fraction(int(n),D);assert val and w not in seen and g==r.ch(w);seen.add(w)
        # The original frame uses the unnormalized 1+C.
        for u,k in r.E.projected(w).items():r.add(delta,u,val*k)
    return delta,len(seen)

def changed(delta):
    parents={r.ch(w) for w in delta};return set().union(*(r.splits(g) for g in parents))

def before(x,relevant):
    h=hashlib.sha256();nr=Counter();b={t:{} for t in relevant}
    for w,c in s.stream_rows(HERE/'negative_Q1.jsonl.gz',h):
        t=r.ch(w);nr[t]+=1
        if t in relevant:assert w not in b[t];b[t][w]=Fraction(c)
    assert sum(nr.values())==517468 and len(nr)==174 and h.hexdigest()==s.INPUT
    im=r.action(x,relevant);return nr,b,{t:r.add_rows(im[t],{w:-v for w,v in b[t].items()}) for t in relevant}

def classify(variants,x,ok,bad,nr,b,rr,columns,exterior,work,full_targets):
    # Each variant is the sum of the base and one stated delta.
    details={};states={};full_polys={};seed=817234
    for name,delta in variants.items():
        affected=changed(delta);cur=r.action(delta,affected);cur={t:r.add_rows(rr[t],a) for t,a in cur.items()};ans={};full_polys[name]={}
        wanted=set(full_targets.get(name,[]));assert wanted<=affected
        for t in sorted(wanted):
            q=r.full({'after':cur[t],'before':rr[t],'rhs':b[t]},columns,work);full_polys[name][t]=q;ans[t]=dict(full_polynomial=q['after'],nonzero=bool(q['after']['nonzero_monomials']))
        pending={(name,t):cur[t] for t in sorted(affected-wanted)};witness=r.witnesses(pending,exterior,work,seed)
        for k,row in pending.items():
            _,t=k
            if k in witness:ans[t]=dict(exterior_nonzero=witness[k],nonzero=True)
            elif not row:ans[t]=dict(formal_terms=0,nonzero=False)
            else:
                assert bytes.fromhex(t).count(0)>=4,('cannot decide with this diagonal scalar',t)
                q=r.full({'residual':row},columns,work)['residual'];ans[t]=dict(full_polynomial=q,nonzero=bool(q['nonzero_monomials']))
        good=(ok-affected)|{t for t in affected if not ans[t]['nonzero']};failed=(bad-affected)|{t for t in affected if ans[t]['nonzero']};assert not good&failed and good|failed==ok|bad
        out=r.add_rows(x,delta);states[name]=(out,good,failed);details[name]=dict(words=len(out),ordered_vector_sha256=r.vector_sha(out),delta_words=len(delta),changed_parent_contents=sorted({r.ch(w) for w in delta}),affected_contents=len(affected),gained=sorted(good-ok),lost=sorted(ok-good),solved=len(good),failed=len(failed),inhomogeneous=len(good&set(nr)),homogeneous=len(good-set(nr)),unresolved=0,evidence=ans)
        print(name,'words',len(out),'solves',len(good),len(good&set(nr)),len(good-set(nr)),'gained/lost',len(good-ok),len(ok-good),flush=True)
    return details,full_polys,states

def read_inputs(files,receipt):
    metadata=s.read(receipt);rec=metadata['receipt'];assert metadata['sacct'].split('|')[2]=='COMPLETED' and metadata['sacct'].endswith('|32|0:0') and rec['base']==BASE_DIGEST and rec['rhs']==s.INPUT
    docs={};ds={}
    for name in (*SINGLES,PAIR):
        a=files/('pair' if name==PAIR else name);pre='all.' if name==PAIR else '';p=a/(pre+'frame_terms.shared.json.gz');rp=a/(pre+'result.json');job=rec['results'][name]
        for e in (rp,p):assert job['files'][str(e.relative_to(files))]==[e.stat().st_size,r.sha(e)]
        d=s.read(rp);x,n=unpack(p);assert job['kind']==d['kind']=='exact_solution' and d['base']==BASE_DIGEST and d['rhs']==s.INPUT and len(x)==d['delta_words'] and n==d['frame_terms']
        ds[name]=x;docs[name]=d
    return rec,ds,docs

def run(files,receipt,columns,exterior,work,out,export=None):
    work.mkdir(exist_ok=True,parents=True);rec,ds,docs=read_inputs(files,receipt);x,ok,bad=base315();T={k:changed(d) for k,d in ds.items()}
    for t in SINGLES:
        assert set(docs[t]['affected_contents'])==T[t] and T[t]<=bad and all(w not in T[t] for w in SINGLES if w!=t)
    assert set(SINGLES).isdisjoint(T[PAIR]) and all(T[t].isdisjoint(docs[PAIR]['targets']) for t in SINGLES)
    singles=r.add_rows(*(ds[t] for t in SINGLES));combined=r.add_rows(singles,ds[PAIR]);variants={'three':singles,'four':ds[PAIR],'seven':combined};un=set.union(*T.values());nr,b,rr=before(x,un);r.probe.controls(exterior,817234,work)
    det,phys,states=classify(variants,x,ok,bad,nr,b,rr,columns,exterior,work,{'three':SINGLES,'four':docs[PAIR]['targets']})
    for t in SINGLES:
        p=phys['three'][t];d=docs[t];assert p['after']==d['physical_after'] and p['after']['nonzero_monomials']==0 and p['before']==d['physical_before']['Q0_O1_plus_Q1_O0'] and p['rhs']==d['physical_before']['negative_Q1'] and p['rhs']['nonzero_monomials']>0
    for t in docs[PAIR]['targets']:assert phys['four'][t]['after']['nonzero_monomials']==0 and (phys['four'][t]['before']['nonzero_monomials']==0)==(t in ok)
    assert set(SINGLES)|set(docs[PAIR]['targets'])<=states['seven'][1]
    if export:
        export.mkdir(parents=True,exist_ok=True)
        for k,(v,good,failed) in states.items():
            s.write(export/(k+'.vector.json.gz'),r.word_rows(v));s.write(export/(k+'.census.json'),dict(vector_sha=det[k]['ordered_vector_sha256'],rhs_sha=s.INPUT,exact_zero=sorted(good),mathematically_nonzero=sorted(failed),unresolved=[]))
    s.write(out,dict(job=rec['job'],rhs_sha=s.INPUT,base_sha=BASE_DIGEST,base_words=len(x),base_successes=len(ok),sign='Q0 O1 + Q1 O0',outside_changed_support='inherits base census',nonzero_probe_prime=r.probe.P,variants=det,full_traceless_three_colour=phys))
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);h=HERE/'joint_70/guarded';a.add_argument('--files',type=Path,default=h);a.add_argument('--receipt',type=Path,default=h/'a8_5665888.json')
    for n in ('columns','exterior','work','out'):a.add_argument('--'+n,type=Path,required=True)
    a.add_argument('--export',type=Path);assert __debug__;run(**vars(a.parse_args()))
