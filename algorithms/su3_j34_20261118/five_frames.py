"""Frozen side frames in SU(3), j=34. Same projected trace in every incident target."""
from concurrent.futures import ProcessPoolExecutor
import argparse, json, math, sys
from pathlib import Path
import numpy as np
import four_source_joint as f

P=f.P
SIDE=(2,3,8,11)
POINTS=80
SEED=934715231
QUOTIENT='ef37638c5392cf9f2aa70ed345d8365393a1fc17ddd6ba6e3b5033d95bbb354c'

def recipes(root,first,last):
    p,w,_=f.v.inputs(root,root/'allcoord')
    assert f.fs.fsha(first/'new_words.json.gz')=='b984cb977c1ae09ca97b65194b6cb3b037bb29900971a60a7373b94ebee5b87b'
    assert f.fs.fsha(last/'result.json')==QUOTIENT
    record=json.loads((last/'result.json').read_bytes());assert f.fs.fsha(last/'words.json.gz')==record['word_file_sha256']
    w+=f.read_words(first/'new_words.json.gz'); nt=len(w); w+=f.read_words(last/'words.json.gz')
    assert nt==6848 and len(w)==len(set(w))==7616
    assert not set(map(str,SIDE)).intersection(p['first_rhs_trace_terms'])
    assert all(not set(SIDE).intersection(p['sources'][s]['targets']) for s in (19,20,21,22))
    return p,w

def indices(p):
    cols={}; off=0
    for i,t in enumerate(p['targets']):
        w=12*math.comb(t['fermions']+2,t['fermions']);cols[i]=list(range(off,off+w));off+=w
    assert off==2244 and [len(cols[t]) for t in SIDE]==[36,36,120,120]
    return cols

def init(runtime,p):
    global alg,project,derivative,frames,charge,patterns,tp,segments
    sys.path.insert(0,str(runtime));import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from g2_full_classical_screen_v1 import contents,sector
    mu=tuple(p['mu']);charge=sector(mu,p['source_length']);alpha=sorted({l for n in (p['source_length'],p['target_length']) for c in contents(mu,n) for l in c})
    frames={t:f.v.make_frame(p['targets'][t],alpha,12,88304773+229*t) for t in SIDE};patterns={s:f.v.maps.pattern(p['sources'][s]) for s in f.SOURCES};tp={t:f.v.maps.pattern(p['targets'][t]) for t in SIDE};start=0;segments={}
    for t,g in frames.items():segments[t]=slice(start,start+g.ncols);start+=g.ncols
    assert start==312

def side_block(task):
    s,ix,ww=task;A=np.zeros((len(ww),312),dtype='<u4');support=set(f.ORDER[s]).intersection(SIDE)
    if support:
        ex=[project(alg.parity_project({w:1},1),charge,P) for w in ww]
        for t in support:A[:,segments[t]]=derivative(ex,frames[t],patterns[s],tp[t])
    return ix,A

def sides(runtime,p,words,workers):
    A=np.empty((len(words),312),dtype='<u4');done=np.zeros(len(words),bool)
    with ProcessPoolExecutor(workers,initializer=init,initargs=(runtime,p)) as pool:
        for ix,b in pool.map(side_block,f.tasks_for(words,12)):
            assert not done[ix].any();A[ix]=b;done[ix]=True
    assert done.all();return A

def direct_checks(runtime,p,words,A):
    """Expand the graded differential as trace words before evaluating at one point."""
    init(runtime,p);from sp_full_profiles import profile
    from native_su_values_v1 import evaluate
    hits={}
    for s in f.SOURCES:
        for t in sorted(set(f.ORDER[s]).intersection(SIDE)):
            g=f.span.part(frames[t],0,1);sl=slice(segments[t].start,segments[t].start+g.ncols)
            for i,(src,w) in enumerate(words):
                if src!=s or not A[i,sl].any():continue
                ex=project(alg.parity_project({w:1},1),charge,P);dd={h:int(n)%P for h,n in alg.differential(ex).items() if n%P and profile(tuple(l for tr in h for l in tr))==tp[t]}
                v=evaluate([dd],g)[0];np.testing.assert_array_equal(v,A[i,sl]);hits[f'{s}:{t}']=[i,f.fs.harray(v)];break
            else:raise AssertionError((s,t))
    return hits

def fresh_frame(root,rt,p,words,workers):
    import four_source_refit as r
    def tasks():
        for batch in f.v.batches(list(f.v.rhs_terms(root,p).items()),192):yield 'rhs',None,batch
        yield from f.tasks_for(words,10)
    A=np.empty((len(words),10*POINTS),dtype='<u4');b=np.zeros(A.shape[1],dtype=np.int64);done=np.zeros(len(words),bool)
    with ProcessPoolExecutor(workers,initializer=r.init,initargs=(rt,p,POINTS,SEED)) as pool:
        for kind,ix,out in pool.map(r.block,tasks()):
            if kind=='rhs':b=(b+out)%P
            else:assert not done[ix].any();A[ix]=out;done[ix]=True
    assert done.all() and b.any();return A,f.fs.packed(b)

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    for k in ('root','first','last','runtime','out'):a.add_argument('--'+k,type=Path,required=True)
    a.add_argument('--workers',type=int,default=2);args=a.parse_args();sys.path.insert(0,str(args.runtime.resolve())); p,ww=recipes(args.root,args.first,args.last);A=sides(args.runtime,p,ww,args.workers); ck=direct_checks(args.runtime,p,ww,A); f.fs.write_array(args.out,A); print(json.dumps({'rows':len(A),'width':A.shape[1],'sha256':f.fs.harray(A),'explicit_checks':ck}),flush=True)
