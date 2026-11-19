"""Exact, rationally liftable tests on the *entire* allowed parent subspace per content.

For each edge of the split graph, check that all pairs of integral coefficient
matrices commute BEFORE reducing modulo p. A nonzero Grassmann coordinate thus
separates the rational residual from every allowed image, without any modular
rank assumption. No conclusion follows from a zero. Ascending support costs
are the number of previously successful contents a parent can reach by splits.
"""
import argparse, hashlib, json, subprocess
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory
import commuting_parent_test as old
s,ft=old.s,old.ft;P=1000003;D=7

def graph(census):
    doc=s.read(census);ok=set(doc['exact_zero']);bad=set(doc['mathematically_nonzero']);assert not doc['unresolved'] and len(ok)==308 and len(bad)==192 and not bad&ok
    back={t:old.split_pairs(t) for t in sorted(ok|bad)};fw=defaultdict(set)
    for t,ps in back.items():
        for g,(k,i,j,_) in ps.items():fw[g].add(t)
    # Independently enumerate the whole outgoing support of every reconstructed parent.
    for g,ts in fw.items():
        src=bytes.fromhex(g);by_forward={bytes(sorted(src[:h]+src[h+1:]+pair)).hex() for h,k in enumerate(src) for _,pair in old.f.SPLITS[k]};assert ts==by_forward
    assert len(fw)==581 and sum(map(len,fw.values()))==3985
    cost={g:len(v&ok) for g,v in fw.items()};return bad,back,cost,doc

def adjsetup(edges):
    adj=defaultdict(set)
    for i,j in edges:adj[i].add(j);adj[j].add(i)
    comps=[];left=set(adj)
    while left:
        part={min(left)};todo=list(part)
        while todo:
            for j in adj[todo.pop()]-part:part.add(j);todo.append(j)
        comps.append(sorted(part));left-=part
    return adj,comps

def mods(edges):
    adj,comps=adjsetup(edges);vertices=sorted(adj)
    # Complements of vertex covers: block matrices may vary independently only on an independent set.
    # Use all maximal independent sets; for disconnected graphs choose one in each component.
    out=[('diagonal',())]
    for n in range(1,len(vertices)+1):
        for sel in combinations(vertices,n):
            S=set(sel)
            if any(adj[i]&S for i in sel) or any(not adj[i]&S for i in set(vertices)-S):continue
            out.append(('block',sel))
    return out

def prod(A,B):return tuple(tuple(sum(x*y for x,y in zip(r,c)) for c in zip(*B)) for r in A)
def transpose(A):return tuple(zip(*A))
def liftmat(r,ids,edges,mode,selected):
    adj,comps=adjsetup(edges); M={};colors={}
    for cc in comps:
        # Integral change of basis (products of determinant-one integer shears).
        U=V=((1,0,0),(0,1,0),(0,0,1))
        for j in range(6):
            i,k=r.sample(range(3),2);a=r.choice([-2,-1,1,2]);add=lambda n:tuple(tuple(int(x==y)+(n if (x,y)==(i,k) else 0) for y in range(3)) for x in range(3));U=prod(U,add(a));V=prod(add(-a),V)
        assert prod(U,V)==((1,0,0),(0,1,0),(0,0,1))
        for k in cc:colors[k]=(U,V)
    for k in ids:
        for h in (range(D) if ft.ODD[k] else [-1]):
            a=ft.traceless(r,3)
            if k in adj:
                if mode=='diagonal':a=tuple(tuple(y if i==j else 0 for j,y in enumerate(row)) for i,row in enumerate(a))
                elif k in selected:a=tuple(tuple(y if (i==2)==(j==2) else 0 for j,y in enumerate(row)) for i,row in enumerate(a))
                else:
                    x=r.choice([-4,-3,-2,-1,1,2,3,4]);a=((x,0,0),(0,x,0),(0,0,-2*x))
                U,V=colors[k];a=prod(prod(U,a),V)
            assert a[0][0]+a[1][1]+a[2][2]==0;M[k,h]=a
    def parts(i):return [M[i,h] for h in (range(D) if ft.ODD[i] else [-1])]
    for i,j in edges:
        for A in parts(i):
            for B in parts(j):assert prod(A,B)==prod(B,A)
    return M

def seedof(t,edges,mode,sel,seed):
    # Use a portable integer seed; the actual witness stores every integer matrix as well.
    return int.from_bytes(hashlib.sha256(json.dumps([t,list(edges),mode,list(sel),seed],separators=(',',':')).encode()).digest()[:12],'big')

def prepare(rhs,shared,census,folder):
    bad,b,c,doc=graph(census);n=s.read(shared);den=int(n['denominator']);assert den>0 and den%P;v,hv=old.vector(shared);assert hv==doc['vector_sha'];rows=defaultdict(dict)
    for w,z in v:
        for a,m in old.Q_multi(s.expanded(w)).items():
            x=s.packed(a);g=bytes(sorted(b''.join(x))).hex();
            if g in bad:rows[g][x]=rows[g].get(x,0)+m*z
    folder.mkdir(parents=True,exist_ok=True); paths={g:folder/(g+'.body') for g in bad};out={g:v.open('w') for g,v in paths.items()};counts=Counter();h=hashlib.sha256();n=0;inv={}
    def emit(g,w,v):
        if not v:return
        q=Fraction(v);assert q.denominator%P;iv=inv.get(q.denominator)
        if iv is None:iv=pow(q.denominator,-1,P);inv[q.denominator]=iv
        out[g].write(f'{q.numerator*iv%P} {len(w)} '+ ' '.join(s.hexes(w))+'\n');counts[g]+=1
    try:
        for w,z in s.stream_rows(rhs,h):
            n+=1;g=bytes(sorted(b''.join(w))).hex()
            if g in bad:emit(g,w,rows[g].pop(w,0)-Fraction(z))
        for g,row in rows.items():
            for w,z in row.items():emit(g,w,z)
    finally:
        for f in out.values():f.close()
    assert n==517468 and h.hexdigest()==s.INPUT==doc['rhs_sha'] and set(counts)==bad
    for g,p in paths.items():
        (folder/(g+'.words')).write_bytes(f'{sum(ft.ODD[k] for k in bytes.fromhex(g))} {counts[g]}\n'.encode()+p.read_bytes());p.unlink()
    return hv,counts

def one(t,costs,parents,exe,folder,seeds,level):
    ids=sorted(set(bytes.fromhex(t)));edges=tuple(sorted({(i,j) for g,(k,i,j,v) in parents.items() if costs[g]<=level}));allow=sorted(g for g in parents if costs[g]<=level);word=folder/(t+'.words');odd=bytes(i for i,v in enumerate(ft.ODD) if v).hex()
    trace={};trials={};modes=mods(edges)
    # Only write the much smaller remaining assignments if the previous one vanishes.
    for seed in seeds:
        for mode,sel in modes:
            ss=seedof(t,edges,mode,sel,seed);m=liftmat(Random(ss),ids,edges,mode,sel);key=f'{seed}_{mode}_'+bytes(sel).hex();buf=[f'{key} {len(m)}']
            for (k,h),A in sorted(m.items()):buf.append(f'{k} {h} '+' '.join(str(x) for row in A for x in row))
            raw=('\n'.join(buf)+'\n').encode();dg=hashlib.sha256(raw).hexdigest();trials[key]=dg
            with TemporaryDirectory(dir=folder) as td:
                path=Path(td)/'mat';path.write_bytes(raw);out=subprocess.check_output([str(exe),str(word),str(path),odd],text=True).strip()
            if out=='none':continue
            test,mask,val=out.split();assert test==key and int(mask).bit_count()==sum(ft.ODD[k] for k in bytes.fromhex(t)) and 0<int(val)<P
            return dict(allowed=allow,edges=edges,trials=trials,witness=dict(mode=mode,selected=list(sel),seed=seed,mask=int(mask),value=int(val),matrices={f'{k:02x}:{h}':A for (k,h),A in sorted(m.items())}))
    return dict(allowed=allow,edges=edges,trials=trials,witness=None)

def run(rhs,shared,census,exe,folder,out,workers,seeds,levels):
    bad,back,cost,doc=graph(census);hv,counts=prepare(rhs,shared,census,folder);targets=sorted(bad);expo=exe.resolve();ls=sorted(set(levels));assert ls and min(ls)>=0 and max(ls)<=17 and len(seeds)==len(set(seeds))
    def task(t):
        d={};seen={}
        for level in ls:
            v=tuple(sorted(g for g in back[t] if cost[g]<=level))
            if v in seen:d[str(level)]={'same_as':seen[v]};continue
            seen[v]=str(level);d[str(level)]=one(t,cost,back[t],expo,folder,seeds,level)
        return t,dict(formal_words=counts[t],parent_success_counts={g:cost[g] for g in sorted(back[t])},tests=d)
    z={};
    with ThreadPoolExecutor(workers) as pool:
        for i,(t,r) in enumerate(pool.map(task,targets)):
            z[t]=r
            if i%16==15: print('finished',i+1,'/',len(targets),flush=True)
    summary={}
    for l in ls:
        def entry(x):
            d=x['tests'][str(l)];return x['tests'][d['same_as']] if 'same_as' in d else d
        rows=[entry(v) for v in z.values()];summary[str(l)]=dict(eligible_parents=sum(v<=l for v in cost.values()),unreachable=sum(not r['allowed'] for r in rows),witnesses=sum(r['witness'] is not None for r in rows),tried_without_witness=sum(r['witness'] is None for r in rows))
    result=dict(rhs_sha=s.INPUT,source_sha=hv,prime=P,exterior_generators=D,seeds=seeds,levels=ls,successes=308,failures=192,summary=summary,targets=z,sign='Q0 O1 + Q1 O0',scope='Nonzero evaluations of the rational residual with the full tested parent image zero over Z. Zero evaluations give no image membership claim. Costs count successful contents reachable by a single parent.')
    s.write(out,result);return result

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--shared',type=Path,default=s.HERE/'joint_70/vector_shared.json.gz');a.add_argument('--census',type=Path,default=s.HERE/'joint_70/census_2724.json');a.add_argument('--exe',type=Path,required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--workers',type=int,default=4);a.add_argument('--seeds',type=int,nargs='+',default=[817234,419071]);a.add_argument('--levels',type=int,nargs='+',default=[0,1,2,17]);v=a.parse_args();assert __debug__;print(run(v.rhs,v.shared,v.census,v.exe,v.work,v.out,v.workers,v.seeds,v.levels)['summary'])
