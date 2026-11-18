"""Screen j=22 highest weights against Pfaffian products and the j=18 class.

All ranks are over the reported prime. The universal closed trace space is
included first. This is a generator exclusion test, not a quantum lift test.
"""
from hw_scan import *
from protected import eval_products, monomials, generators
from collections import Counter


def family_generators(jmax=22):
    # Two quartic Pfaffians cost at least j=16. Keep reflection-odd
    # generators in the separately tagged Pf family below.
    gs=[g for g in generators(jmax-16) if g['reflection']==0]
    for inds in combinations_with_replacement(range(3),4):
        s,w=canon(tuple((i,) for i in inds),())
        for n in range(4):
            for acts in combinations(range(3),n):
                e={w:s}
                for i in acts:e=action_expr(e,str(i))
                if e:
                    gs.append(dict(name=('Pf',inds,acts),expr=e,sector=expr_sec(e),
                                   odd=n%2,reflection=1,family='Pf'))
    core=deserialize(json.loads((ROOT/'so8_j18_representative.json').read_text())['terms'])
    for n in range(4):
        for acts in combinations(range(3),n):
            e=core
            for i in acts:e=action_expr(e,str(i))
            if e:
                gs.append(dict(name=('O18',acts),expr=e,sector=expr_sec(e),
                               odd=(n+1)%2,reflection=0,family='O18'))
    for g in gs:g.setdefault('family','trace')
    return gs


def expr_sec(e):
    secs={charges_of_content(count(w)) for w in e}
    assert len(secs)==1
    return secs.pop()


class Raising:
    def __init__(self,gs,p):
        self.gs,self.p=gs,p
        self.solvers={};self.cache={}

    def solve(self,e,sector,family):
        p=self.p;key=(sector,family)
        if key not in self.solvers:
            piv={}
            for k,g in enumerate(self.gs):
                if g['sector']!=sector or g['family']!=family:continue
                v={w:c%p for w,c in g['expr'].items() if c%p};recipe={k:1}
                while v:
                    w=min(v);c=v[w]
                    if w not in piv:
                        inv=pow(c,-1,p)
                        piv[w]=({u:a*inv%p for u,a in v.items()},
                                {u:a*inv%p for u,a in recipe.items()})
                        break
                    row,r=piv[w]
                    for u,a in row.items():
                        v[u]=(v.get(u,0)-c*a)%p
                        if not v[u]:del v[u]
                    for u,a in r.items():
                        recipe[u]=(recipe.get(u,0)-c*a)%p
                        if not recipe[u]:del recipe[u]
            self.solvers[key]=piv
        v={w:c%p for w,c in e.items() if c%p};answer={}
        while v:
            w=min(v);c=v[w]
            if w not in self.solvers[key]:
                raise AssertionError(('generator raising not spanned',sector,family,w))
            row,r=self.solvers[key][w]
            for u,a in row.items():
                v[u]=(v.get(u,0)-c*a)%p
                if not v[u]:del v[u]
            for u,a in r.items():answer[u]=(answer.get(u,0)+c*a)%p
        return {u:a for u,a in answer.items() if a}

    def generator(self,k,i):
        if (k,i) in self.cache:return self.cache[k,i]
        g=self.gs[k];e=defaultdict(int)
        for w,c in g['expr'].items():
            for v,d in raising(w,i).items():e[v]+=c*d
        e={w:c for w,c in e.items() if c}
        if not e:result={}
        else:result=self.solve(e,expr_sec(e),g['family'])
        self.cache[k,i]=result
        return result

    def product(self,m,i):
        out=defaultdict(int)
        for pos,k in enumerate(m):
            for target,c in self.generator(k,i).items():
                new=m[:pos]+(target,)+m[pos+1:]
                odds=[k for k in new if self.gs[k]['odd']]
                if len(odds)!=len(set(odds)):continue
                sign=(-1)**sum(odds[j]>odds[l] for j in range(len(odds)) for l in range(j+1,len(odds)))
                out[tuple(sorted(new))]+=sign*c
        return {w:c%self.p for w,c in out.items() if c%self.p}


def screen(sec,p,census_dir=None):
    level=2*sum(sec[1:])
    started=time.monotonic();gs=family_generators(level);raising_maps=Raising(gs,p)
    ms=[m for m in monomials(gs,sec,0)
        if any(gs[k]['family'] in ('Pf','O18') for k in m)]
    images=[{(i,m):c for i in (1,2) for m,c in raising_maps.product(m0,i).items()} for m0 in ms]
    P,_=formal_closed_kernel(images,p)
    trace_basis=basis(sec,0)
    universal,_=formal_closed_kernel([action_expr(e) for e in trace_basis],p)
    exact=exact_hw(sec,0);parts=[];uparts=[];controls={}
    words=sorted({w for e in trace_basis for w in e})
    degrees={w:sum(count(w)[3:6]) for w in words}
    for d,dim in sorted(exact.items()):
        wd=[w for w in words if degrees[w]==d];ix={w:i for i,w in enumerate(wd)}
        rr=[];cc=[];vv=[]
        for i,e in enumerate(trace_basis):
            for w,c in e.items():
                if w in ix:rr.append(i);cc.append(ix[w]);vv.append(c%p)
        A=csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(trace_basis),len(wd)))
        blocks=[];gblocks=[];used=0;k=0 if not d else d+2;perpoint=comb(k,d)
        for attempt in range(4):
            target=dim+12 if not attempt else 2*used
            while used<target:
                npts=min(max(1,32//perpoint),8,ceil((target-used)/perpoint))
                fr=Frame(d,npts,23001+100003*d+7919*used,p,k=k)
                blocks.append((A@fr.values(wd)%p).astype(np.int32))
                gblocks.append((P@eval_products(gs,ms,fr,d)%p).astype(np.int32))
                used+=fr.ncols
            E=np.concatenate(blocks,axis=1).astype(np.int64)
            if rank_mod(E,p)==dim:break
        assert rank_mod(E,p)==dim
        G=np.concatenate(gblocks,axis=1).astype(np.int64)
        assert rank_mod(np.concatenate([E,G]),p)==dim
        controls[d]={'exact_dimension':dim,'faithful_evaluation_rank':dim,'columns':used}
        uparts.append(universal@E%p);parts.append(G)
    U=np.concatenate(uparts,axis=1);G=np.concatenate(parts,axis=1)
    ru=rank_mod(U,p);rg=rank_mod(np.concatenate([U,G]),p)
    path=(ROOT/'hw_runs' if census_dir is None else pathlib.Path(census_dir))/('j%d_n%d_q%d_%d_%d_r0_p%d.json'%((level,)+sec+(p,)))
    census=json.loads(path.read_text())
    assert ru==census['universal_closed_evaluated_rank']
    kernel=sum(exact.values())-census['rank_out']
    assert rg<=kernel
    return {'sector':sec,'prime':p,'family_product_words':len(ms),
            'formal_highest_weight_products':len(P),'classical_kernel_dimension':kernel,
            'universal_closed_trace_rank':ru,'rank_with_Pfaffian_and_O18_products':rg,
            'unexplained_classical_kernel':kernel-rg,'controls':controls,
            'seconds':time.monotonic()-started,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Finite-field generator screening. Zero residual excludes a new generator modulo the recorded families over this prime only.'}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--sector',required=True)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--output',required=True)
    ap.add_argument('--seconds',type=int,default=120)
    a=ap.parse_args();signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('wall limit')))
    signal.alarm(a.seconds)
    r=screen(tuple(map(int,a.sector.split(','))),a.prime)
    pathlib.Path(a.output).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
