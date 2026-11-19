"""Test the first equation with evaluation functionals killing every restricted incoming column.

For target 10 a second family solves both zero- and one-jet equations on the scalar/fermion letters.
Only use a nonzero modular evaluation as rational evidence when the supplied integer matrices annihilate the classical map exactly.
"""
import argparse,gzip,hashlib,json,sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
import projected_q0 as maps
from coupled_first_replay import PRIME,RHS_SHA,coeff
from coupled_kernel_fit import one_frame

def kind(l):return len(l[2]),l[0]+l[1]
def pattern(ws):return tuple(sorted(Counter(kind(l) for l in ws).items()))
def matching(whole,t,s):
    z=Counter(dict(s));z.subtract(dict(t));keep=[]
    for source in whole:
        rules=[]
        for c,l,r in Q_letter(source):
            delta=z.copy();delta[kind(source)]-=1;delta[kind(l)]+=1;delta[kind(r)]+=1
            if not any(delta.values()):rules.append((c,l,r))
        if rules:keep.append((source,rules))
    return keep

def int_frame(alphabet,degree,tag,n,seed):
    # At target 10 the two missing equations are the (two different) commuting-scalar Koszul conditions on psi and its jets.
    f=Frame(3,n,degree+2,[degree],alphabet,p=PRIME,seed=seed)
    rng=np.random.default_rng(seed); diag=lambda a:np.diag([int(a[0]),int(a[1]),-int(a[0]+a[1])]);k=degree+2
    def rnd(shape):return rng.integers(-4,5,shape,dtype=np.int64)
    def traceless():
        M=rnd((3,3));M[-1,-1]=-M[0,0]-M[1,1];return M
    v={l:np.stack([[traceless() for _ in range(k if is_ferm(l) else 1)] for _ in range(n)]) for l in alphabet}
    for p in range(n):
        for l in alphabet:
            if kind(l)[0]==1:
                mat=diag(rnd(2));v[l][p,0]=mat
        if tag=='diagonal_scalars':continue
        if tag=='nilpotent_scalars':
            E=np.zeros((3,3),dtype=np.int64);E[0,1]=1;G=np.zeros((3,3),dtype=np.int64);G[0,2]=1
            # Commuting base and jets tangent to a common conjugation orbit.
            for i in range(3):
                l=(0,0,(i,));a,b=rnd(2);v[l][p,0]=a*E+b*G
            for j in [(0,1),(1,0)]:
                H=traceless()
                for i in range(3):
                    l=(*j,(i,));a,b=rnd(2);M=v[(0,0,(i,))][p,0];v[l][p,0]=a*E+b*G+H@M-M@H
            continue
        raise ValueError(tag)
    for l,x in v.items():f.L[l]=(int(is_ferm(l)),x%PRIME)
    return f,v

def calc(a):
    global Q_letter,is_ferm,Frame
    sys.path.insert(0,str(a.runtime.resolve()));import sp_environment
    from fortuity.full.letters import Q_letter,is_ferm
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents
    from native_su_values_v1 import evaluate
    root=a.root;plan=json.loads((root/'classical_order1.json').read_bytes())['incoming'];mu=tuple(plan['mu']);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});pats={maps.pattern(plan['targets'][t]):t for t in (6,10)};tags=('diagonal_scalars','nilpotent_scalars');frames={};refs={};ranks={};facts={};sums={};direct=Counter()
    for t in (6,10):
        tp=maps.pattern(plan['targets'][t]);deg=plan['targets'][t]['fermions'];sources=[s for s in plan['sources'] if t in s['targets']]
        for it,tag in enumerate(tags):
            seed=641007+3391*t+11*it;fr,integer=int_frame(alpha,deg,tag,a.points,seed);key=t,tag;frames[key]=fr;refs[key]=one_frame(fr);sums[key]=np.zeros(fr.ncols,dtype=np.int64)
            restricted=[];res={}
            for s in sources:
                rr=matching(alpha,tp,maps.pattern(s));assert rr
                cnt=0
                for l,rules in rr:
                    # All edges of the restricted space reaching these targets are boson-boson replacements. Others are recorded only as possible incoming support.
                    if all(not is_ferm(u) and not is_ferm(v) for _,u,v in rules):
                        M=sum(c*(integer[u][:,0]@integer[v][:,0]) for c,u,v in rules);cnt+=np.count_nonzero(M)
                    else:cnt=None;break
                if cnt==0:restricted.append(s['index'])
                res[str(s['index'])]=None if cnt is None else int(cnt)
            assert restricted==({6:[12,13],10:[16,17]}[t]);assert set(restricted)=={s['index'] for s in sources if 12<=s['index']<=22}
            nint=max(int(np.max(np.abs(x))) for x in integer.values());assert 3*nint**2*len(alpha)**2<2**63
            facts[key]=dict(target=t,frame=tag,points=a.points,seed=seed,prime=PRIME,integer_zero_on_sources=restricted,remaining_possible_sources=[s['index'] for s in sources if s['index'] not in restricted],matrices_sha256=hashlib.sha256(b''.join(integer[l].astype('<i8').tobytes() for l in alpha)).hexdigest(),max_integer_entry=nint)
    b={6:{},10:{}};ns=Counter();h=hashlib.sha256();total=0
    def update(t):
        for tag in tags:
            key=t,tag
            if direct[key]<1:
                small=refs[key];w=next(iter(b[t]));assert np.array_equal(small.eval_words([w],dtype=np.int64)[0]*b[t][w]%PRIME,evaluate([{w:b[t][w]}],small)[0]);direct[key]+=1
            sums[key]=(sums[key]+evaluate([b[t]],frames[key])[0])%PRIME
        b[t]={}
    with gzip.open(root/'negative_Q1.jsonl.gz','rb') as compressed:
        for line in compressed:
            h.update(line);total+=1;c,ww=json.loads(line);w=maps.decode(ww);t=pats.get(pattern(l for tr in w for l in tr));
            if t is None:continue
            ns[t]+=1;b[t][w]=coeff(c)
            if len(b[t])==512:update(t)
    assert h.hexdigest()==RHS_SHA and total==517468 and ns=={6:31720,10:27780}
    for t in (6,10):
        if b[t]:update(t)
    results=[]
    for key,v in sums.items():
        z=np.flatnonzero(v);facts[key].update(rhs_terms=ns[key[0]],width=len(v),rhs_nonzero=len(z),rhs_first_nonzero=None if not len(z) else [int(z[0]),int(v[z[0]])],rhs_values=v.tolist());results.append(facts[key]);print(key,len(z),facts[key]['rhs_first_nonzero'],flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(dict(rhs_stream_sha256=h.hexdigest(),records=results),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=maps.HERE);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--out',type=Path,required=True);p.add_argument('--points',type=int,default=6);calc(p.parse_args())
