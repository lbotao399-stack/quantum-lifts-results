"""Extract a rational spanning set of three-colour parent words by content.

Each word is a polynomial on an arbitrary diagonal scalar and generic other
sl(3) matrices, including independent odd entries. The semisimple regular
locus is dense in the even scalar; its conjugates detect invariant identities.
A square integer minor nonzero modulo the stated prime proves at least its
order in characteristic zero. Accept a spanning claim only upon attaining the
independently computed, unprojected character dimension. If fewer words have
been tried, failure to attain the dimension is expressly inconclusive.
"""
import argparse,hashlib,random,subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from flint import nmod_mat
import content_dimensions as dims
import content_parents
import scalar_filtration as l
from physical_image_solve import pivots

s,t=l.s,l.tri;P=t.PRIME

def export(words,exe,odd,folder):
    with TemporaryDirectory(dir=folder) as z:
        d=Path(z);pin,pout=d/'in',d/'out';pin.write_text(str(len(words))+'\n'+''.join(f"{len(w)} {' '.join(s.hexes(w))}\n" for w in words));subprocess.run([str(exe.resolve()),str(pin),odd,str(pout),'input10'],check=True);rows=[]
        with pout.open() as f:
            for line in f:
                n,*a=map(int,line.split());r=dict(zip(a[::2],a[1::2]));assert len(a)==2*n==2*len(r) and all(r.values());rows.append(r)
        assert len(rows)==len(words);return rows

def selected_matrix(rows,keys):
    idx={k:i for i,k in enumerate(keys)};A=np.zeros((len(rows),len(keys)),dtype=np.int64)
    for i,row in enumerate(rows):
        for k,v in row.items():
            j=idx.get(k)
            if j is not None:A[i,j]=v
    return A

def control(g):
    # Graded Burnside + graded symmetric powers, independent of recursive traces.
    from itertools import product
    from math import gcd,comb,factorial
    from collections import Counter
    counts=Counter(g);ls=sorted(counts);target=tuple(counts[i] for i in ls);zero=tuple(0 for _ in ls);h={zero:1}
    parity=lambda u:sum(n*l.f.ODD[x] for n,x in zip(u,ls))%2
    for u in product(*(range(t+1) for t in target)):
        n=sum(u)
        if n<2:continue
        total=0
        for rot in range(n):
            periods=gcd(rot,n);r=n//periods
            if any(v%r for v in u):continue
            base=tuple(v//r for v in u);mult=factorial(periods)
            for v in base:mult//=factorial(v)
            total+=mult*(-1)**((r-1)*parity(base))
        assert total>=0 and total%n==0;c=total//n
        if not c:continue
        nxt=Counter(h)
        for v,k in h.items():
            m=min((t-a)//b for t,a,b in zip(target,v,u) if b)
            for j in range(1,m+1):
                w=tuple(a+j*b for a,b in zip(v,u));nxt[w]+=k*(comb(c,j) if parity(u) else comb(c+j-1,j))
        h=nxt
    return h.get(target,0)

def frame(g,exe,folder,seed,check=None):
    folder.mkdir(parents=True,exist_ok=True);odd=bytes(i for i,v in enumerate(l.f.ODD) if v).hex();dim=dims.dimension(g);assert dim>0 and len(g)==10 and g.count(0)>=4
    if check:
        z=s.read(check);assert z['content']==g.hex() and z['dimension']==dim and z['prime']==P;ws=[s.unpack_hex(w) for w in z['words']];assert len(set(ws))==dim and all(bytes(sorted(b''.join(w)))==g for w in ws);rows=export(ws,exe,odd,folder);A=selected_matrix(rows,z['polynomial_keys']);det=int(nmod_mat((A%P).tolist(),P).det());assert det and det==z['determinant_mod_prime'];return dict(content=g.hex(),rank=dim,determinant=det)
    pool=list(content_parents.all_words(g));burn=control(g);assert len(pool)==burn,(len(pool),burn);random.Random(seed).shuffle(pool);selected=[];rows=[];pos=0;steps=[]
    while pos<len(pool):
        new=pool[pos:pos+max(192,2*(dim-len(rows)))];pos+=len(new);vs=export(new,exe,odd,folder);words=selected+new;allrows=rows+vs;freq=CounterKey(allrows);ks=sorted(freq,key=lambda k:hashlib.blake2b(k.to_bytes(8,'big'),digest_size=8).digest());take=min(len(ks),max(dim*3,300));A=selected_matrix(allrows,ks[:take]);rl=pivots(A.T);steps.append([pos,take,len(rl)]);assert len(rl)>=len(rows) or take<len(ks)
        # A new sampled column set can in principle miss a previous pivot.
        if len(rl)<len(rows) or (pos==len(pool) and len(rl)<dim):
            A=selected_matrix(allrows,ks);rl=pivots(A.T)
        rows=[allrows[i] for i in rl];selected=[words[i] for i in rl];print(g.hex(),pos,'rank',len(rows),'of',dim,flush=True)
        if len(rows)==dim:
            I=pivots(A[rl,:]);assert len(I)==dim;matrix=A[rl,:][:,I];det=int(nmod_mat((matrix%P).tolist(),P).det());assert det;return dict(content=g.hex(),dimension=dim,prime=P,formal_words=len(pool),independent_formal_count=burn,tested=pos,sampling_seed=seed,samples=steps,determinant_mod_prime=det,words=[s.hexes(w) for w in selected],polynomial_keys=[ks[i] for i in I],scope='Exact full parent span at three colours, before charge-conjugation projection.')
    raise ValueError(('Parent dimension not saturated by sampled rows',g.hex(),len(rows),dim))

def CounterKey(rows):return {k for row in rows for k in row}
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--group',required=True);a.add_argument('--exe',type=Path,required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--check',type=Path);a.add_argument('--seed',type=int,default=340647);b=a.parse_args();g=bytes.fromhex(b.group);z=frame(g,b.exe,b.work,b.seed,b.check);s.write(b.out,z);print('verified' if b.check else 'complete',g.hex(),flush=True)
