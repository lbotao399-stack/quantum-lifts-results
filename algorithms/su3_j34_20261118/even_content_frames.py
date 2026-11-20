"""Rational bases and detecting polynomial coefficients of a fixed C-even content.

A row is the archived unnormalized (1+C) applied to a canonical trace word.
Weights with repeated odd letters use exterior powers. The dimension is obtained
from the published-runtime twined character, and its sum with the C-odd dimension
is checked against a separate A2 constant-term calculation. A square minor
nonzero modulo 1000003 certifies independence over Q and injectivity of the
selected monomial coefficients on this entire rational content.
"""
import argparse,random,hashlib,subprocess
from collections import Counter
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from flint import nmod_mat
import physical_parent_basis as base
from physical_image_solve import pivots
from su3_charge_dimensions_v1 import parity_parts

s,f,l=base.s,base.l.f,base.l; P=base.P

@lru_cache(20000)
def reverse_trace(tr):
    n=sum(f.ODD[x] for x in tr);sg,t=f.canon(tr[::-1]);return sg*(-1)**(len(tr)+n*(n-1)//2),t

def conjugate(w):
    a=1;ts=[]
    for tr in w:
        sg,u=reverse_trace(tr);a*=sg;ts.append(u)
    sg,u=f.multi(tuple(ts));assert a*sg in (-1,1);return a*sg,u

def projected(w):
    sg,u=conjugate(w)
    if u==w:return {w:2} if sg==1 else {}
    return {w:1,u:sg}

def parity_frame_dims(g):
    counts=Counter(g);ev,od=parity_parts(*(tuple(sorted((k for i,k in counts.items() if f.ODD[i]==v),reverse=True)) for v in (False,True)));assert ev+od==base.dims.dimension(g);return ev,od

def poly_iter(words,exe,folder):
    # Only one length per invocation. The C++ exporter uses independent ordered
    # Grassmann coordinates; no numerical matrix configurations are sampled.
    assert words;size=sum(map(len,words[0]));assert size in (10,11)
    odd=bytes(i for i,v in enumerate(f.ODD) if v).hex()
    with TemporaryDirectory(dir=folder) as nm:
        d=Path(nm);a,b=d/'in',d/'out';a.write_text(str(len(words))+'\n'+''.join(f"{len(w)} {' '.join(s.hexes(w))}\n" for w in words));args=[str(exe.resolve()),str(a),odd,str(b)];subprocess.run(args+(['input10'] if size==10 else []),check=True)
        with b.open() as file:
            for w,line in zip(words,file,strict=True):
                n,*v=map(int,line.split());assert len(v)==2*n;yield zip(v[::2],v[1::2])

def polys(words,exe,folder):
    ps=[projected(w) for w in words];pool=sorted({v for row in ps for v in row});raw=dict(zip(pool,(dict(r) for r in poly_iter(pool,exe,folder)),strict=True));out=[]
    for p in ps:
        row=Counter()
        for w,a in p.items():
            for k,v in raw[w].items():row[k]+=a*v
        out.append({k:v for k,v in row.items() if v})
    return out

def make(g,exe,folder,seed=340819):
    folder.mkdir(parents=True,exist_ok=True);dim,odd=parity_frame_dims(g);assert dim>0 and len(g) in (10,11) and len(g)-g.count(0)<=7
    words=base.content_parents.all_words(g);assert len(words)==base.control(g);ws=[]
    for w in words:
        sign,v=conjugate(w)
        if w<v or (v==w and sign>0):ws.append(w)
    # Compare the packed involution with the independent symbolic runtime.
    import su_trace_adapter_v1 as a
    for w in ws[::max(1,len(ws)//16)]:
        z=a.parity_project({s.expanded(w):1},1);assert {s.packed(m):n for m,n in z.items()}==projected(w)
    rnd=random.Random(seed);rnd.shuffle(ws);keep=[];rows=[];steps=[];pos=0;keys=[]
    while pos<len(ws):
        trial=ws[pos:pos+max(128,min(2*(dim-len(rows)),max(dim+128,384)))];pos+=len(trial);ks=keep+trial;vs=rows+polys(trial,exe,folder);avail=sorted(set().union(*(r.keys() for r in vs)));rnd.shuffle(avail);preferred=list(dict.fromkeys(keys+avail[:max(2*dim,256)]));A=base.selected_matrix(vs,preferred);rr=pivots(A.T)
        if len(rr)<len(keep) or (len(rr)<dim and pos==len(ws)):
            preferred=avail;A=base.selected_matrix(vs,preferred);rr=pivots(A.T)
        assert len(keep)<=len(rr)<=dim;cc=pivots(A[rr,:]);keys=[preferred[j] for j in cc];keep=[ks[i] for i in rr];rows=[vs[i] for i in rr];steps.append([pos,len(rr)]);print(g.hex(),pos,len(rr),'of',dim,flush=True)
        if len(rr)==dim:
            det=int(nmod_mat((A[np.ix_(rr,cc)]%P).tolist(),P).det());assert det;return dict(content=g.hex(),length=len(g),C_even_dimension=dim,C_odd_dimension=odd,formal_unprojected_words=len(words),formal_even_orbits=len(ws),prime=P,seed=seed,steps=steps,words=[s.hexes(w) for w in keep],polynomial_keys=keys,determinant_mod_prime=det,projector='unnormalized 1+C')
    raise ValueError(('did not reach character dimension',g.hex(),dim,len(keep)))

def check(z,exe,folder):
    g=bytes.fromhex(z['content']);dim,od=parity_frame_dims(g);ws=[s.unpack_hex(w) for w in z['words']];kk=z['polynomial_keys'];assert z['prime']==P and len(set(ws))==len(ws)==len(set(kk))==len(kk)==dim==z['C_even_dimension'] and od==z['C_odd_dimension'];assert all(bytes(sorted(b''.join(w)))==g for w in ws);vs=polys(ws,exe,folder);A=base.selected_matrix(vs,kk);det=int(nmod_mat((A%P).tolist(),P).det());assert det and det==z['determinant_mod_prime'];return det

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--content',required=True);a.add_argument('--exe',type=Path,required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--check',action='store_true');b=a.parse_args();b.work.mkdir(parents=True,exist_ok=True);path=b.work/(b.content+'.json.gz');g=bytes.fromhex(b.content)
    if b.check:print(check(s.read(path),b.exe,b.work))
    else:z=make(g,b.exe,b.work);s.write(path,z)
