"""Antisymmetric derivative exchange in the SU(3) first-correction equation.

Check the complete signed RHS, including cyclic and graded product signs.
At bidegree (1,1), a fixed content either has one (1,1) letter, or one
(1,0) letter and its (0,1) partner. In the latter case, on the two adjoint
colour slots exchange is ordinary transposition multiplied by the common
letter's Koszul sign. Count its two eigenspaces using the existing
three-colour symmetric/exterior character (including charge conjugation).
"""
import argparse,hashlib
from collections import Counter
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory
import replay_guarded as v
import parent_frontier as front
from su3_charge_dimensions_v1 import parity_parts
r,s,f=v.r,v.s,v.r.f
TABLE=bytes([f.ID[b,a,fl] for a,b,fl in f.ALPH]+list(range(len(f.ALPH),256)))
def flip(t):return bytes(sorted(bytes.fromhex(t).translate(TABLE))).hex()
@lru_cache(160000)
def trace(t):return f.canon(t.translate(TABLE))
def word(w):
    ts=[];sg=1
    for t in w:
        a,n=trace(t);sg*=a;ts.append(n)
    a,z=f.multi(tuple(ts));assert a*sg in (1,-1);return a*sg,z
def exchange_vector(x):
    out={}
    for w,c in x.items():
        a,z=word(w);assert z not in out;out[z]=a*c
    return out
def control():
    for i,ss in enumerate(f.SPLITS):
        assert TABLE[TABLE[i]]==i and f.ODD[TABLE[i]]==f.ODD[i]
        assert sorted((c,p.translate(TABLE)) for c,p in ss)==sorted(f.SPLITS[TABLE[i]])
def subspaces(g):
    dim=r.E.parity_frame_dims(g)[0];c=Counter(g);moved=[k for k in c if TABLE[k]!=k]
    if not moved:return 0,dim,'one_mixed_derivative'
    i,j=sorted(moved);assert TABLE[i]==j and TABLE[j]==i and c[i]==c[j]==1 and f.ODD[i]==f.ODD[j]
    del c[i],c[j]
    vals=[]
    for anti in (True,False):
        bosons=[n for k,n in c.items() if not f.ODD[k]];fermions=[n for k,n in c.items() if f.ODD[k]]
        take_exterior=anti ^ bool(f.ODD[i]);(fermions if take_exterior else bosons).append(2)
        vals.append(parity_parts(tuple(sorted(bosons,reverse=True)),tuple(sorted(fermions,reverse=True)))[0])
    assert sum(vals)==dim
    return *vals,'two_exchanged_derivatives'
def dimensions(contents):
    universe={g.hex() for g in contents};assert {flip(t) for t in universe}==universe
    pairs=Counter();cat=Counter();dim=Counter();full=0;even=0;odd=0
    for t in sorted(universe):
        g=bytes.fromhex(t);n=r.E.parity_frame_dims(g)[0];full+=n;u=flip(t)
        if u<t:continue
        if u==t:
            minus,plus,k=subspaces(g);cat[k]+=1;dim[k]+=minus
        else:
            assert n==r.E.parity_frame_dims(bytes.fromhex(u))[0];minus=plus=n;cat['two_contents']+=1;dim['two_contents']+=minus
        odd+=minus;even+=plus
        if minus:pairs['nonempty_negative']+=1
    assert odd+even==full
    return dict(full_contents=len(universe),orbits=dict(cat),negative_dimension_by_orbit_type=dict(dim),active_negative_orbits=pairs['nonempty_negative'],original_dimension=full,positive_dimension=even,negative_dimension=odd)
def rhs_eigen(work):
    h=hashlib.sha256();n=0
    with TemporaryDirectory(dir=work) as tmp:
        files=[(Path(tmp)/str(i)).open('w') for i in range(32)]
        def emit(w,c,k):
            t=bytes(sorted(b''.join(w)));i=int.from_bytes(hashlib.sha256(t).digest()[:2],'big')%32
            files[i].write(f'{k} {c} '+','.join(s.hexes(w))+'\n')
        for w,c in s.stream_rows(v.HERE/'negative_Q1.jsonl.gz',h):
            sg,z=word(w);u,ww=word(z);assert sg*u==1 and ww==w;n+=1;num=f.numerator(c);emit(w,num,0);emit(z,num*sg,1)
        for stream in files:stream.close()
        assert h.hexdigest()==s.INPUT and n==517468;bad=0
        for i in range(32):
            b=[{},{}]
            with (Path(tmp)/str(i)).open() as st:
                for line in st:
                    q,c,k=line.split();d=b[int(q)];assert k not in d;d[k]=int(c)
            bad+=sum(b[0].get(k,0)!=-b[1].get(k,0) for k in set(b[0])|set(b[1]))
        assert bad==0;return dict(uncompressed_sha256=s.INPUT,formal_words=n,negative_exchange_mismatches=bad)
def run(work,out):
    work.mkdir(exist_ok=True,parents=True);control();F,B=front.incidence();assert all({flip(t.hex()) for t in F[g]}=={t.hex() for t in F[bytes.fromhex(flip(g.hex()))]} for g in F)
    z=dict(rhs=rhs_eigen(work),commuting_Q0_letters=len(f.ALPH),verified_split_terms=sum(map(len,f.SPLITS)),length_10_parents=dimensions(F),length_11_equations=dimensions(B));s.write(out,z);print(z)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--work',required=True,type=Path);a.add_argument('--out',required=True,type=Path);assert __debug__;run(**vars(a.parse_args()))
