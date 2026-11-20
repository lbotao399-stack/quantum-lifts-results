"""Independent integral spot checks of the exhaustive SU(3) matrix expansion.

Decode each archived leading monomial and extract its coefficient by a different
trace-index path enumeration: enumerate colours at every position, force scalar
edges diagonal and never construct the polynomial or use its C++ monomial keys.
Use one nonzero control from each of the 15 homogeneous contents.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
from verify_source_closure import prepare
from second_action import HERE

def unrank(key,n):
    e=key&15;u=key>>4;ids=[]
    for j in range(n,0,-1):
        a=next(i for i in range(64,j-2,-1) if comb(i,j)<=u);u-=comb(a,j);ids.append(a-(j-1))
    ids=list(reversed(ids));assert not u and all(0<=i<56 for i in ids);return e,tuple(ids)

def pure_trace(tr,axis,odd,needed):
    # (entry codes, x power, y power) -> signed coefficient.
    d=Counter()
    for indices in product(range(3),repeat=len(tr)):
        terms=[((),0,0,1)]
        for i,l in enumerate(tr):
            row,col=indices[i],indices[(i+1)%len(tr)]
            if l==axis:
                if row!=col:terms=[];break
                choices=([(-1,1,0,1)] if row==0 else [(-1,0,1,1)] if row==1 else [(-1,1,0,-1),(-1,0,1,-1)])
            else:
                ent=row*3+col;choices=[(8*l+ent,0,0,1)] if ent!=8 else [(8*l,0,0,-1),(8*l+4,0,0,-1)]
                choices=[z for z in choices if z[0] in needed]
            if not choices:terms=[];break
            terms=[(w+(a,) if a>=0 else w,x+b,y+c,k*v) for w,x,y,k in terms for a,b,c,v in choices]
        for k,x,y,sg in terms:
            od=[i for i in k if i//8 in odd]
            if len(set(od))!=len(od):continue
            sign=(-1)**sum(i>j for h,i in enumerate(od) for j in od[h+1:]);d[tuple(sorted(k)),x,y]+=sg*sign
    return {k:v for k,v in d.items() if v}

def coefficient(w,axis,odd,keys,x,y):
    need=Counter(keys);row={((),0,0):1}
    for tr in w:
        pol=pure_trace(tr,axis,odd,need);nextrow=Counter()
        for (a,i,j),k in row.items():
            for (b,m,n),c in pol.items():
                if i+m>x or j+n>y:continue
                data=a+b
                if Counter(data)-need:continue
                o=[u for u in data if u//8 in odd]
                if len(o)!=len(set(o)):continue
                sign=(-1)**sum(u>v for u in a if u//8 in odd for v in b if v//8 in odd)
                nextrow[tuple(sorted(data)),i+m,j+n]+=k*c*sign
        row={k:v for k,v in nextrow.items() if v}
    return row.get((keys,x,y),0)

def run(source,certificate,receipt,exe,work,out):
    import subprocess
    from verify_source_closure import write_words
    p=Path(certificate);doc=json.loads(p.read_bytes());rec=json.loads(receipt.read_bytes())['receipt'];assert rec['files'][p.name]==[p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()] and doc['complete'] and doc['nonzero_monomials']==0
    labels,os,inputs=prepare(source);odd=set(bytes.fromhex(os));alphabet=sorted(labels,key=labels.get);assert set(doc['contents'])=={g.hex() for g in inputs};work.mkdir(exist_ok=True,parents=True);tests={};total=0
    for g,raw in sorted(inputs.items()):
        z=doc['contents'][g.hex()];axis=z['diagonal_letter'];w=tuple(bytes.fromhex(x) for x in z['control_word']);c=raw[w];pin=work/'control.txt';pout=work/'control.json'
        D=write_words(pin,{w:c},axis,os);subprocess.run([str(exe.resolve()),str(pin),str(pout)],check=True);ctrl=json.loads(pout.read_bytes());assert ctrl['nonzero_monomials']==z['control_nonzero_monomials'] and ctrl['first_integer_coefficient']==z['control_first_integer_coefficient']
        x,ks=unrank(ctrl['first_key'],12-g.count(axis));y=g.count(axis)-x;poly=coefficient(w,axis,odd,ks,x,y);assert int(c*D)*poly==ctrl['first_integer_coefficient']!=0;total+=len(raw)
        tests[g.hex()]=dict(control_packed_key=ctrl['first_key'],entry_codes=ks,diagonal_powers=[x,y],coefficient_without_word_weight=poly,cleared_word_coefficient=int(c*D),integer_control_coefficient=ctrl['first_integer_coefficient'])
    assert total==doc['tested_words']==25299 and len(tests)==15
    out.write_text(json.dumps(dict(job=rec['job'],contents=len(tests),rational_residuals_verified_by_exhaustive_polynomial=doc['monomials'],coefficient_controls=tests),indent=2)+'\n')
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);d=HERE/'classical_source';a.add_argument('--source',type=Path,default=HERE/'O34_C_even.json.gz');a.add_argument('--certificate',type=Path,default=d/'certificate.json');a.add_argument('--receipt',type=Path,default=d/'a8_5666124.json')
    for k in ('exe','work','out'):a.add_argument('--'+k,type=Path,required=True)
    assert __debug__;run(**vars(a.parse_args()))
