"""Allow simultaneous changes of the high-content correction. Every candidate
is an exact signed combination of original parents; its classical variations
cancel formally in 69 contents before any three-colour identity is used.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction
from math import gcd,lcm
from pathlib import Path
import scalar_filtration as layer
from scalar_filtration import s,f,p,tri

def integral(vec):
    den=lcm(*(c.denominator for c in vec.values()));nums={u:int(c*den) for u,c in vec.items()};g=gcd(*nums.values());return {u:c//g for u,c in nums.items()}

def null_columns(ws,other,mid):
    # Column elimination in the complement of the middle content. Operations
    # are carried through on source coordinates and on the middle image.
    freq=Counter(v for w in ws for v in other[w]);coord={v:i for i,v in enumerate(sorted(freq,key=lambda k:(freq[k],k)))};bas={};out={}
    for w in sorted(ws,key=lambda u:(bool(mid[u]),len(other[u]),u)):
        v={coord[k]:Fraction(c) for k,c in other[w].items()};src={w:Fraction(1)}
        while v:
            i=min(v);c=v[i]
            if i not in bas:
                bas[i]=({k:t/c for k,t in v.items()},{k:t/c for k,t in src.items()});break
            b,u=bas[i];p.axpy(v,b,-c);p.axpy(src,u,-c)
        else:
            src=integral(src);image={};side={}
            for u,c in src.items():p.axpy(image,mid[u],c);p.axpy(side,other[u],c)
            assert not side
            if image:out[w]=(src,image)
    return len(bas),out

def prepare(args):
    h=layer.load_vector(args.work/'top_vector.json.gz');targets={c for c in layer.contents(11,0,layer.CHARGES) if c.count(0)>=6}|set(tri.CONTENTS);assert len(targets)==70
    rhs=p.rhs_rows(args.rhs,targets);b={w:c for w,c in rhs.items() if p.content(w)==tri.MID};allmid=set(tri.CONTENTS)
    for w,c in h.items():p.axpy(b,{v:t for v,t in f.q(w).items() if p.content(v) in allmid},-c)
    assert {p.content(w) for w in b}=={tri.MID}
    selected=sorted({w for v in set(b)|set(rhs) for w in f.parents(v)});M={};O={};blocks=defaultdict(list)
    for w in selected:
        im={v:t for v,t in f.q(w).items() if p.content(v) in targets};
        if not im:continue
        M[w]={v:t for v,t in im.items() if p.content(v)==tri.MID};O[w]={v:t for v,t in im.items() if p.content(v)!=tri.MID};blocks[f.key(w)].append(w)
    ids={};images={};rank=0
    for ws in blocks.values():
        r,sol=null_columns(ws,O,M);rank+=r
        for k,(u,y) in sol.items():ids[k]=u;images[k]=y
    loc=args.work/'coupled';z=layer.recipe(loc,images,b,{tri.MID},'Middle image of signed combinations preserving all 69 other projected contents; source labels refer to kernel_columns.json.gz.');s.write(loc/'kernel_columns.json.gz',[[s.hexes(w),[[s.hexes(u),c] for u,c in sorted(ids[w].items())]] for w in sorted(ids)]);s.write(loc/'preparation.json',dict(contents=sorted(c.hex() for c in targets),selected=len(selected),nonzero_parents=len(M),formal_side_rank=rank,side_kernel=len(M)-rank,nonzero_middle_columns=len(images),target_trace_words=len(z['word']),target_rhs_words=len(b)));print('parents; side rank, kernel; nonzero middle columns, words',len(M),rank,len(M)-rank,len(images),len(z['word']),flush=True)

def combine(args):
    loc=args.work/'coupled';co=layer.load_vector(args.work/'top_vector.json.gz');data=dict((s.unpack_hex(w),{s.unpack_hex(u):int(v) for u,v in terms}) for w,terms in s.read(loc/'kernel_columns.json.gz'))
    sol=s.read(loc/'solution.json.gz');
    for k,c in sol['solution']:p.axpy(co,data[s.unpack_hex(k)],Fraction(c))
    layer.save_vector(args.work/'joint_vector.json.gz',co);s.write(args.work/'joint_receipt.json',dict(contents=70,used=len(co),kernel_coefficients=sum(1 for _ in sol['solution']),physical_middle_coordinates=sol['verified_polynomial_coordinates'],mode='preserve side components on the original parent-word expansions'));print('joint rational vector',len(co),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('mode',choices=('prepare','combine'));a.add_argument('--work',required=True,type=Path);a.add_argument('--rhs',default=s.HERE/'negative_Q1.jsonl.gz',type=Path);b=a.parse_args();assert __debug__; (prepare if b.mode=='prepare' else combine)(b)
