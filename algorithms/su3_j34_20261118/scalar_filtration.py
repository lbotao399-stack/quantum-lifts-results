"""Solve the first correction on successive levels of the number of copies of letter 00.

At length 11 classical trace words in charge (1,1,7,4,3) have 2..7 copies.
A parent term with k copies contributes to target terms with k or k+1 copies.
Filter on the number of copies before imposing finite-colour trace identities.
The top calculation enforces even those target contents absent from -Q1 O0.
"""
import argparse, gzip, hashlib, json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import lcm
from pathlib import Path
import three_content_physical as tri
from three_content_physical import s,f,p,ft

CHARGES=(1,1,7,4,3)
@lru_cache(maxsize=None)
def contents(n,minimum,charges):
    if not n:return () if any(charges) else (b'',)
    found=[]
    for i in range(minimum,len(f.ALPH)):
        rem=tuple(c-d for c,d in zip(charges,f.DEGREES[i]))
        if min(rem)<0 or max(rem)>n-1:continue
        found.extend(bytes([i])+v for v in contents(n-1,i,rem))
    return tuple(found)

def save_vector(path,d):s.write(path,[[s.hexes(w),str(c)] for w,c in sorted(d.items()) if c])
def load_vector(path):return {s.unpack_hex(t):Fraction(c) for t,c in s.read(path)}
def denominator(d):return lcm(*(c.denominator for c in d.values()))
def merge(x,images,co):
    for w,c in co.items():p.axpy(x,images[w],c)

def recipe(work,images,rhs,cc,scope):
    work.mkdir(parents=True,exist_ok=True);ws=sorted(set(rhs).union(*(m.keys() for m in images.values())));ix={w:i for i,w in enumerate(ws)};src=sorted(images);den=denominator(rhs)
    def encode(row):
        assert all(int(v)==v for v in row.values());return sorted([ix[w],int(v)] for w,v in row.items())
    z=dict(rhs_sha256=s.INPUT,contents=sorted(c.hex() for c in cc),source=[s.hexes(w) for w in src],word=[s.hexes(w) for w in ws],den=den,image=[encode(images[w]) for w in src],rhs=encode({w:c*den for w,c in rhs.items()}),odd=bytes(i for i,v in enumerate(ft.ODD) if v).hex(),scope=scope)
    s.write(work/'recipe.json.gz',z);(work/'words.txt').write_text(str(len(ws))+'\n'+''.join(f"{len(w)} {' '.join(s.hexes(w))}\n" for w in ws));return z

def search(work,images,rhs,cc,exe,jobs,scope):
    from physical_image_solve import solve,verify
    z=recipe(work,images,rhs,cc,scope);tri.columns(work,z,exe,jobs);A,keys=tri.load_columns(work,z);solve(work,z,A,keys)
    if not (work/'solution.json.gz').exists():return None
    verify(work,z,A,keys);return {s.unpack_hex(w):Fraction(c) for w,c in s.read(work/'solution.json.gz')['solution']}

def system(rhs,keep,parent):
    parents=sorted({w for v in rhs for w in f.parents(v) if parent(w)});images={}
    for w in parents:
        m={v:c for v,c in f.q(w).items() if p.content(v) in keep};
        if m:images[w]=m
    return images

def connected(images,rhs):
    inv=defaultdict(set);adj={}
    for w,m in images.items():
        adj[w]={p.content(v) for v in m}
        for c in adj[w]:inv[c].add(w)
    left=set(p.content(v) for v in rhs)
    while left:
        cc={min(left)};src=set()
        while True:
            src=set().union(*(inv[c] for c in cc));more=cc|set().union(*(adj[w] for w in src));
            if more==cc:break
            cc=more
        left-=cc;yield cc,{w:images[w] for w in sorted(src)},{w:c for w,c in rhs.items() if p.content(w) in cc}

def formal(images,rhs):
    ws=sorted(images);cert,coords,rem=p.rref([images[w] for w in ws],rhs);sol={ws[j]:v for j,v in coords.items()};diff=dict(rhs);merge(diff,images,{w:-v for w,v in sol.items()});assert diff=={w:Fraction(v) for w,v in rem.items()};return cert,sol

def top(args):
    allowed=contents(11,0,CHARGES);assert len(set(allowed))==len(allowed)==500
    cc={c for c in allowed if c.count(0)>=6};assert len(cc)==67
    b=p.rhs_rows(args.rhs,cc);images=system(b,cc,lambda w:p.content(w).count(0)==5);ans={};report=[]
    for t,(part,M,rhs) in enumerate(connected(images,b)):
        cert,individual=formal(M,rhs);entry=dict(contents=sorted(c.hex() for c in part),source_words=len(M),rhs_words=len(rhs),formal=cert);print('component',t,len(part),len(M),cert,flush=True)
        if cert['residual_terms']:
            individual=search(args.work/('component-%02d'%t),M,rhs,part,args.exe,args.jobs,'Projection onto this connected set of maximum-multiplicity contents, using parents with five copies of 00.')
            if individual is None:raise RuntimeError('Selected physical parent span inconclusive')
            d=s.read(args.work/('component-%02d'%t)/'solution.json.gz');entry.update(monomials=d['verified_polynomial_coordinates'],minor=d['selected_minor_determinant'],minor_order=len(d['pivot_word_indices']))
        entry['solution_terms']=len(individual);assert not set(ans).intersection(individual);ans.update(individual);report.append(entry)
    save_vector(args.work/'top_vector.json.gz',ans);s.write(args.work/'top_receipt.json',dict(charges=CHARGES,length=11,all_allowed=500,allowed_top=len(cc),rhs_sha256=s.INPUT,top_rhs_words=len(b),top_rhs_contents=len({p.content(w) for w in b}),distinct_parents=len(images),used=len(ans),components=report,scope='Every length-eleven content with at least six copies of 00, jointly.'))

def middle(args):
    allowed=contents(11,0,CHARGES);topcc={c for c in allowed if c.count(0)>=6};keep=set(tri.CONTENTS);raw=p.rhs_rows(args.rhs,keep);rhs=dict(raw);high=load_vector(args.work/'top_vector.json.gz')
    for w,v in high.items():p.axpy(rhs,{t:c for t,c in f.q(w).items() if p.content(t) in keep},-v)
    assert {p.content(t) for t in rhs}=={tri.MID};images=system(rhs,keep,lambda w:p.content(w).count(0)==4)
    # These sources have zero image in all 67 higher contents by letter count. Keep only those with zero image in the two side contents as well.
    images={w:m for w,m in images.items() if {p.content(v) for v in m}=={tri.MID}};r=args.work/'middle';co=search(r,images,rhs,keep,args.exe,args.jobs,'Residual after the common top correction; selected parents cannot enter the higher 67 contents or either side content.')
    if co is None:return
    assert not (set(high)&set(co));high.update(co);save_vector(args.work/'joint_vector.json.gz',high);s.write(args.work/'joint_receipt.json',dict(allowed=len(keep|topcc),top=len(topcc),additional=len(keep-topcc),additional_rhs=len(raw),adjusted_rhs=len(rhs),additional_sources=len(images),additional_used=len(co),joint_used=len(high)))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=('top','middle'));parser.add_argument('--work',type=Path,required=True);parser.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');parser.add_argument('--exe',type=Path,required=True);parser.add_argument('--jobs',type=int,default=2);args=parser.parse_args();assert __debug__;args.work.mkdir(parents=True,exist_ok=True);(top if args.mode=='top' else middle)(args)
