"""Find stress-multiplet products saturating the three Cartan-visible physical grades."""
import argparse
from collections import Counter
from functools import cache
from itertools import product
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import sys


def rows_by_family(weight,length,graviton):
    # The three family labels are the numbers of internal dotted slots,
    # external dotted translations, and external odd translations.
    inside=tuple((p,q,a,b,2-p-q-a-b) for p,q in product(range(2),repeat=2)
                 for a in range(3-p-q) for b in range(3-p-q-a))
    options={}
    for s in inside:
        for d,e,mask in product(range(weight[0]+1),range(weight[1]+1),range(8)):
            g=(*s,d,e,mask)
            w=graviton.weight(g)
            if any(x>y for x,y in zip(w,weight)): continue
            f=(sum(s[:2]),d+e,mask.bit_count())
            options.setdefault(f,[]).append((g,w))
    families=sorted(options)
    max_d=sum(weight[:2])
    excess=sum(weight)-length
    # Enumerate only compatible degree choices, before enumerating components.
    @cache
    def patterns(start,left,dotted,outer):
        if not left: return [()] if dotted==outer==0 else []
        found=[]
        for j in range(start,len(families)):
            h,d,e=families[j]
            if h+d<=dotted and d+e<=outer:
                found.extend((families[j],)+tail for tail in patterns(j,left-1,dotted-h-d,outer-d-e))
        return found
    def words(pattern):
        @cache
        def build(i,left,first):
            if i==len(pattern): return [()] if not any(left) else []
            out=[]
            f=pattern[i]
            for k,(g,w) in enumerate(options[f][first:],first):
                if any(v>n for v,n in zip(w,left)): continue
                rem=tuple(n-v for v,n in zip(w,left))
                lo=(k+graviton.parity(g)) if i+1<len(pattern) and pattern[i+1]==f else 0
                out.extend((g,)+tail for tail in build(i+1,rem,lo))
            return out
        result=[]
        for raw in build(0,weight,0):
            sign,word=graviton.canonical(raw,graviton.parity)
            if sign: result.append(word)
        build.cache_clear()
        return sorted(result)
    # Small families first. This ordering only chooses witnesses.
    available=patterns(0,length//2,max_d,excess)
    for f in sorted(available,key=lambda f:(math.prod(len(options[i]) for i in f),f)):
        yield f,words(f)


def build(repo,t,first_length,stop):
    sys.path.insert(0,str(repo.resolve()/"algorithms/dual_singlets_20261118"))
    import sp3_omega_products as g
    from mixed33_exact import exact,apply
    from flint import fmpz_mat
    weight=(t,t,6-t,6-t,6-t)
    length=first_length-t
    rows=[]; recipes=[]; pivots={}
    for family,words in rows_by_family(weight,length,g):
        if not words: continue
        images=[g.raising(w) for w in words]
        for co in exact(images)[2]:
            assert not apply(co,images)
            v=Counter()
            for i,n in co.items(): v.update({w:n*c for w,c in g.cartan(words[i],1).items()})
            residue={key:c%1000003 for key,c in v.items() if c%1000003}
            while residue:
                key=min(residue)
                if key not in pivots:
                    inv=pow(residue[key],-1,1000003)
                    pivots[key]={k:c*inv%1000003 for k,c in residue.items()}
                    rows.append(dict(v))
                    recipes.append([[words[i],int(n)] for i,n in sorted(co.items())])
                    break
                v0=residue[key]
                for k,n in pivots[key].items():
                    new=(residue.get(k,0)-v0*n)%1000003
                    if new:
                        residue[k]=new
                    else:
                        residue.pop(k,None)
            if len(rows)>=stop: break
        if len(rows)>=stop: break
    columns=sorted(pivots)
    matrix=[[v.get(k,0) for k in columns] for v in rows]
    det=int(fmpz_mat(matrix).det()) if rows else 1
    assert det%1000003
    return dict(t=t,length=length,weight=weight,rank=len(rows),minor=str(det),columns=columns,recipes=recipes)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",required=True,type=Path)
    p.add_argument("--lease", action="store_true")
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    p.add_argument("--first-length",type=int,choices=(8,9),default=9)
    p.add_argument("grades",type=int,nargs="+",choices=(1,3,4,5))
    a=p.parse_args()
    assert __debug__
    bus=Path("~/.agent-bus/bin/agent-bus").expanduser()
    if a.lease:
        subprocess.run([bus,"lease",str(os.getpid()),"--agent","codex","--purpose",
                        "Build SU(2) level-36 protected product minors","--hours","2"],check=True)
    try:
        for t in a.grades:
            ranks={(9,1):1,(9,3):58,(9,5):6,(8,4):3}
            data=build(a.repo,t,a.first_length,ranks[a.first_length,t])
            (a.out/f"products-{t}.json.gz").write_bytes(gzip.compress(json.dumps(data,separators=(",",":")).encode(),mtime=0))
            print(t,data["rank"],data["minor"],flush=True)
    finally:
        if a.lease: subprocess.run([bus,"release",str(os.getpid())],check=True)
