"""Compare rank-three Cartan product minors with both full-letter cohomology bounds."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import cache
from itertools import product
from pathlib import Path
from time import monotonic
import gzip
import json

from graviton import cartan, factor_polynomial, parity

HERE = Path(__file__).resolve().parent
DATA = HERE
PRIME = 1009


def subtract(a,b):
    return tuple(x-y for x,y in zip(a,b))


def mod_insert(row, echelon, p):
    row = {k:v%p for k,v in row.items() if v%p}
    while row:
        pivot = min(row)
        if pivot not in echelon:
            k = row[pivot]
            inv = pow(k,-1,p)
            echelon[pivot] = {j:n*inv%p for j,n in row.items()}
            return pivot,k
        k = row[pivot]
        for j,n in echelon[pivot].items():
            r = (row.get(j,0)-k*n) % p
            if r:
                row[j] = r
            else:
                row.pop(j,None)
    return None


@cache
def factors(n,w):
    """Select independent single-colour components of the polarized trace."""
    a,b,*t = w
    choices = []
    frame = {}
    for h,k,mask in product(range(min(1,a)+1), range(min(1,b)+1), range(8)):
        m = tuple(c-int(bool(mask & (1<<i))) for i,c in enumerate(t))
        if min(m)<0 or h+k+sum(m)!=n:
            continue
        g = (h,k,*m,a-h,b-k,mask)
        if mod_insert(factor_polynomial(g),frame,PRIME):
            choices.append(g)
    return tuple(choices)


def read_census(group,level,length=8):
    files = gzip.decompress((DATA/f'{group}_full_j{level}.jsonl.gz').read_bytes()).splitlines()
    raw = [json.loads(line) for line in files]
    assert all(r['check'] for r in raw)
    rows = {(tuple(r['sec']),r['Y']):r for r in raw}
    assert len(rows)==len(raw)
    weight = Counter()
    for (sec,y),r in rows.items():
        if y!=length:
            continue
        spin = sec[:2]
        down = tuple(i-1 for i in spin)+sec[2:]
        h = r['m']-r['r_up']-rows.get((down,y+1),{}).get('r_up',0)
        if not h:
            continue
        d = level-2*y-sum(spin)
        a = (d+(spin[0]-spin[1])//2)//2
        t = tuple((i-d)//2 for i in sec[2:])
        assert 2*a == d+(spin[0]-spin[1])//2 and 0<=a<=d
        assert list(t)==sorted(t,reverse=True) and min(t)>=0 and h>0
        for u in range(t[1],t[0]+1):
            for v in range(t[2],t[1]+1):
                for k in range(v,u+1):
                    weight[a,d-a,k,u+v-k,sum(t)-u-v] += h
    return weight


def patterns(length,maximum):
    if length==0:
        yield ()
    for k in range(min(maximum,length),1,-2):
        for tail in patterns(length-k,k):
            yield (k,)+tail


def products(target,length):
    choices = defaultdict(list)
    for w in product(*(range(n+1) for n in target)):
        for n in range(2,length+1,2):
            if n<=sum(w)<=sum(target)-length+n:
                choices[n].extend((g,w) for g in factors(n,w))
    for items in choices.values():
        items.sort()
    def assign(pattern,rem,start,prefix):
        n,*tail = pattern
        if not tail:
            for g in factors(n,rem):
                if start is None or g>start or g==start and not parity(g):
                    yield tuple(prefix+[g])
            return
        for g,w in choices[n]:
            if start is not None and (g<start or g==start and parity(g)):
                continue
            delta = subtract(rem,w)
            if min(delta)<0 or sum(delta)<sum(tail):
                continue
            yield from assign(tail,delta,g if tail[0]==n else None,prefix+[g])
    for pattern in sorted(patterns(length,length),key=lambda t:(len(t),tuple(-k for k in t))):
        yield from assign(pattern,target,None,[])


def run(w,limit,do_rank,length=8,progress=True):
    counts = Counter()
    frame = {}
    ids = {}
    pivots = []
    cutoff = time0 = monotonic()
    for n,g in enumerate(products(w,length)):
        counts[len(g)] += 1
        if not do_rank:
            continue
        row = {ids.setdefault(k,len(ids)):v for k,v in cartan(g).items()}
        chosen = mod_insert(row,frame,PRIME)
        if chosen:
            pivots.append([n,*chosen])
            if len(pivots)==limit:
                break
        if progress and monotonic()-cutoff>30:
            print('tested',n,'independent',len(pivots),'target',limit,flush=True)
            cutoff=monotonic()
    det=1
    for _,_,v in pivots:
        det = det*v % PRIME
    return dict(weight=w,length=length,prime=PRIME,product_counts=dict(counts),limit=limit,cartan_rank=len(pivots),
                selected_rows=[k for k,_,_ in pivots],selected_columns=[k for _,k,_ in pivots],
                determinant=det)


if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--weight',type=int,nargs=5,required=True)
    p.add_argument('--counts-only',action='store_true')
    p.add_argument('--length',type=int,default=8)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    w=tuple(a.weight)
    level=3*sum(w[:2])+2*sum(w[2:])
    upper=read_census('usp6',level,a.length)[w]
    print('upper',upper,'weight',w,flush=True)
    r=run(w,upper,not a.counts_only,a.length)
    print(r['product_counts'],r['cartan_rank'],r['determinant'],flush=True)
    a.out.write_text(json.dumps(r,separators=(',',':'))+'\n')
