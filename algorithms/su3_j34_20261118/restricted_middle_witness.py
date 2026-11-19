"""An integer test on the middle component of -Q1 O0 - Q0(top).

A numerical letter label is the index in the archived sorted upstream alphabet.
Use the coefficient of the ordered product of five independent Grassmann generators.
The first scalar is diag(1,1,-2), where the centralizer contains a full 2x2 block.
This is a necessary condition on corrections to the fixed top vector, not a
functional on the unrestricted classical quotient. All arithmetic is exact.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import hashlib
import replay_formal_blocks as s
from fortuity.full.letters import Q_letter
from fortuity.full.words import Q_multi

MID=bytes.fromhex('0000000000010505050710')
# Each triple uses a row, a column and a Grassmann mask; the value is an integer.
def diagonal(x,y,mask):return {(i,i,mask):v for i,v in enumerate((x,y,-x-y)) if v}
M={0:diagonal(1,1,0),1:diagonal(1,1,1),5:{**diagonal(1,0,2),(0,1,4):1,(1,0,8):1},7:diagonal(1,0,16),16:diagonal(1,0,0)}

def wedge(a,b):
    if a&b:return 0
    return (-1)**sum((b&((1<<j)-1)).bit_count() for j in range(5) if a>>j&1)

def mul(a,b):
    r=Counter()
    for (i,j,u),x in a.items():
        for (k,h,v),y in b.items():
            if j==k:r[i,h,u|v]+=x*y*wedge(u,v)
    return {k:v for k,v in r.items() if v}

@lru_cache(None)
def trace(tr):
    a={(i,i,0):1 for i in range(3)}
    for k in tr:a=mul(a,M[k])
    r=Counter()
    for (i,j,t),v in a.items():
        if i==j:r[t]+=v
    return {k:v for k,v in r.items() if v}

def value(w):
    u={0:1}
    for tr in w:
        r=Counter()
        for i,x in u.items():
            for j,y in trace(tr).items():r[i|j]+=x*y*wedge(i,j)
        u=r
    return u.get(31,0)

def run(rhs,vector):
    hd=hashlib.sha256();num=0;seen=0;b=Fraction();d=Fraction();words=s.read(vector);table={};charges=Counter(MID)
    for w,c in s.stream_rows(rhs,hd):
        num+=1
        if bytes(sorted(b''.join(w)))==MID:b+=Fraction(c)*value(w);seen+=1
    assert hd.hexdigest()==s.INPUT and num==517468
    for w,c in words:
        for u,n in Q_multi(s.expanded(s.unpack_hex(w))).items():
            v=s.packed(u)
            if bytes(sorted(b''.join(v)))==MID:d+=Fraction(c)*n*value(v)
    for i,letter in enumerate(s.ALPH):
        for (a,h) in sorted({tuple(sorted((s.LETTER[u],s.LETTER[v]))) for n,u,v in Q_letter(letter)}):
            if (Counter((a,h))-charges):continue
            rem=charges-Counter((a,h));cc=bytes(sorted(list(rem.elements())+[i]));ex=Counter();pairs=[]
            for n,u,v in Q_letter(letter):
                j,k=s.LETTER[u],s.LETTER[v]
                if sorted((j,k))!=[a,h]:continue
                rr=mul(M[j],M[k]);pairs.append([n,f'{j:02x}',f'{k:02x}'])
                for coord,coef in rr.items():ex[coord]+=n*coef
            assert cc not in table;table[cc]=dict(letter=f'{i:02x}',splits=pairs,matrix=[[r+1,c+1,m,v] for (r,c,m),v in sorted(ex.items()) if v])
    assert len(table)==6 and {r['letter'] for r in table.values() if r['matrix']}=={'12','0d'} and b!=d
    return dict(input_sha256=s.INPUT,vector_sha256=hashlib.sha256(__import__('json').dumps(words,separators=(',',':')).encode()).hexdigest(),input_words=num,middle_content=MID.hex(),middle_rhs_words=seen,top_vector_words=len(words),letters={f'{k:02x}':[[i+1,j+1,t,v] for (i,j,t),v in sorted(a.items())] for k,a in M.items()},mask_order=[1,2,4,8,16],entry_order='one-based row, one-based column, Grassmann mask, signed coefficient',rhs=str(b),image_of_top=str(d),rhs_minus_image=str(b-d),parent_contents={k.hex():r for k,r in sorted(table.items())},scope='Every correction of this particular top vector on the middle content must use at least one predecessor from one of the two surviving contents. Other contents of the charge sector are not constrained by this test.')
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--vector',type=Path,default=s.HERE/'scalar_layers/top_vector.json.gz');a.add_argument('--out',type=Path,required=True);v=a.parse_args();r=run(v.rhs,v.vector);s.write(v.out,r);print(r['rhs'],r['image_of_top'],r['rhs_minus_image'],[(k,u['matrix']) for k,u in r['parent_contents'].items() if u['matrix']],flush=True)
