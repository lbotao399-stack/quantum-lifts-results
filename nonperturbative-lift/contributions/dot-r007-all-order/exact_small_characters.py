"""Weyl coefficient of graded boson/exterior Taylor letters, integer arithmetic."""
from functools import lru_cache
from itertools import permutations, product
from collections import Counter
from math import comb

roots=((1,0),(0,1),(1,1))
orb=((0,0),(0,0),*roots,*tuple((-i,-j) for i,j in roots))
rho=(1,0,-1)
# For each total power n, coefficient for Torus weight (a,b).
@lru_cache(None)
def colour(n,ferm):
    def bins(j,n,q):
        if j==len(orb):
            if not n: yield q
            return
        a,b=orb[j]
        for k in range(min(n,1 if ferm else n)+1):
            yield from bins(j+1,n-k,(q[0]+k*a,q[1]+k*b))
    return tuple(Counter(bins(0,n,(0,0))).items())

def dims(Theta,r,Ys):
    letters=sorted(((a,b,*t) for a,b in product(range(r+2),range(r+1)) for t in product((0,1),repeat=3) if a+b+sum(t)),key=lambda v:(sum(v),v),reverse=True)
    @lru_cache(None)
    def at(i,n,W,g,h):
        if n==0:return int(not(any(W) or g or h))
        if i==len(letters) or sum(W)<n or max(abs(g),abs(h))>n:return 0
        letter=letters[i];ferm=(sum(letter[2:])+1)%2
        m=min(n,8 if ferm else n,*(x//y for x,y in zip(W,letter) if y))
        total=0
        for k in range(m+1):
            charges=tuple(x-k*y for x,y in zip(W,letter))
            for (a,b),c in colour(k,ferm):
                total+=c*at(i+1,n-k,charges,g-a,h-b)
        return total
    result={}
    for Y in Ys:
        out=0
        for w in permutations(range(3)):
            s=(-1)**sum(w[i]>w[j] for i in range(3) for j in range(i+1,3))
            f=tuple(Theta[i]+rho[i]-rho[w[i]] for i in range(3))
            if min(f)<0:continue
            for a,b,c in product((0,1),repeat=3):
                z=(-1)**(a+b+c)
                # 1 - e^{-α}; need weight + sum α in the character.
                g,h=(a+c,b+c)
                out+=s*z*(at(0,Y,(r,r,*f),g,h)-at(0,Y,(r+1,r-1,*f),g,h))
        result[Y]=out
    return result
do=[((2,1,0),1,range(3,6),[6,12,5]),((2,2,1),1,range(4,7),[64,101,66]),((1,1,0),2,range(3,6),[9,28,20])]
for t,r,Ys,expected in do:
    d=dims(t,r,Ys)
    assert list(d.values())==expected,(t,r,d,expected)
    print(t,r,d)
