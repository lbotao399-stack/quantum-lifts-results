"""Exact signed-dihedral trace algebra and classical differential Q0."""
from fractions import Fraction
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import product
from math import comb
from sys import argv
def add(d,e,c):
    for k,v in e.items():
        s=d.get(k,0)+v*c
        if s:d[k]=s
        else:d.pop(k,None)

def canonical(seq, odd):
    odd_items=[x for x in seq if odd(x)]
    if len(set(odd_items))!=len(odd_items):
        return 0,()
    return (-1)**sum(x>y for i,x in enumerate(odd_items) for y in odd_items[i+1:]),tuple(sorted(seq))
W=(5,2,2,1,1)
def wt(a):
 m,p,q=a;return tuple((m>>i)&1 for i in range(3))+(p,q)
def odd(a):return (a[0].bit_count()+1)%2
def op(t):return sum(map(odd,t))%2
L=[(m,p,q) for m,p,q in product(range(8),range(2),range(2)) if m or p+q]
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def count(t):return tuple(sum(wt(a)[j] for a in t) for j in range(5))
@lru_cache(None)
def cyc(t):
 d=[];k=sum(map(odd,t));r=t[::-1]
 for w,s in [(t,1),(r,(-1)**(len(t)+k*(k-1)//2))]:
  for i in range(len(w)):
   a=sum(map(odd,w[:i]))
   d.append((w[i:]+w[:i],s*(-1)**(a*(k-a))))
 v=min(w for w,s in d);ss={s for w,s in d if w==v}
 return v,ss.pop() if len(ss)==1 else 0
@lru_cache(None)
def words(n):
 D=defaultdict(list)
 for t in product(L,repeat=n):
  w=count(t)
  if min(sub(W,w))<0:continue
  z,s=cyc(t)
  if s and z==t:D[w].append(t)
 return D
def shape(lens):
 F=[]
 def go(j,r,t):
  if j==len(lens):
   if not any(r):
    s,u=canonical(t,op)
    if s:F.append(u)
   return
  for v,K in words(lens[j]).items():
   k=sub(r,v)
   if min(k)<0:continue
   for z in K:
    if j and len(t[-1])==lens[j] and t[-1]>z:continue
    go(j+1,k,t+(z,))
 go(0,W,())
 return sorted(set(F))
def ra(a,k):
 m,p,q=a
 if k==0: g,s={2:(1,1),6:(5,-1)}.get(m,(m,0));return [((g,p,q),s)] if s else []
 if k==1: g,s={4:(2,1),5:(3,-1)}.get(m,(m,0));return [((g,p,q),s)] if s else []
 return [((m,p+1,q-1),q)] if q else []
@lru_cache(None)
def Q(a):
 m,p,q=a; out=Counter()
 ab=[(0,m,-1)]
 if m==0:ab=[(0,0,Fraction(-1,2))]
 if m in (3,5,6):
  i={6:0,5:1,3:2}[m];ab.append((1<<((i+1)%3),1<<((i+2)%3),1))
 if m==7:ab += [(h,7^h,-1) for h in (1,2,4)]
 # exact ghost c,c coefficient, which combines with the doubling of odd bracket
 for a,b,c in ab:
  c=Fraction(c)
  for i,j in product(range(p+1),range(q+1)):
   A=(a,i,j);B=(b,p-i,q-j)
   if not(any(A) and any(B)):continue
   r=c*comb(p,i)*comb(q,j)
   out[A,B]+=r;out[B,A]-=r*(-1)**(odd(A)*odd(B))
 return dict(out)
@lru_cache(None)
def deriv(t,k):
 z=Counter()
 for i,a in enumerate(t):
  f=Q(a).items() if k==3 else (((b,),s) for b,s in ra(a,k))
  for v,b in f:
   u,s=cyc(t[:i]+v+t[i+1:])
   if s:z[u]+=b*s*(-1)**((k==3)*sum(map(odd,t[:i])))
 return {k:v for k,v in z.items() if v}
def move(w,k):
 z=Counter()
 for i,t in enumerate(w):
  for u,c in deriv(t,k).items():
   s,v=canonical(w[:i]+(u,)+w[i+1:],op)
   if s:z[v]+=s*c*(-1)**((k==3)*sum(map(op,w[:i])))
 return {k:v for k,v in z.items() if v}
