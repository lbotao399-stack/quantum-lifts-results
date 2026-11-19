"""Standalone Sp(3) Cartan graviton-weight replay; standard library only.

Default: rerun one exact weight span and check all recorded Weyl arithmetic.
Pass --full to recompute all 12 exact polynomial spans; takes longer.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement as cwr, combinations, product
from time import monotonic
from pathlib import Path
import argparse
import json

p=1000003
target=(7,2,2,1,1)
dp={(0,0,0):1,(1,-1,0):-1,(0,1,-1):-1,(2,-1,-1):1,(1,1,-2):1,(2,0,-2):-1}
phi=[(0,0,(i,)) for i in range(3)]
def par(w):
 return (1+len(w[2]))%2


def parts(n,m=2):
 if n==0:yield ()
 for i in range(m,n+1):
  for s in parts(n-i,i):yield (i,)+s

def keyof(W):
 a=[x for x in W if par(x[0])]
 if len(set(a))!=len(a):return 0,()
 return (-1)**sum(x>y for i,x in enumerate(a) for y in a[i+1:]),tuple(sorted(W))

@lru_cache(maxsize=25000)
def cartan_trace(t):
 if len(t)%2:return ()
 R=defaultdict(int)
 for i in range(3):
  s,k=keyof(tuple((x,i) for x in t))
  if s:R[k]+=2*s
 return tuple(R.items())

def red(T,U,c=1):
 for w,a in U.items():
  z=(T.get(w,0)+c*a)%p
  if z:T[w]=z
  else:T.pop(w,None)

def cf(w,i,agent):
 if agent=='z':
  a,b,S=w
  return (a+int(i==0),b+int(i==1),S),1
 if agent=='t':
  a,b,S=w
  if i in S:return w,0
  h=(-1)**(len(S)+sum(t<i for t in S))
  return (a,b,tuple(sorted(S+(i,)))),h
 if agent=='f':
  u,v=i;a,b,S=w
  if v not in S or u in S:return w,0
  T=tuple(u if x==v else x for x in S)
  return (a,b,tuple(sorted(T))),(-1)**sum(x>y for j,x in enumerate(T) for y in T[j+1:])
 if agent=='J':
  a,b,S=w;k=b if i else a
  return ((a+1,b-1,S) if i else (a-1,b+1,S)),k
 raise ValueError(agent)

def act(P,i,k):
 R={}
 for w,c in P.items():
  x=0
  for j,(v,h) in enumerate(w):
   L,s=cf(v,i,k)
   if s:
    if k=='t':s*=(-1)**x
    f,W=keyof(w[:j]+((L,h),)+w[j+1:])
    if f:red(R,{W:f*s*c})
   x+=par(v)
 return R

def norm(P):
 if not P:return P
 a=pow(P[min(P)],-1,p)
 return {w:c*a%p for w,c in P.items()}

def cmult(A,B):
 R={}
 for w,c in A.items():
  for t,d in B.items():
   s,U=keyof(w+t)
   if s:red(R,{U:s*c*d})
 return R

def recipe(u,e,f,D,i,j):
 B=dict(cartan_trace(tuple(phi[k] for k in u)+((1,0,()),)*e+((0,1,()),)*f))
 for d in D:B=act(B,d,'t')
 for _ in range(i):B=act(B,0,'z')
 for _ in range(j):B=act(B,1,'z')
 return norm({w:c%p for w,c in B.items() if c%p})

def generic(T,n):
 R=[];W=set();h=sum(T)-10
 for e,f in product(range(2),repeat=2):
  if e>T[3] or f>T[4]:continue
  for u in cwr(range(3),n-e-f):
   q=[u.count(k) for k in range(3)]
   if any(q[z]>T[z] for z in range(3)):continue
   for i,j in product(range(T[3]-e+1),range(T[4]-f+1)):
    for t in range(min(3,h-i-j)+1):
     for D in combinations(range(3),t):
      Q=tuple(q[z]+int(z in D) for z in range(3))+(e+i,f+j)
      if any(a>b for a,b in zip(Q,T)):continue
      P=recipe(u,e,f,D,i,j)
      if not P:continue
      A=tuple(sorted(P.items()))
      if A in W:continue
      W.add(A);R.append((Q,P))
 D0=defaultdict(dict);B=[]
 for Q,V in R:
  A={tuple(w0 for w0,c0 in w):v for w,v in V.items() if all(c0==0 for w0,c0 in w)}
  D=D0[Q]
  while A:
   k=min(A);c=A[k]
   if k not in D:
    m=pow(c,-1,p);D[k]={w:x*m%p for w,x in A.items()};B.append((Q,V));break
   red(A,D[k],-c)
 return B

@lru_cache(maxsize=300000)
def canonical_color(w):
 G=[[] for i in range(3)]
 for x,i in w:G[i].append(x)
 return G[0]<=G[1]<=G[2]

def rank_weight(T):
 A={n:generic(T,n) for n in (2,4,6,8,10)}
 BY={n:defaultdict(list) for n in A}
 for n,G in A.items():
  for j,(q,P) in enumerate(G):BY[n][q].append(j)
 E={};loss=0;n0=0;tm=monotonic()
 def row(P):
  P={w:c for w,c in P.items() if canonical_color(w)}
  while P:
   w=min(P);c=P[w]
   if w not in E:
    a=pow(c,-1,p);E[w]={k:v*a%p for k,v in P.items()};return
   red(P,E[w],-c)
 for sh in sorted(parts(10),key=lambda x:(len(x),x)):
  if any(n not in A for n in sh):continue
  def rec(i,q,P,k0):
   nonlocal loss,n0
   if i==len(sh):
    assert q==T;row(P);n0+=1;return
   s=sh[i]
   I=BY[s].get(tuple(a-b for a,b in zip(T,q)),[]) if i+1==len(sh) else range(len(A[s]))
   for k in I:
    if i and sh[i-1]==s and k<k0:continue
    A0,B=A[s][k];Z=tuple(a+b for a,b in zip(A0,q))
    if any(a>b for a,b in zip(Z,T)) or sum(T)-sum(Z)<sum(sh[i+1:]):continue
    U=cmult(P,B)
    if U:rec(i+1,Z,U,k)
    else:loss+=1
  rec(0,(0,0,0,0,0),{():1},-1)
  print(T,sh,len(E),n0,loss,round(monotonic()-tm,2),flush=True)
 return len(E)


def main():
 cmd=argparse.ArgumentParser(description=__doc__)
 cmd.add_argument('--full',action='store_true')
 opt=cmd.parse_args()
 receipt=json.loads(Path(__file__).with_name('sp3_weight_ranks.json').read_text())
 C={tuple(v):int(x) for v,x in receipt['ranks']}
 assert receipt['prime']==p and receipt['completed']
 expected={(*tuple(a+b for a,b in zip(u,target[:3])),*v) for u in dp for v in ((1,1),(2,0))}
 assert set(C)==expected
 W=sorted(C) if opt.full else [(9,2,0,2,0)]
 for w in W:
  assert rank_weight(w)==C[w],w
 def contribution(u):
  z=tuple(a+b for a,b in zip(u,target[:3]))
  return C[(*z,1,1)]-C[(*z,2,0)]
 r=sum(s*contribution(u)for u,s in dp.items())
 assert r==receipt['exact_HW_mod_rank']==454
 print('Verified Weyl rank',r,'and',len(W),'full Cartan spans')

if __name__=='__main__':main()
