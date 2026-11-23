"""Exact Sp finite-colour trace identities from an antisymmetrizer on 2n+2 indices.
The adjoint matrix X is encoded by the symmetric fundamental tensor J X.
Marked half-edges are antisymmetrized; unmarked half-edges are contracted with J^-1.
This vanishes identically in fundamental dimension 2n, with integral coefficients.
"""
from pathlib import Path
import os,sys,math,pickle,json,time,random
from collections import Counter
from functools import lru_cache
HERE=Path(__file__).resolve().parent;ROOT=Path(os.environ.get('N4_WORKSPACE_ROOT',str(HERE.parents[2])))
os.environ.setdefault('N4_TRACE_RUNTIME',str(ROOT/'quantum_lifts_frontier_source/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync'))
sys.path.insert(0,os.environ.get('N4_DUAL_SOURCE',str(ROOT/'quantum_lifts_frontier_work/root_audit/so11_pin/algorithms/dual_singlets_20261118')))
import mixed_incoming as m
from fortuity.full.letters import is_ferm

@lru_cache(None)
def pairings(points):
 if not points:return ((1,()),)
 a=points[0];rest=points[1:];out=[]
 for j,b in enumerate(rest):
  for sign,pairs in pairings(rest[:j]+rest[j+1:]):out.append(((-1)**j*sign,((a,b),)+pairs))
 return tuple(out)

def graph_traces(length,edges):
 adjacency={}
 for a,b in edges:
  assert a<b
  adjacency[a]=b;adjacency[b]=a
 assert len(adjacency)==2*length
 seen=set();traces=[];phase=1;order=[]
 for start in range(length):
  if start in seen:continue
  vertex=start;incoming=2*start;cycle=[]
  while vertex not in seen:
   seen.add(vertex);cycle.append(vertex);order.append(vertex)
   outgoing=incoming^1;nxt=adjacency[outgoing]
   phase*=1 if outgoing<nxt else -1
   incoming=nxt;vertex=nxt//2
  assert vertex==start
  traces.append(tuple(cycle))
 return phase,tuple(traces),tuple(order)

@lru_cache(maxsize=32)
def skeleton(length,rank,remaining):
 marked=tuple(2*i for i in range(2*rank+2))
 result=[]
 for sign,edges in pairings(marked):
  phase,traces,order=graph_traces(length,edges+remaining)
  result.append((sign*phase,traces,order))
 return tuple(result)

def identity(letters,rank,remaining):
 odd={i for i,l in enumerate(letters) if is_ferm(l)}
 out=Counter()
 for sign,traces,order in skeleton(len(letters),rank,remaining):
  fermions=[i for i in order if i in odd]
  sign*=(-1)**sum(a>b for i,a in enumerate(fermions) for b in fermions[i+1:])
  sg,v=m.so_canon(tuple(tuple(letters[i] for i in tr) for tr in traces))
  if sg:out[v]+=sign*sg
 return {v:c for v,c in out.items() if c}

def spin_action(row,raise_spin):
 out=Counter()
 for word,c in row.items():
  for t,tr in enumerate(word):
   for i,(a,b,theta) in enumerate(tr):
    fac=(a+1 if b else 0) if raise_spin else (b+1 if a else 0)
    if not fac:continue
    l=(a+1,b-1,theta) if raise_spin else (a-1,b+1,theta)
    sg,v=m.so_canon(word[:t]+(tr[:i]+(l,)+tr[i+1:],)+word[t+1:])
    if sg:out[v]+=c*fac*sg
 return {v:c for v,c in out.items() if c}

def primitive(row):
 if not row:return {}
 g=math.gcd(*row.values())
 if row[min(row)]<0:g=-g
 return {v:c//g for v,c in row.items()}

def singlet(row,spin=3):
 for j in range(1,spin+1):
  returned=spin_action(spin_action(row,True),False)
  ans=Counter({v:j*(j+1)*c for v,c in row.items()})
  for v,c in returned.items():ans[v]-=c
  row=primitive({v:c for v,c in ans.items() if c})
 assert not spin_action(row,True)
 check=Counter()
 for w,c in row.items():
  for v,d in m.raise_word(w,3).items():check[v]+=c*d
 assert not any(check.values())
 return row

class Echelon:
 def __init__(self,p=1000003):self.p=p;self.rows={};self.words={}
 def add(self,row):
  import numpy as np
  entries=[(self.words.setdefault(v,len(self.words)),c%self.p) for v,c in row.items() if c%self.p]
  r=np.zeros(len(self.words),dtype=np.int64)
  for i,c in entries:r[i]=c
  for pivot,(indices,values) in sorted(self.rows.items()):
   fac=int(r[pivot])
   if fac:r[indices]=(r[indices]-fac*values)%self.p
  indices=np.flatnonzero(r)
  if not len(indices):return False
  pivot=int(indices[0]);inv=pow(int(r[pivot]),-1,self.p)
  values=r[indices]*inv%self.p
  self.rows[pivot]=(indices,values)
  return True

def random_remaining(length,rank,rng):
 base=2*rank+2;extra=length-base
 templates={0:[([],[])],1:[([1],[])],2:[([2],[]),([1,1],[]),([],[2])],3:[([3],[]),([2,1],[]),([1,1,1],[]),([1],[2]),([],[3])]}
 chains,cycles=rng.choice(templates[extra]);edges=[];v=base
 for i in range(base//2):
  a,b=4*i+1,4*i+3;n=chains[i] if i<len(chains) else 0
  left=a
  for j in range(n):
   edges.append((left,2*v));left=2*v+1;v+=1
  edges.append((left,b))
 for n in cycles:
  vertices=list(range(v,v+n));v+=n
  for x,y in zip(vertices,vertices[1:]+vertices[:1]):edges.append((2*x+1,2*y))
 assert v==length
 return tuple(sorted(tuple(sorted(e)) for e in edges))

def profile_word(w):
 return tuple(sorted((a+b,len(theta)) for tr in w for a,b,theta in tr))

def profile_components(row):
 result={}
 for w,c in row.items():result.setdefault(profile_word(w),{})[w]=c
 return result

def profile_inventory(rank,w,length):
 from sp_dimensions import multiplicity_parts
 result=Counter();by_profile={}
 for shifted,sign in [(w,1),((w[0]+1,w[1]-1,*w[2:]),-1)]:
  for content in m.contents(shifted,length):
   profile=tuple(sorted((a+b,len(theta)) for a,b,theta in content))
   counts=Counter(content)
   bos=tuple(sorted((n for l,n in counts.items() if not is_ferm(l)),reverse=True))
   ferm=tuple(sorted((n for l,n in counts.items() if is_ferm(l)),reverse=True))
   result[profile]+=sign*(multiplicity_parts(bos,ferm,length+1)-multiplicity_parts(bos,ferm,rank))
   if sign==1:by_profile.setdefault(profile,[]).append(content)
 assert all(n>=0 for n in result.values())
 return {p:n for p,n in result.items() if n},by_profile

def generate(rank,w,length,expected,seconds=240):
 start=time.monotonic();rng=random.Random(1837321+rank*10000+length)
 contents=list(m.contents(w,length));ech=Echelon();rows=[];recipes=[]
 profile_dimensions,by_profile=profile_inventory(rank,w,length)
 assert sum(profile_dimensions.values())==expected
 path=HERE/f'exact_relations_sp{rank}_L{length}.pkl'
 report={'rank':rank,'weight':w,'length':length,'exact_character_relation_dimension':expected,'relations':0,'trials':0}
 if path.exists():
  rows,recipes=pickle.loads(path.read_bytes())
  for r in rows:
   assert not spin_action(r,True)
   assert ech.add(r)
 else:
  recipes_seen=set()
  for old in (HERE/f'relations_sp{rank}_L{length}.pkl',HERE/f'initial_unnormalized_relations_sp{rank}_L{length}.pkl'):
   if not old.exists():continue
   oldrows,oldrecipes=pickle.loads(old.read_bytes())
   for oldrow,(letters,remaining) in zip(oldrows,oldrecipes):
    key=(letters,remaining)
    if key in recipes_seen:continue
    recipes_seen.add(key)
    row=oldrow if not spin_action(oldrow,True) else singlet(identity(letters,rank,remaining),w[1])
    if row and ech.add(row):rows.append(row);recipes.append(key)
  path.write_bytes(pickle.dumps((rows,recipes),protocol=4))
  print('clean reconstruction',rank,length,len(rows),'seconds',time.monotonic()-start,flush=True)
 previous=HERE/f'exact_relations_sp{rank}_L{length-1}.pkl'
 if previous.exists():
  preceding,_=pickle.loads(previous.read_bytes())
  for i,source in enumerate(preceding):
   row=m.q(source)
   assert not m.q(row) and not spin_action(row,True)
   for profile,part in profile_components(row).items():
    if part and ech.add(part):rows.append(primitive(part));recipes.append(('Qcomponent',length-1,i,profile))
  path.write_bytes(pickle.dumps((rows,recipes),protocol=4))
  print('Q-seeded',rank,length,len(rows),'seconds',time.monotonic()-start,flush=True)
 for i,row in list(enumerate(rows)):
  for profile,part in profile_components(row).items():
   if part and ech.add(part):rows.append(primitive(part));recipes.append(('component',i,profile))
 profile_echelons={p:Echelon() for p in profile_dimensions}
 for row in rows:
  for profile,part in profile_components(row).items():
   assert profile in profile_dimensions
   profile_echelons[profile].add(part)
 assert sum(len(e.rows) for e in profile_echelons.values())==len(rows)
 assert all(len(profile_echelons[p].rows)<=n for p,n in profile_dimensions.items())
 path.write_bytes(pickle.dumps((rows,recipes),protocol=4))
 print('profile seed',rank,length,len(rows),'missing',[(p,n-len(profile_echelons[p].rows)) for p,n in profile_dimensions.items() if n>len(profile_echelons[p].rows)],flush=True)
 rng.seed(284239+rank*1000+length+len(rows))
 try:
  while len(rows)<expected and time.monotonic()-start<seconds:
   missing=[p for p,n in profile_dimensions.items() if len(profile_echelons[p].rows)<n]
   chosen=rng.choice(missing)
   letters=list(rng.choice(by_profile[chosen]));rng.shuffle(letters)
   remaining=random_remaining(length,rank,rng)
   raw=identity(tuple(letters),rank,remaining)
   row=singlet(raw,w[1])
   report['trials']+=1
   if row and profile_echelons[chosen].add(row):
    assert ech.add(row)
    rows.append(row);recipes.append((tuple(letters),remaining))
    report.update(relations=len(rows),seconds=time.monotonic()-start,total_terms=sum(map(len,rows)))
    temporary=path.with_suffix('.part')
    temporary.write_bytes(pickle.dumps((rows,recipes),protocol=4));temporary.replace(path)
    (HERE/f'exact_relations_sp{rank}_L{length}.json').write_text(json.dumps(report,indent=2))
    print('relation',rank,length,len(rows),'trial',report['trials'],'terms',len(row),'seconds',report['seconds'],flush=True)
 finally:
  report.update(relations=len(rows),seconds=time.monotonic()-start,total_terms=sum(map(len,rows)),complete=len(rows)==expected)
  (HERE/f'exact_relations_sp{rank}_L{length}.json').write_text(json.dumps(report,indent=2))
 return rows,recipes

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--rank',type=int,default=4);p.add_argument('--length',type=int,default=10);p.add_argument('--expected',type=int,default=67);p.add_argument('--seconds',type=int,default=240)
 a=p.parse_args();generate(a.rank,(3,3,2*a.rank-1,0,0),a.length,a.expected,a.seconds)
