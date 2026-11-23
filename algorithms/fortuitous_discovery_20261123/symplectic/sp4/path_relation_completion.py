from paired_relation_basis import *
from itertools import combinations,permutations
import numpy as np

def path_odd(path):return sum(map(is_ferm,path))%2

def canonical_paths(paths):
 phase=1;out=[]
 for path in paths:
  reverse=path[::-1];odd=sum(map(is_ferm,path));sign=(-1)**(len(path)+odd*(odd-1)//2)
  if path==reverse and sign==-1:return 0,()
  if reverse<path:path=reverse;phase*=sign
  out.append(path)
 sg,v=canonical_items(out,path_odd)
 return phase*sg,v

@lru_cache(None)
def one_triple_content(letters):
 out=set()
 for inds in combinations(range(len(letters)),3):
  selected=tuple(letters[i] for i in inds);rest=tuple(l for i,l in enumerate(letters) if i not in inds)
  for triple in set(permutations(selected)):
   for pairs in paired_content(rest):
    sg,v=canonical_paths((triple,)+pairs)
    if sg:out.add(v)
 return tuple(sorted(out))

def raising_paths(paths):
 out=Counter()
 for i,path in enumerate(paths):
  for j,(a,b,theta) in enumerate(path):
   if not b:continue
   new=path[:j]+((a+1,b-1,theta),)+path[j+1:]
   sg,v=canonical_paths(paths[:i]+(new,)+paths[i+1:])
   if sg:out[v]+=(a+1)*sg
 return {v:c for v,c in out.items() if c}

def expand_paths(paths,rank):
 endpoints=[(i,j) for i,path in enumerate(paths) for j in (0,len(path)-1)]
 internal=[(i,j) for i,path in enumerate(paths) for j in range(1,len(path)-1)]
 tokens=endpoints+internal;positions={v:i for i,v in enumerate(tokens)}
 letters=tuple(paths[i][j] for i,j in tokens);marked=len(endpoints);edges=[]
 for i,path in enumerate(paths):
  left=2*positions[i,0]+1
  for j in range(1,len(path)-1):
   vertex=positions[i,j];edges.append((left,2*vertex));left=2*vertex+1
  edges.append((left,2*positions[i,len(path)-1]+1))
 remaining=tuple(sorted(tuple(sorted(e)) for e in edges))
 original=[positions[i,j] for i,path in enumerate(paths) for j in range(len(path)) if is_ferm(path[j])]
 phase=(-1)**sum(a>b for i,a in enumerate(original) for b in original[i+1:])
 return {v:phase*c for v,c in identity(letters,rank,remaining).items()}

def complete(rank=4):
 start=time.monotonic();w=(3,3,2*rank-1,0,0);L=2*rank+3
 cache=HERE/f'exact_relations_sp{rank}_L{L}.pkl';rows,recipes=pickle.loads(cache.read_bytes());e=Echelon()
 for row in rows:assert e.add(row)
 dims,contents=profile_inventory(rank,w,L);ep={p:Echelon() for p in dims}
 for row in rows:
  for p,part in profile_components(row).items():ep[p].add(part)
 missing=[p for p,n in dims.items() if len(ep[p].rows)<n]
 print('initial',len(rows),'missing',[(dims[p]-len(ep[p].rows)) for p in missing],'seconds',time.monotonic()-start,flush=True)
 for prof in missing:
  basis=sorted({v for content in contents[prof] for v in one_triple_content(tuple(sorted(content)))})
  images=[raising_paths(v) for v in basis];target=sorted({v for r in images for v in r});ids={v:i for i,v in enumerate(target)}
  A=fmpz_mat(len(target),len(basis))
  for j,row in enumerate(images):
   for v,c in row.items():A[ids[v],j]=c
  K,dim=A.nullspace();expanded={}
  print('path space',len(basis),len(target),'HW',dim,'seconds',time.monotonic()-start,flush=True)
  for j in range(dim):
   co=[int(K[i,j]) for i in range(len(basis))];g=math.gcd(*co);co=[c//g for c in co];r=Counter();recipe=[]
   for i,c in enumerate(co):
    if not c:continue
    if i not in expanded:expanded[i]=expand_paths(basis[i],rank)
    for v,n in expanded[i].items():r[v]+=c*n
    recipe.append((basis[i],c))
   r=primitive({v:n for v,n in r.items() if n});assert not spin_action(r,True)
   if r and e.add(r):
    rows.append(r);recipes.append(('paths',recipe));assert ep[prof].add(r)
    temp=cache.with_suffix('.part');temp.write_bytes(pickle.dumps((rows,recipes),protocol=4));temp.replace(cache)
    print('new',len(rows),'terms',len(r),'seconds',time.monotonic()-start,flush=True)
   if len(ep[prof].rows)==dims[prof]:break
  assert len(ep[prof].rows)==dims[prof]
 assert len(rows)==188
 (HERE/f'exact_relations_sp{rank}_L{L}.json').write_text(json.dumps({'rank':rank,'weight':w,'length':L,'relations':len(rows),'exact_character_relation_dimension':188,'complete':True,'total_terms':sum(map(len,rows)),'seconds':time.monotonic()-start},indent=2))
 print('COMPLETE',len(rows),'seconds',time.monotonic()-start,flush=True)
if __name__=='__main__':
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('path completion cap')));signal.alarm(180)
 complete()
