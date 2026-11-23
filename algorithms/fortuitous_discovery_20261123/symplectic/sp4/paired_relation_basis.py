from pfaffian_relations import *
from flint import fmpz_mat
import signal

def canonical_items(items,odd):
 odditems=[v for v in items if odd(v)]
 if len(set(odditems))<len(odditems):return 0,()
 return (-1)**sum(x>y for i,x in enumerate(odditems) for y in odditems[i+1:]),tuple(sorted(items))

def canonical_pairs(pairs):
 sign=1;out=[]
 for pair in pairs:
  sg,q=canonical_items(pair,is_ferm)
  if not sg:return 0,()
  sign*=sg;out.append(q)
 sg,q=canonical_items(out,lambda v:sum(map(is_ferm,v))%2)
 return sign*sg,q

@lru_cache(None)
def paired_content(letters):
 if not letters:return ((),)
 first=letters[0];out=set()
 for other in sorted(set(letters[1:])):
  j=letters.index(other,1)
  for rest in paired_content(letters[1:j]+letters[j+1:]):
   sg,pairs=canonical_pairs(((first,other),)+rest)
   if sg:out.add(pairs)
 return tuple(sorted(out))

def raising_pairs(pairs):
 out=Counter()
 for i,pair in enumerate(pairs):
  for j,(a,b,theta) in enumerate(pair):
   if not b:continue
   nextpair=pair[:j]+((a+1,b-1,theta),)+pair[j+1:]
   sg,v=canonical_pairs(pairs[:i]+(nextpair,)+pairs[i+1:])
   if sg:out[v]+=(a+1)*sg
 return {v:c for v,c in out.items() if c}

def run(rank):
 start=time.monotonic();w=(3,3,2*rank-1,0,0);L=2*rank+2
 basis=sorted({p for content in m.contents(w,L) for p in paired_content(tuple(sorted(content)))})
 images=[raising_pairs(p) for p in basis];target=sorted({v for row in images for v in row});ids={v:i for i,v in enumerate(target)}
 A=fmpz_mat(len(target),len(basis))
 for j,row in enumerate(images):
  for v,c in row.items():A[ids[v],j]=c
 K,dim=A.nullspace();print('pair space',len(basis),len(target),'HW',dim,'seconds',time.monotonic()-start,flush=True)
 remaining=tuple((4*i+1,4*i+3) for i in range(L//2));expanded={};rows=[];recipes=[];ech=Echelon()
 for j in range(dim):
  coefficients=[int(K[i,j]) for i in range(len(basis))];g=math.gcd(*coefficients);coefficients=[v//g for v in coefficients]
  row=Counter();recipe=[]
  for i,c in enumerate(coefficients):
   if not c:continue
   if i not in expanded:expanded[i]=identity(tuple(l for pair in basis[i] for l in pair),rank,remaining)
   for v,n in expanded[i].items():row[v]+=c*n
   recipe.append((basis[i],c))
  row=primitive({v:n for v,n in row.items() if n});assert not spin_action(row,True)
  if row and ech.add(row):rows.append(row);recipes.append(('paired',recipe))
  print('expanded',j+1,'rank',len(rows),'terms',len(row),'seconds',time.monotonic()-start,flush=True)
 assert len(rows)==dim
 path=HERE/f'paired_relations_sp{rank}_L{L}.pkl';path.write_bytes(pickle.dumps((rows,recipes),protocol=4))
 (HERE/f'paired_relations_sp{rank}_L{L}.json').write_text(json.dumps({'rank':rank,'weight':w,'length':L,'input_pair_space':len(basis),'raising_rank':len(basis)-dim,'relation_dimension':dim,'complete':True,'seconds':time.monotonic()-start},indent=2))
 return rows,recipes
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--rank',type=int,default=4);args=a.parse_args()
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('paired relation cap')));signal.alarm(180)
 run(args.rank)
