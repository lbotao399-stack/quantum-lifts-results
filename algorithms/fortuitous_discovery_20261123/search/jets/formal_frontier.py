exec(open(__file__.replace('formal_frontier.py','character_scan.py')).read().split('start=time.monotonic()')[0])
from trace_descent import highest_weights
import mixed_incoming as m
from collections import defaultdict
from flint import nmod_mat
import pickle,hashlib
start=time.monotonic();w=(3,3,7,0,0);length=11
signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('bounded formal computation')))
signal.alarm(240)
rows={};report={'weight':w,'length':length,'prime':m.P}
def block_key(word):
 return tuple(sorted(tuple(sum(l[j] if j<2 else int(j-1 in l[2]) for l in tr) for j in range(5)) for tr in word))
def rank(images):
 blocks=defaultdict(list)
 owners={}
 parent=list(range(len(images)))
 def find(i):
  while parent[i]!=i:
   parent[i]=parent[parent[i]];i=parent[i]
  return i
 for i,r in enumerate(images):
  for key in {block_key(v) for v in r}:
   if key in owners:parent[find(i)]=find(owners[key])
   else:owners[key]=i
 for i,r in enumerate(images):
  if r:blocks[find(i)].append(r)
 total=0;largest=0
 for k,rs in blocks.items():
  ws=sorted({v for r in rs for v in r});ids={v:i for i,v in enumerate(ws)}
  print('rank block',len(rs),len(ws),flush=True)
  mat=nmod_mat(len(rs),len(ws),m.P)
  for i,r in enumerate(rs):
   for v,c in r.items():mat[i,ids[v]]=c
  total+=mat.rank();largest=max(largest,len(rs))
 return {'rank':total,'blocks':len(blocks),'largest_block_rows':largest,'rows':len(images)}
try:
 for ell in [10,11]:
  cache=HERE/f'formal_337_L{ell}.pkl'
  if cache.exists():rows[ell]=pickle.loads(cache.read_bytes())
  else:
   rows[ell]=[r for v in highest_weights(w,ell).values() for r in v]
   cache.write_bytes(pickle.dumps(rows[ell],protocol=4))
  print('basis',ell,len(rows[ell]),'seconds',time.monotonic()-start,flush=True)
  qcache=HERE/f'q_337_L{ell}.pkl'
  if qcache.exists():images=pickle.loads(qcache.read_bytes())
  else:
   images=[m.q(r) for r in rows[ell]]
   assert all(not m.q(r) for r in images)
   qcache.write_bytes(pickle.dumps(images,protocol=4))
  report[ell]=rank(images)
  report['seconds']=time.monotonic()-start
  (HERE/'formal_frontier.json').write_text(json.dumps(report,indent=2))
  print('map',ell,report[ell],'seconds',report['seconds'],flush=True)
 if 10 in report and 11 in report:
  report['stable_H_upper']=len(rows[11])-report[10]['rank']-report[11]['rank']
  (HERE/'formal_frontier.json').write_text(json.dumps(report,indent=2))
except TimeoutError:
 print('TIME LIMIT',time.monotonic()-start,flush=True)
finally:signal.alarm(0)
