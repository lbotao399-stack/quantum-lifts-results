"""Replay every Sp(5) identity recipe and complete exact source character deficit."""
from path_relation_completion import *
import argparse,hashlib

def main(output,checkpoint,seconds):
 start=time.monotonic();all_rows={};qcache={};checks=[]
 saved=json.loads(checkpoint.read_text()) if checkpoint.exists() else {'verified':{}}
 def qrow(length,index):
  key=(length,index)
  if key not in qcache:qcache[key]=m.q(all_rows[length][index])
  return qcache[key]
 def rebuild(recipe,length,index):
  if not isinstance(recipe[0],str):
   letters,remaining=recipe
   return singlet(identity(letters,5,remaining),3)
  kind=recipe[0]
  if kind=='paired':
   out=Counter();remaining=tuple((4*i+1,4*i+3) for i in range(length//2))
   for pairs,c in recipe[1]:
    for w,n in identity(tuple(l for pair in pairs for l in pair),5,remaining).items():out[w]+=c*n
  elif kind=='paths':
   out=Counter()
   for paths,c in recipe[1]:
    for w,n in expand_paths(paths,5).items():out[w]+=c*n
  elif kind=='Q':
   _,previous,j=recipe;assert previous==length-1
   out=qrow(previous,j)
  elif kind=='Qcomponent':
   _,previous,j,profile=recipe;assert previous==length-1
   out=profile_components(qrow(previous,j)).get(profile,{})
  elif kind=='component':
   _,j,profile=recipe;assert j<index
   out=profile_components(all_rows[length][j]).get(profile,{})
  else:raise ValueError(('unknown recipe',kind))
  return primitive({w:c for w,c in out.items() if c})
 for length,name,expected in [(12,'exact_relations_sp5_L12.pkl',67),(13,'exact_relations_sp5_L13.pkl',188),(14,'target_sp5_snapshot_181.pkl',181)]:
  rows,recipes=pickle.loads((HERE/name).read_bytes());assert len(rows)==len(recipes)==expected
  all_rows[length]=rows;e=Echelon();digest=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
  old=saved['verified'].setdefault(str(length),{'file_sha256':digest,'rows':[]})
  assert old['file_sha256']==digest
  for i,(row,recipe) in enumerate(zip(rows,recipes)):
   if i not in old['rows']:
    assert rebuild(recipe,length,i)==row,('identity recipe mismatch',length,i)
    old['rows'].append(i)
    checkpoint.write_text(json.dumps(saved,indent=2)+'\n')
   assert not spin_action(row,True),('not highest weight',length,i)
   assert e.add(row),('dependent identity',length,i)
   if (i+1)%20==0:print('replayed',length,i+1,'seconds',time.monotonic()-start,flush=True)
   if time.monotonic()-start>seconds:
    print('CHECKPOINT time cap; rerun to continue.',flush=True);return
  checks.append({'length':length,'identities':expected,'exact_recipes_verified':expected,
   'highest_weight_verified':expected,'independent_modulo':1000003,'file_sha256':digest})
 from trace_descent import dimensions
 char={}
 for length in (12,13,14):
  finite=sum(dimensions((3,3,9,0,0),length,10,'Sp').values())
  stable=sum(dimensions((3,3,9,0,0),length,32,'Sp').values())
  char[length]={'stable':stable,'finite':finite,'identities':stable-finite}
 assert char[12]['identities']==67 and char[13]['identities']==188 and char[14]['identities']==183
 result={'passed':True,'checks':checks,'exact_character_receipts':char,
  'complete_source_identity_space':True,'target_subspace_is_genuine_but_partial':True,
  'seconds_this_invocation':time.monotonic()-start}
 output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=HERE/'identity_replay_result.json')
 a.add_argument('--checkpoint',type=Path,default=HERE/'identity_replay_checkpoint.json')
 a.add_argument('--seconds',type=int,default=300);args=a.parse_args()
 main(args.out,args.checkpoint,args.seconds)
