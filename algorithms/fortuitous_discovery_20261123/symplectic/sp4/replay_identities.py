"""Rebuild every exact finite-colour identity from its algebraic recipe.
A seed antisymmetrizes ten fundamental indices in dimension eight, so it is zero
as an exact tensor identity. The graph expansion tracks antisymmetric symplectic
edge orientation and the Grassmann permutation sign. Every subsequent operation
is an explicit linear substitution, Taylor-profile projection, or actual Q0.
"""
from path_relation_completion import *
import argparse,hashlib

def rebuild(recipe,length,index,all_rows):
 if not isinstance(recipe[0],str):
  letters,remaining=recipe
  return singlet(identity(letters,4,remaining),3)
 kind=recipe[0]
 if kind=='paired':
  out=Counter();remaining=tuple((4*i+1,4*i+3) for i in range(length//2))
  for pairs,c in recipe[1]:
   for w,n in identity(tuple(l for pair in pairs for l in pair),4,remaining).items():out[w]+=c*n
 elif kind=='paths':
  out=Counter()
  for paths,c in recipe[1]:
   for w,n in expand_paths(paths,4).items():out[w]+=c*n
 elif kind=='Q':
  _,previous,j=recipe;assert previous==length-1
  out=m.q(all_rows[previous][j])
 elif kind=='Qcomponent':
  _,previous,j,profile=recipe;assert previous==length-1
  out=profile_components(m.q(all_rows[previous][j])).get(profile,{})
 elif kind=='component':
  _,j,profile=recipe;assert j<index
  out=profile_components(all_rows[length][j]).get(profile,{})
 else:raise ValueError(('unknown identity recipe',kind))
 return primitive({w:c for w,c in out.items() if c})

def main(output):
 start=time.monotonic();all_rows={};checks=[]
 for length,name,expected in [(10,'exact_relations_sp4_L10.pkl',67),(11,'exact_relations_sp4_L11.pkl',188),(12,'target_snapshot_181.pkl',181)]:
  rows,recipes=pickle.loads((HERE/name).read_bytes());assert len(rows)==len(recipes)==expected;all_rows[length]=rows;e=Echelon()
  for i,(row,recipe) in enumerate(zip(rows,recipes)):
   actual=rebuild(recipe,length,i,all_rows)
   assert actual==row,('identity recipe mismatch',length,i)
   assert not spin_action(row,True),('not highest weight',length,i)
   assert e.add(row),('dependent identity',length,i)
   if (i+1)%20==0:print('replayed',length,i+1,'seconds',time.monotonic()-start,flush=True)
  checks.append({'length':length,'identities':expected,'exact_recipes_verified':expected,'highest_weight_verified':expected,'independent_modulo':1000003,'file_sha256':hashlib.sha256((HERE/name).read_bytes()).hexdigest()})
 from trace_descent import dimensions
 char={}
 for length in (10,11,12):
  finite=sum(dimensions((3,3,7,0,0),length,8,'Sp').values());stable=sum(dimensions((3,3,7,0,0),length,28,'Sp').values())
  char[length]={'stable':stable,'finite':finite,'identities':stable-finite}
 assert char[10]['identities']==67 and char[11]['identities']==188 and char[12]['identities']==183
 result={'passed':True,'checks':checks,'exact_character_receipts':char,'complete_source_identity_space':True,'target_subspace_is_genuine_but_partial':True,'seconds':time.monotonic()-start}
 output.write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=HERE/'identity_replay_result.json');args=a.parse_args();main(args.out)
