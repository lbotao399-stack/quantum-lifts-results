from pfaffian_relations import *
import hashlib,signal
start=time.monotonic();e4=Echelon()
for q in pickle.loads((HERE/'exact_q_relations_sp4_L11.pkl').read_bytes()):
 for part in profile_components(q).values():e4.add(part)
assert len(e4.rows)==151
oldrows,oldrecipes=pickle.loads((HERE/'target_snapshot_181.pkl').read_bytes());selected=[]
for i,(row,recipe) in enumerate(zip(oldrows,oldrecipes)):
 if not isinstance(recipe[0],str) and e4.add(row):selected.append(i)
assert len(e4.rows)==181 and len(selected)==30
(HERE/'transport_target_selected.json').write_text(json.dumps({'source_rank':4,'target_rank':5,'selected_source_recipe_indices':selected,'role':'Each graph is re-expanded as a rank-five Pfaffian identity. No cohomology transport is assumed.'},indent=2))
path=HERE/'exact_relations_sp5_L14.pkl';rows,recipes=pickle.loads(path.read_bytes());e5=Echelon()
for row in rows:assert e5.add(row)
finished_path=HERE/'transport_target_progress.json'
finished=set(json.loads(finished_path.read_text()).get('processed',[]) if finished_path.exists() else [])
for i in selected:
 if i in finished:continue
 letters,remaining=oldrecipes[i];newletters=letters[:10]+((0,0,(0,)),(0,0,(0,)))+letters[10:]
 shift=lambda a:a+4 if a>=20 else a
 newremaining=tuple(sorted([tuple(sorted((shift(a),shift(b)))) for a,b in remaining]+[(21,23)]))
 row=singlet(identity(newletters,5,newremaining),3)
 if row and e5.add(row):rows.append(row);recipes.append((newletters,newremaining))
 finished.add(i)
 tmp=path.with_suffix('.part');tmp.write_bytes(pickle.dumps((rows,recipes),protocol=4));tmp.replace(path)
 finished_path.write_text(json.dumps({'processed':sorted(finished),'target_rows':len(rows),'seconds':time.monotonic()-start},indent=2))
 print('expanded recipe',i,'target rank',len(rows),'seconds',time.monotonic()-start,flush=True)
print('DONE',len(rows),'seconds',time.monotonic()-start,flush=True)
if __name__=='__main__':pass
