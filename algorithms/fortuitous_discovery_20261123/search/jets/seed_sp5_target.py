from pfaffian_relations import *
import hashlib
start=time.monotonic();budget=150;rank=5;length=14
path=HERE/'exact_relations_sp5_L14.pkl';source_path=HERE/'exact_relations_sp5_L13.pkl';checkpoint=HERE/'sp5_q13_seed_progress.json'
qdir=HERE/'sp5_q13_identity_rows';qdir.mkdir(exist_ok=True)
if path.exists():rows,recipes=pickle.loads(path.read_bytes())
else:rows,recipes=[],[]
e=Echelon()
for row in rows:assert e.add(row)
source,_=pickle.loads(source_path.read_bytes());done=set(json.loads(checkpoint.read_text()).get('seeded_source_indices',[]) if checkpoint.exists() else [])
def save():
 tmp=path.with_suffix('.part');tmp.write_bytes(pickle.dumps((rows,recipes),protocol=4));tmp.replace(path)
 checkpoint.write_text(json.dumps({'rank':5,'source_length':13,'target_length':14,'source_rows_observed':len(source),'seeded_source_indices':sorted(done),'target_relation_rows':len(rows),'seconds':time.monotonic()-start},indent=2))
for i,row in enumerate(source):
 if i in done:continue
 if time.monotonic()-start>budget:break
 h=hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest();f=qdir/f'row_{i}.pkl'
 if f.exists():old,image=pickle.loads(f.read_bytes());assert old==h
 else:
  image=m.q(row);f.write_bytes(pickle.dumps((h,image),protocol=4))
 for profile,part in profile_components(image).items():
  if part and e.add(part):rows.append(primitive(part));recipes.append(('Qcomponent',13,i,profile))
 done.add(i)
 if i%5==0 or i==len(source)-1:save();print('seeded',i+1,'target',len(rows),'seconds',time.monotonic()-start,flush=True)
save();print('checkpoint target',len(rows),'seeded',len(done),'seconds',time.monotonic()-start,flush=True)
