from relation_quotient import *
start=time.monotonic();source,_=pickle.loads((HERE/'exact_relations_sp5_L12.pkl').read_bytes());folder=HERE/'sp5_q12_identity_rows';folder.mkdir(exist_ok=True);images=[]
for i,row in enumerate(source):
 path=folder/f'row_{i}.pkl';h=hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest()
 if path.exists():old,q=pickle.loads(path.read_bytes());assert old==h
 else:
  q=m.q(row);path.write_bytes(pickle.dumps((h,q),protocol=4))
 images.append(q)
 if i%5==0:print('Q identity',i+1,'seconds',time.monotonic()-start,flush=True)
 if time.monotonic()-start>150:break
if len(images)==67:
 (HERE/'exact_q_relations_sp5_L12.pkl').write_bytes(pickle.dumps(images,protocol=4))
 rank=row_rank(images,2147483647);assert rank==67
 (HERE/'sp5_preceding_relation_rank.json').write_text(json.dumps({'rank':rank,'source_dimension':67,'prime':2147483647,'exact_rational_rank':67,'seconds':time.monotonic()-start},indent=2))
 print('COMPLETE rank',rank,'seconds',time.monotonic()-start,flush=True)
