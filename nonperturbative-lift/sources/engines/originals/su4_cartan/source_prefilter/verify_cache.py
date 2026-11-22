"""Read back the saved source minor and three complete source-polynomial sketches."""
from pathlib import Path
import hashlib,json,pickle,time
from flint import nmod_mat,ctx
ROOT=Path(__file__).resolve().parent;ctx.threads=1;start=time.monotonic()
with (ROOT/'source_basis_snapshot.pkl').open('rb') as stream:data=pickle.load(stream)
report=data['report'];records=data['records'];prime=report['prime'];dimension=report['dimension']
assert report['status']=='SNAPSHOT_SOURCE_INDEPENDENCE_PASS'
assert len(records)==report['snapshot_columns']==report['source_sketch_rank_mod_p']
for i,record in enumerate(records):
    with (ROOT/('vector_%04d.pkl'%i)).open('rb') as stream:saved=pickle.load(stream)
    assert saved=={k:v for k,v in record.items() if k!='source_pivot'}
    stat=Path(record['source_column']).stat();entry=record['snapshot_entry']
    assert (stat.st_size,stat.st_mtime_ns)==(entry['size_bytes'],entry['mtime_ns'])
vectors=[record['source_vector'] for record in records]
matrix=nmod_mat(vectors,prime);assert matrix.rank()==len(records)
coordinates=[record['source_pivot'] for record in records]
minor=nmod_mat([[row[i] for i in coordinates] for row in vectors],prime)
determinant=int(minor.det());assert determinant==report['source_minor_determinant_mod_p']!=0
checks=[]
for index in sorted({0,len(records)//2,len(records)-1}):
    record=records[index];raw=Path(record['source_column']).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==record['source_column_sha256']
    source=pickle.loads(raw);del raw
    assert source['descriptor']==record['descriptor']
    vector=[0]*dimension;coefficients={}
    for word,value in source['polynomial'].items():
        digest=hashlib.blake2b(repr(word).encode(),digest_size=12).digest()
        position=int.from_bytes(digest[:8],'little')%dimension
        sign=1 if digest[8]&1 else -1
        key=(value.numerator,value.denominator)
        if key not in coefficients:coefficients[key]=key[0]*pow(key[1],-1,prime)%prime
        vector[position]=(vector[position]+sign*coefficients[key])%prime
    assert vector==record['source_vector']
    checks.append({'index':index,'source_terms':len(source['polynomial']),'complete_source_sketch_match':True})
    del source,vector
result={'passed':True,'cached_vectors_read_back':len(records),'independent_FLINT_rank':matrix.rank(),
    'source_minor_determinant_mod_p':determinant,'complete_source_sketch_checks':checks,
    'all_original_snapshot_stats_unchanged':True,'seconds':time.monotonic()-start,
    'new_Q0_images_computed':0}
(ROOT/'readback.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
