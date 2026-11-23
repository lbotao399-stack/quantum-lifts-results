"""Cold, same-input coefficient checks against all 13!! matchings."""
from fast_paths import *
import hashlib,json,pickle,time

cache_native_canonicalization()
_,recipes=pickle.loads((pc.HERE/'exact_relations_sp5_L12.pkl').read_bytes())
A=(0,0,(0,));records=[]
for index in (0,33,66):
    paths=tuple(recipes[index][1][0][0])+((A,A),)
    pf.pairings.cache_clear();pf.skeleton.cache_clear()
    pf.m.so_canon.cache_clear();pf.m.fw.canon_word.cache_clear()
    start=time.perf_counter();old=pc.expand_paths(paths,6);old_seconds=time.perf_counter()-start
    pf.m.so_canon.cache_clear();pf.m.fw.canon_word.cache_clear();cycle.cache_clear()
    start=time.perf_counter();new=expand_paths(paths,6);new_seconds=time.perf_counter()-start
    assert old==new
    record=dict(seed_recipe_index=index,path_recipe=paths,matching_count=135135,terms=len(new),
                exact_integer_coefficients_equal=True,endpoint_seconds=old_seconds,
                cycle_seconds=new_seconds,cold_native_canonicalization_caches=True,
                semantic_sha256=hashlib.sha256(json.dumps([(c,pf.m.encode(w)) for w,c in sorted(new.items())],separators=(',',':')).encode()).hexdigest())
    records.append(record)
    print(json.dumps({k:v for k,v in record.items() if k!='path_recipe'}),flush=True)
(HERE/'benchmark_sp6.json').write_text(json.dumps(dict(records=records),indent=2)+'\n')
