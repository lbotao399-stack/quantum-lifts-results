from fast_paths import *
import argparse
import hashlib
import json
import pickle
import time

parser=argparse.ArgumentParser()
parser.add_argument('--rank',type=int,choices=[4,5],required=True)
parser.add_argument('--samples',type=int,default=12)
args=parser.parse_args()
rank=args.rank
start=time.monotonic()
records=[]
for length in range(2*rank+2,2*rank+5):
    filename=pc.HERE/f'exact_relations_sp{rank}_L{length}.pkl'
    rows,recipes=pickle.loads(filename.read_bytes())
    raw=[(i,r) for i,r in enumerate(recipes) if not isinstance(r[0],str)]
    paths=[(i,r[1][0][0]) for i,r in enumerate(recipes) if r[0] in ('paired','paths')]
    for index,recipe in raw[:args.samples]:
        letters,remaining=recipe
        t=time.monotonic();old=pf.identity(letters,rank,remaining);t1=time.monotonic()
        new=identity(letters,rank,remaining);t2=time.monotonic()
        assert old==new,('raw',rank,length,index)
        records.append(dict(length=length,index=index,kind='raw',terms=len(new),
                            endpoint_seconds=t1-t,cycle_seconds=t2-t1))
    for index,paths in paths[:args.samples]:
        t=time.monotonic();old=pc.expand_paths(paths,rank);t1=time.monotonic()
        new=expand_paths(paths,rank);t2=time.monotonic()
        assert old==new,('paths',rank,length,index)
        records.append(dict(length=length,index=index,kind='paths',terms=len(new),
                            endpoint_seconds=t1-t,cycle_seconds=t2-t1))
    print('compared',rank,length,len(records),'seconds',time.monotonic()-start,flush=True)
result=dict(rank=rank,exact_comparisons=len(records),all_integer_coefficients_equal=True,
            endpoint_seconds=sum(r['endpoint_seconds'] for r in records),
            cycle_seconds=sum(r['cycle_seconds'] for r in records),records=records)
(HERE/f'frontier_validation_sp{rank}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)
