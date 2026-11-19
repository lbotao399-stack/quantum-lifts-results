from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,os,subprocess,sys
out=Path('output');base=Path.home()/'codex_n4_lifts_20261118';root=base/'su3_prime1000033_v1_20261118';rt=root/'runtime/work/sp3-native.nosync';shard=int(os.environ['SLURM_ARRAY_TASK_ID']);assert 0<=shard<8 and int(os.environ['SLURM_CPUS_PER_TASK'])==32
work=[]
for eq,n in [(1,'su3_correction_profiles_v1.json'),(2,'su3_second_correction_profiles_v1.json')]:
 plan=json.loads((rt/n).read_bytes());work.extend((t['dimension'],eq,t['index']) for t in plan['targets'])
# Greedy by cost, then run the smaller profiles first within a shard to return useful controls promptly.
lanes=[[] for _ in range(32)];cost=[0]*32
for d,e,i in sorted(work,reverse=True):
 k=min(range(32),key=lambda j:cost[j]);cost[k]+=d*d;lanes[k].append((e,i))
assert len({t for l in lanes for t in l})==len(work)==50

def lane(j):
 for eq,i in reversed(lanes[j]):
  print(json.dumps(dict(begin=[eq,i],lane=j)),flush=True)
  subprocess.run([sys.executable,'-u','reproject.py','--root',str(root),'--first-zip',str(base/'scratch/su3_j30_first_recipes_20261119.ssh.zip'),'--out',str(out),'--equation',str(eq),'--index',str(i),'--workers','8','--batch','16'],check=True)
with ThreadPoolExecutor(4) as pool:list(pool.map(lane,range(shard*4,shard*4+4)))
