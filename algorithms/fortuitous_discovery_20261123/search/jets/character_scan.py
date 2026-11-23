from pathlib import Path
import os,sys,json,time,signal
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
os.environ['N4_TRACE_RUNTIME']=str(ROOT/'upstream/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync')
sys.path.insert(0,str(ROOT/'upstream/algorithms/dual_singlets_20261118'))
from trace_descent import dimensions
start=time.monotonic()
signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('bounded character scan')))
signal.alarm(240)
records=[]
weights=[]
for s in range(6,12):
 for a in range(1,5):
  for b in range(1,a+1):
   k=s-a-b
   if 1<=k<=7:weights.append((a,b,k,0,0))
try:
 for w in weights:
  row={'weight':w,'cochains':{}}
  for g,n in [('Sp',4),('SO',7),('Sp',6),('Sp',20)]:
   ds={L:sum(dimensions(w,L,n,g).values()) for L in range(max(2,sum(w)-3),sum(w)+1)}
   row['cochains'][g+str(n)]=ds
  records.append(row)
  print(json.dumps(row),flush=True)
  (HERE/'character_scan.json').write_text(json.dumps({'records':records,'seconds':time.monotonic()-start},indent=2))
except TimeoutError:
 print('TIME LIMIT',flush=True)
finally:
 signal.alarm(0)
