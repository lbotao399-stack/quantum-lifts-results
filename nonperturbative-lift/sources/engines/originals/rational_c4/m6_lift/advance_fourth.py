"""Read back the fourth joint primitive and prepare the last full-Q0 target."""
from pathlib import Path
import sys,os,pickle,json,hashlib,time,gc,runpy,shutil
SCRIPT_ROOT=Path(__file__).resolve().parent;ROOT=SCRIPT_ROOT/'compact_first/lift'
JOINT=ROOT/'joint4'
sys.path.insert(0,str(SCRIPT_ROOT.parent))
import exact_product as E,integer_qn
start=time.time()
def load(path):
    with path.open('rb') as f:return pickle.load(f)
def save(path,data):
    tmp=path.with_suffix('.tmp')
    with tmp.open('wb') as f:pickle.dump(data,f)
    tmp.replace(path)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(1<<20):h.update(block)
    return h.hexdigest()
report=json.loads((JOINT/'summary.json').read_text())
assert report['passed']
source=JOINT/'component4_exact.pkl';digest=sha(source);assert digest==report['sha256']
p=load(source);targets=load(JOINT/'targets_exact.pkl')
assert not E.sub(integer_qn.q0(E.T,p),targets['Q0'])
assert not load(JOINT/'Q0_full_residual.pkl')
assert not load(JOINT/'cartan_Q1_full_residual.pkl')
adopt=ROOT/'component4_exact.pkl'
if not adopt.exists():shutil.copyfile(source,adopt)
assert sha(adopt)==digest
terms=len(p);del p;gc.collect()
os.environ.update(M6_RUN_ROOT=str(ROOT),M6_DIRECT_FANS='1',M6_INTEGER_FANS='1')
sys.argv=[str(SCRIPT_ROOT/'component_images.py'),'4']
runpy.run_path(sys.argv[0],run_name='__main__')
q=load(ROOT/'Q1_component4.pkl')
assert not E.sub(E.cartan(q),targets['cartan_Q1'])
del targets;gc.collect()
(ROOT/'level4_independent_readback.json').write_text(json.dumps(dict(
    passed=True,terms=terms,sha256=digest,full_Q0_residual_terms=0,
    full_Q1_image_Cartan_residual_terms=0,seconds=time.time()-start),indent=2)+'\n')
target=E.sc(q,-1);del q
for name in ('Q2_component3.pkl','Q3_component2.pkl','Q4_component1.pkl','Q5_component0.pkl'):
    for word,c in load(ROOT/name).items():target[word]=target.get(word,0)-c
target={w:c for w,c in target.items() if c}
save(ROOT/'level5_target.pkl',target)
closure=integer_qn.q0(E.T,target);save(ROOT/'level5_Q0_closure.pkl',closure)
cartan=E.cartan(target);save(ROOT/'level5_Cartan_residual.pkl',cartan)
assert not closure and not cartan,(len(closure),len(cartan))
summary=dict(passed=True,terms=len(target),Q0_closure_terms=0,Cartan_terms=0,
    preceding_sha256=digest,target_sha256=sha(ROOT/'level5_target.pkl'),seconds=time.time()-start,
    scope='Last target is Q0 exact by the proved six-dimensional Cartan isomorphism in (Nb,Nc)=(0,6); explicit primitive pending')
(ROOT/'level5_target_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('FOURTH READBACK AND LAST TARGET PASS',json.dumps(summary),flush=True)
