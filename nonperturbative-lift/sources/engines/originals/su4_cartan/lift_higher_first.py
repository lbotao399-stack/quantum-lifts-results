"""Bounded exact full-Q0 lift of SU4 scalar Cartan candidates k4,5,6."""
import os,sys,pickle,json,itertools,time,random,hashlib
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from flint import fmpq,fmpq_mat
os.environ['SUN']='4'
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'rational_c4'))
import exact_product as E
sys.path.insert(0,str(ROOT/'su3_rational'))
import fast_wick_edges,fast_selected_direct
fast_wick_edges.install(E.T);fast_selected_direct.install(E.T)
T=E.T;k=int(os.environ.get('SU4_K','4'));assert k in (4,5,6)
from full_q0_stream import verify as verify_full_q0
import sampler_q0
dest=OUT/f'k{k}_first';dest.mkdir(exist_ok=True)
print('PID',os.getpid(),'k',k,flush=True)
start=time.time();rng=random.Random(840+k)
search_seconds=float(os.environ.get('SU4_SEARCH_SECONDS','1600'))
saved=pickle.load(open(OUT/'leading_radial.pkl','rb'))
C0,a1=saved['generators'][:2];x,y,z=saved['casimirs']
R=E.sub(saved['generators'][2],E.sc(saved['generators'][3],F(1,4)))
if k==4:
    leading=E.add(E.mul(x,R),E.mul(E.add(E.sc(E.mul(x,x),F(-5,16)),E.sc(z,F(1,2))),C0),E.sc(E.mul(y,a1),F(-2,3)))
elif k==5:
    leading=E.add(E.mul(y,R),E.mul(E.add(E.sc(E.mul(x,x),F(3,16)),E.sc(z,F(-3,4))),a1))
else:
    leading=E.add(E.mul(E.add(E.sc(E.mul(x,x),3),E.sc(z,-8)),R),
        E.mul(E.add(E.mul(x,z),E.sc(E.mul(x,x,x),F(-11,24)),E.sc(E.mul(y,y),F(-4,9))),C0),
        E.sc(E.mul(x,y,a1),F(4,3)))
verify_full_q0(T,leading,{})
target_path=dest/'target.pkl'
leading_bytes=pickle.dumps(leading,protocol=4)
if (dest/'leading.pkl').exists():assert pickle.load(open(dest/'leading.pkl','rb'))==leading
if target_path.exists():target=pickle.load(open(target_path,'rb'))
else:
    if k==6:
        gate=json.loads((OUT/'integer_Q1_gate.json').read_text());assert gate['status']=='PASS'
        sys.path.insert(0,str(ROOT/'rational_c4/m6_lift'))
        import integer_qn
        target=integer_qn.qn(T,1,leading)
    else:target=T.qn(1,leading)
    target_path.write_bytes(pickle.dumps(target,protocol=4))
print('saved full Q1 target',len(target),'seconds',time.time()-start,flush=True)
assert not E.cartan(target)
closure_path=dest/'target_closure.json'
target_hash=hashlib.sha256(target_path.read_bytes()).hexdigest()
leading_hash=hashlib.sha256(leading_bytes).hexdigest()
if closure_path.exists():
    closure=json.loads(closure_path.read_text())
    assert closure['target_sha256']==target_hash and closure['leading_sha256']==leading_hash
    assert closure['status']=='PASS'
    print('REUSED certified full Q0 target closure',flush=True)
else:
    checks=verify_full_q0(T,target,{},verbose=True)
    closure={'status':'PASS','target_sha256':target_hash,'leading_sha256':leading_hash,'full_Q0_channels':checks}
    temporary=closure_path.with_suffix('.tmp');temporary.write_text(json.dumps(closure,indent=2)+'\n');temporary.replace(closure_path)
print('leading',len(leading),'target',len(target),flush=True)
(dest/'leading.pkl').write_bytes(leading_bytes)

def scalar(p):
    j=1
    while True:
        raised=E.Jp(p)
        if not raised:return p
        p=E.sub(p,E.sc(E.S.Jm(raised),F(1,j*(j+1))));j+=1
        assert j<12

def types(Nb):
    jets=list(itertools.product(range(3),repeat=2))
    for bs in itertools.combinations_with_replacement(jets,Nb):
        left=(2-sum(j[0] for j in bs),2-sum(j[1] for j in bs))
        if min(left)<0:continue
        for cs in itertools.combinations_with_replacement([j for j in jets if sum(j)],2):
            if sum(j[0] for j in cs)==left[0] and sum(j[1] for j in cs)==left[1]:
                yield [('b',)+j for j in bs]+[('c',)+j for j in cs]

def cycles(permutation):
    seen=set();out=[]
    for i in range(len(permutation)):
        if i in seen:continue
        cycle=[];j=i
        while j not in seen:seen.add(j);cycle.append(j);j=permutation[j]
        out.append(tuple(cycle))
    return tuple(out)

@lru_cache(maxsize=4)
def trace(letters):return E.cyc_poly(letters)
@lru_cache(maxsize=None)
def gaps(n,total):return tuple(t for t in itertools.product(range(4),repeat=n) if sum(t)==total)
@lru_cache(maxsize=32)
def power(i,j,l):return E.mul(*([x]*i+[y]*j+[z]*l))
options=[]
for letters in types(k-1):
    special=tuple(t for t in letters if t!=('b',0,0));bare=len(letters)-len(special)
    for i in range(bare//2+1):
        for j in range(bare//3+1):
            for l in range(bare//4+1):
                remaining=bare-2*i-3*j-4*l
                if remaining>=0 and gaps(len(special),remaining):options.append((special,(i,j,l),remaining))

prime=2147483647;dimension=1021;hashes={}
def sketch(poly):
    out=[F(0)]*dimension
    for word,c in poly.items():
        if word not in hashes:
            h=hashlib.blake2b(repr(word).encode(),digest_size=12).digest()
            if len(hashes)>250000:hashes.clear()
            hashes[word]=(int.from_bytes(h[:8],'little')%dimension,1 if h[8]&1 else -1)
        i,sign=hashes[word];out[i]+=sign*c
    return out

def canonical_descriptor(desc):
    special,powers,gs,cycs=desc;factors=[]
    for cycle in cycs:
        word=tuple(u for i in cycle for u in ([('b',0,0)]*gs[i]+[special[i]]))
        rotations=[]
        for cut in range(len(word)):
            left=sum(u[0]=='c' for u in word[:cut]);right=sum(u[0]=='c' for u in word[cut:])
            rotations.append((word[cut:]+word[:cut],(-1)**(left*right)))
        minimum=min(w for w,sign in rotations)
        signs={sign for w,sign in rotations if w==minimum}
        if len(signs)>1:return None
        factors.append(minimum)
    factors.sort()
    for left,right in zip(factors,factors[1:]):
        if left==right and sum(u[0]=='c' for u in left)%2:return None
    return powers,tuple(factors)
def modular(v):return np.array([c.numerator*pow(c.denominator,-1,prime)%prime for c in v],dtype=np.int64)
rows=[];coordinates=[];vectors=[];descriptors=[]
def reduce(v):
    v=v.copy()
    for pivot,row in rows:
        co=int(v[pivot])
        if co:v=(v-co*row)%prime
    return v
target_vector=sketch(target);target_mod=modular(target_vector)
for path in sorted(dest.glob('column_*.pkl')):
    previous=pickle.load(open(path,'rb'));vector=previous['vector'];remainder=reduce(modular(vector))
    support=np.flatnonzero(remainder);assert len(support)
    pivot=int(support[0]);assert pivot==previous['pivot']
    rows.append((pivot,remainder*pow(int(remainder[pivot]),-1,prime)%prime))
    coordinates.append(pivot);vectors.append(vector);descriptors.append(previous['descriptor'])
if vectors:print('RESUMED independent columns',len(vectors),flush=True)

def exact_solve():
    matrix=fmpq_mat([[fmpq(v[i].numerator,v[i].denominator) for v in vectors] for i in coordinates])
    rhs=fmpq_mat([[fmpq(target_vector[i].numerator,target_vector[i].denominator)] for i in coordinates])
    result=matrix.solve(rhs)
    coefficients=[F(int(result[i,0].numerator),int(result[i,0].denominator)) for i in range(len(vectors))]
    (dest/'exact_coefficients_checkpoint.pkl').write_bytes(pickle.dumps({'coefficients':coefficients,'descriptors':descriptors},protocol=4))
    primitive={}
    for i,co in enumerate(coefficients):
        if co:
            p=pickle.load(open(dest/f'column_{i:04d}.pkl','rb'))['polynomial']
            for word,c in p.items():
                value=primitive.get(word,F(0))+co*c
                if value:primitive[word]=value
                else:primitive.pop(word,None)
    (dest/'primitive_candidate.pkl').write_bytes(pickle.dumps(primitive,protocol=4))
    checks=verify_full_q0(T,primitive,target,require_zero=False,verbose=True)
    residual={i:row['residual_terms'] for i,row in enumerate(checks) if row['residual_terms']}
    print('EXACT SOLVE primitive',len(primitive),'residual channel sizes',residual,flush=True)
    return coefficients,primitive,residual

seen=set(descriptors);seen_canonical={canonical_descriptor(desc) for desc in descriptors};solution=None
state_path=dest/'sampler_state.pkl'
if state_path.exists():
    state=pickle.load(open(state_path,'rb'));seen.update(state['seen']);rng.setstate(state['rng'])
    seen_canonical.update(state.get('seen_canonical',()))
if vectors and not np.any(reduce(target_mod)):
    coefficients,primitive,residual=exact_solve()
    if not residual:solution=(coefficients,primitive)
candidate_pool=None
inventory_path=dest/'descriptor_inventory.pkl'
source_span=None
if k==6 and os.environ.get('SU4_SOURCE_PREFILTER')=='1':
    sys.path.insert(0,str(OUT/'source_prefilter'))
    from cache_snapshot import ModSpan,source_sketch,atomic_pickle
    source_state_path=dest/'source_span_checkpoint.pkl'
    source_data=pickle.load(open(source_state_path if source_state_path.exists() else OUT/'source_prefilter/source_basis_snapshot.pkl','rb'))
    source_span=ModSpan(source_data['rows']);source_records=source_data['records'][:];source_hashes={}
    source_keys={canonical_descriptor(row['descriptor']) for row in source_records}
    for i,desc in enumerate(descriptors):
        if canonical_descriptor(desc) in source_keys:continue
        polynomial=pickle.load(open(dest/f'column_{i:04d}.pkl','rb'))['polynomial']
        source_vector=source_sketch(polynomial,source_hashes)
        assert source_span.insert(source_vector)
        source_records.append({'descriptor':desc,'source_vector':source_vector})
        source_keys.add(canonical_descriptor(desc))
    def save_source_state():
        atomic_pickle(source_state_path,{'rows':source_span.rows,'records':source_records})
    save_source_state()
    print('SOURCE PREFILTER SEED',source_span.rank,'of 248',flush=True)
if k==6 and inventory_path.exists():
    inventory=pickle.load(open(inventory_path,'rb'))['all']
    assert all(canonical_descriptor(desc)==key for key,desc in inventory.items())
    excluded=source_keys if source_span is not None else seen_canonical
    candidate_pool=[desc for key,desc in inventory.items() if key not in excluded]
    rng.shuffle(candidate_pool)
    print('FINITE UNSEEN CANONICAL CANDIDATES',len(candidate_pool),flush=True)
for trial in range(20000):
    if solution is not None or time.time()-start>search_seconds:break
    if source_span is not None and source_span.rank==248:break
    if candidate_pool is not None:
        if not candidate_pool:break
        desc=candidate_pool.pop();special,powers,gs,cycs=desc
    else:
        special,powers,remaining=rng.choice(options);gs=rng.choice(gaps(len(special),remaining))
        cycs=cycles(tuple(rng.sample(range(len(special)),len(special))))
        if any(len(c)==1 and gs[c[0]]==0 for c in cycs):continue
        desc=(special,powers,gs,cycs)
    if source_span is None and desc in seen:continue
    seen.add(desc)
    canonical=canonical_descriptor(desc)
    if canonical is None or (source_span is None and canonical in seen_canonical):continue
    seen_canonical.add(canonical)
    factors=[]
    for cycle in cycs:
        letters=[]
        for i in cycle:letters.extend([('b',0,0)]*gs[i]+[special[i]])
        factors.append(trace(tuple(letters)))
    poly=E.mul(power(*powers),scalar(E.mul(*factors)))
    if not poly:continue
    if source_span is not None:
        source_vector=source_sketch(poly,source_hashes)
        if not source_span.insert(source_vector):continue
        assert source_span.rank<=248
        source_records.append({'descriptor':desc,'source_vector':source_vector})
        source_keys.add(canonical)
    if (OUT/'sampler_q0_gates.json').exists():
        assert json.loads((OUT/'sampler_q0_gates.json').read_text())['status']=='PASS'
        image=sampler_q0.q0(T,poly)
    else:image=E.Q0(poly)
    if source_span is not None:
        print('SOURCE RANK',source_span.rank,'Q0 terms',len(image),flush=True)
    if not image:
        if source_span is not None:save_source_state()
        continue
    vector=sketch(image);remainder=reduce(modular(vector))
    support=np.flatnonzero(remainder)
    if not len(support):
        if source_span is not None:save_source_state()
        continue
    pivot=int(support[0]);rows.append((pivot,remainder*pow(int(remainder[pivot]),-1,prime)%prime))
    coordinates.append(pivot);vectors.append(vector);descriptors.append(desc)
    i=len(vectors)-1
    column_path=dest/f'column_{i:04d}.pkl';column_tmp=column_path.with_suffix('.tmp')
    column_tmp.write_bytes(pickle.dumps({'polynomial':poly,'descriptor':desc,'vector':vector,'pivot':pivot},protocol=4));column_tmp.replace(column_path)
    if source_span is not None:save_source_state()
    temp=state_path.with_suffix('.tmp');temp.write_bytes(pickle.dumps({'seen':seen,'seen_canonical':seen_canonical,'rng':rng.getstate()},protocol=4));temp.replace(state_path)
    (dest/'progress.json').write_text(json.dumps({'rank':len(vectors),'trial':trial,'seconds':time.time()-start})+'\n')
    print('independent',len(vectors),'trial',trial,'poly',len(poly),'image',len(image),'seconds',time.time()-start,flush=True)
    if not np.any(reduce(target_mod)):
        coefficients,primitive,residual=exact_solve()
        if not residual:solution=(coefficients,primitive);break
    if time.time()-start>search_seconds:break
temp=state_path.with_suffix('.tmp')
temp.write_bytes(pickle.dumps({'seen':seen,'seen_canonical':seen_canonical,'rng':rng.getstate()},protocol=4))
temp.replace(state_path)
if solution is None:
    record={'status':'PENDING_BOUNDED_SEARCH','k':k,'selected_Q0_images':len(vectors),'seconds':time.time()-start,'full_rational_primitive_certified':False,'state_preserved':True}
    (dest/'bounded_status.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    sys.exit(3)
coefficients,primitive=solution
assert not E.Jp(primitive) and not E.S.Jm(primitive)
record={'status':'FULL_RATIONAL_FIRST_LIFT_PASS','k':k,'weight':[k+1,k+1],
        'leading_terms':len(leading),'target_terms':len(target),'primitive_terms':len(primitive),
        'selected_Q0_images':len(vectors),'seconds':time.time()-start,
        'full_Q0_primitive_equals_Q1_leading':True,'scalar_E2_dimension':1,'higher_survival':'OPEN'}
(dest/'first_lift.pkl').write_bytes(pickle.dumps({'leading':leading,'Q1_leading':target,'primitive':primitive,
    'coefficients':coefficients,'descriptors':descriptors,'record':record},protocol=4))
(dest/'first_lift.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2),flush=True)
