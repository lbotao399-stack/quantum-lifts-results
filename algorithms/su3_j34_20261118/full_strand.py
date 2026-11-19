"""Put the four strand maps and the 3,900 reduced source rows on identical full target frames.

Input rows are trace recipes. Field arithmetic is exact. Select full-frame columns,
retain the actual coefficients of the *original* pivot recipes, save residuals on
all unselected coordinates, and certify independent rows of the quotient.
"""
import argparse, gzip, hashlib, json, sys
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as maps

P=1000003
PLAN='1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c'
# task, strand source, original compressed recipe archive, original square-matrix hash
STRANDS={7:(4,13,'914e51faa3108c0a049885238e4783b2e7da432bb0f270deff437ab2b55e00f3', 'b0a61e7d5dc4f6b8bb4e8bb0fae487a9128af78e57309d52724ca2533f0afe9e'),
         8:(5,14,'3997d22f80ffac3b1de3ad8e283c9714e3a84ae1a3a48664c5fb9e80673ba98f','fdbea0e20bff3e7974359208acc4d33a179538eed289e4af02438cf9751ad526'),
         9:(6,15,'2fddffd042635fd879509e07f23fd153d9a11e03510b76ca1109581ca7f29dfc','16849e4643dfb291d12f870a145907365c7ad6af5e2e00a6a8caaf659708b741'),
         12:(7,18,'3bb461ff5d6b5cb4b0fb8f0e66c30f29ca545489825d14b81841690f2d894230','e04ff3b823e44530105e40a51213b6dafad711eb85c9713c6c4fe6386dd91cc6')}
EXTRA='7720fa7593cbb6c47789914a4b611073fcb22f764a254a8ec0d8fa5ba53b4527'

def fsha(f):
 h=hashlib.sha256()
 with f.open('rb') as b:
  for v in iter(lambda:b.read(2**20),b''): h.update(v)
 return h.hexdigest()
def packed(a):return np.ascontiguousarray(a,dtype='<u4')
def harray(a):return hashlib.sha256(packed(a).tobytes()).hexdigest()
def putjson(f,o):
 k=f.with_suffix('.part');k.write_text(json.dumps(o,indent=2)+'\n');k.replace(f)
def write_array(path,a):
 k=path.with_suffix('.part.npy');np.save(k,packed(a),allow_pickle=False);k.replace(path);return fsha(path)
def readwords(file,expected):
 assert fsha(file)==expected;c=json.loads(gzip.decompress(file.read_bytes()));w=[maps.decode(row) for i,row in c['rows'] if i==0];assert len(w)==len(c['rows'])==len(set(w))==len(c['columns'])==len(set(c['columns'])) and c['determinant'] and c['map']['prime']==P and c['map']['plan_sha256']==PLAN;return c,w

def init(runtime,plan,t,s):maps.setup(Path(runtime).resolve(),plan,[t],[s],P)
def explicit(w,fr,s,t):
 from native_su_values_v1 import evaluate
 from sp_full_profiles import profile
 expr=maps.project(maps.algebra.parity_project({w:1},1),maps.charge,P)
 out=maps.algebra.differential(expr);keep={v:int(c)%P for v,c in out.items() if c%P and profile(tuple(l for trace in v for l in trace))==maps.pattern(t)}
 return evaluate([keep],fr)[0]

def build(path,words,plan,t,s,runtime,workers):
 from concurrent.futures import ProcessPoolExecutor
 init(runtime,plan,t,s);f=maps.frames[0];width=f.ncols;args=(runtime,plan,t,s);tmp=path.with_suffix('.part.npy');rows=8
 if path.exists():
  receipt=json.loads(path.with_suffix('.json').read_bytes());assert receipt['sha256']==fsha(path) and receipt['recipes_sha256']==harray_words(words) and receipt['source']==s and receipt['target']==t and receipt['plan_sha256']==PLAN
 else:
  data=np.lib.format.open_memmap(tmp,mode='w+',dtype='<u4',shape=(len(words),width))
  with ProcessPoolExecutor(max_workers=workers,initializer=init,initargs=args) as pool:
   for k in range(0,len(words),rows*2*workers):
    ids=list(range(k,min(len(words),k+rows*2*workers),rows)); tasks=[[(0,w) for w in words[j:j+rows]] for j in ids]
    for j,v in zip(ids,pool.map(maps.values,tasks)):data[j:j+len(v)]=v
  data.flush();del data
  saved=np.load(tmp,mmap_mode='r',allow_pickle=False);checks=[]
  for j,ip in [(0,0),(len(words)//2,f.P//2),(len(words)-1,f.P-1)]:
   one=copy(f);one.P=1;one.ncols=f.cpp;one.L={l:(d,a[ip:ip+1]) for l,(d,a) in f.L.items()};np.testing.assert_array_equal(saved[j,ip*f.cpp:(ip+1)*f.cpp],explicit(words[j],one,plan['sources'][s],plan['targets'][t]));checks.append([j,ip,f.cpp])
  del saved;tmp.replace(path);receipt=dict(source=s,target=t,plan_sha256=PLAN,prime=P,rows=len(words),width=width,recipes_sha256=harray_words(words),sha256=fsha(path),direct_word_checks=checks);putjson(path.with_suffix('.json'),receipt)
 return np.load(path,mmap_mode='r',allow_pickle=False),receipt

def harray_words(words):return hashlib.sha256(json.dumps(words,separators=(',',':')).encode()).hexdigest()

def product(l,r):
 from exact_blas_rows import matmul_mod
 return matmul_mod(l,r,P)

def factor(x,preferred,panel=80):
 """Pivot columns in given row order; store multipliers of forward row operations."""
 n,w=x.shape; z=packed(x).copy();L=np.zeros((n,n),dtype='<u4');chosen=[];panels=[];det=1
 for i in range(0,n,panel):
  j=min(i+panel,n);h=j-i; block=z[i:j];used=set(chosen);cols=[k for k in preferred[i:j] if k not in used];M=nmod_mat(block[:,cols].tolist(),P) if len(cols)==h else None;d=int(M.det()) if M is not None else 0
  if not d:
   available=[k for k in range(w) if k not in used];a,r=nmod_mat(block[:,available].tolist(),P).rref();assert r==h,'the full target frame does not separate the stored recipes';cols=[];k=0
   for row in range(h):
    while not a[row,k]:k+=1
    cols.append(available[k]);k+=1
   M=nmod_mat(block[:,cols].tolist(),P);d=int(M.det());assert d
  T=packed(M.inv().tolist());u=packed(product(T,block));np.testing.assert_array_equal(u[:,cols],np.eye(h,dtype='<u4'));z[i:j]=u
  for b in range(j,n,320):
   v=z[b:b+320];c=v[:,cols].copy();L[b:b+len(v),i:j]=c;v[:]=(v.astype(np.int64)-product(c,u))%P
  chosen.extend(cols);panels.append((i,j,cols,T,d));det=det*d%P
 assert len(set(chosen))==n and det;return z,L,panels,chosen,det

def project(data,U,L,panels):
 y=packed(data).copy();total=U.shape[0];beta=np.zeros((len(y),total),dtype='<u4')
 for i,j,cols,T,d in panels:
  b=y[:,cols].copy();beta[:,i:j]=b;y[:]=(y.astype(np.int64)-product(b,U[i:j]))%P
 for i,j,cols,T,d in reversed(panels):
  c=packed(product(beta[:,i:j],T));beta[:,i:j]=c
  if i:beta[:,:i]=(beta[:,:i].astype(np.int64)-product(c,L[i:j,:i]))%P
 return beta,y

def minor(data,cols):
 h=hashlib.sha256()
 for b in range(0,len(data),32):h.update(packed(data[b:b+32,cols]).tobytes())
 return h.hexdigest()

def quotient_rank(arr,complement,labels):
 from su_modular_rowspace_v1 import RowSpace
 space=RowSpace(len(complement),P);selected=[]
 for b in range(0,len(arr),32):selected.extend(b+k for k in space.extend(arr[b:b+32,complement]))
 cols=[complement[int(k)] for k in space.pivots];sub=arr[np.ix_(selected,cols)];d=int(nmod_mat(sub.tolist(),P).det()) if selected else 1;assert d
 return dict(order=len(selected),determinant=d,columns=cols,rows=[labels[i] for i in selected],entries_sha256=harray(sub))

def run(a):
 t=a.target;task,s,sh,oldhash=STRANDS[t];assert fsha(a.plan)==PLAN;p=json.loads(a.plan.read_bytes())['incoming'];assert t not in map(int,p['first_rhs_trace_terms']);c,words=readwords(a.base/'minor.json.gz',sh);assert c['map']['sources']==[s] and c['map']['targets']==[t]
 out=a.out/f'target_{t}';out.mkdir(exist_ok=True,parents=True);x,built=build(out/'images.npy',words,p,t,s,a.runtime,a.workers);n,w=x.shape;oldwidth=c['map']['frame_width'];narrow=minor(x,c['columns']);assert (narrow==oldhash) if oldwidth==w else True
 U,L,panels,cols,det=factor(x,c['columns']);assert not(w==oldwidth and cols==c['columns']) or det==c['determinant'];selected_set=set(cols);complement=[k for k in range(w) if k not in selected_set]
 # Read the first Schur stage in its fixed source-row order. A structural zero is a real zero row.
 base=a.reduced/f'target_{t}';f=base/'result.json';old=json.loads(f.read_bytes());assert(old['target'],old['width'],old['prime'],old['plan_sha256'])==(t,w,P,PLAN) and old['minor_sha256']=='fd22fbcd0ae0ada903acca1831b3ad70dc17fbb9570c833a8dd96544d01040e1';blocks=[];labels=[]
 for z in old['rows']:
  sz=z['free'];s0=z['source'];assert sz=={1:501,2:0,4:1845,10:1554}[s0]
  labels.extend([s0,j] for j in range(sz))
  if z['structural_zero']:blocks.append((sz,None));continue
  file=base/z['path'];assert fsha(file)==z['sha256'];v=np.load(file,mmap_mode='r',allow_pickle=False);assert v.shape==(sz,w);blocks.append((sz,v))
 assert len(labels)==3900 and [z['source'] for z in old['rows']]==[1,2,4,10]
 extra=None
 if t==9:
  ce,we=readwords(a.extra,EXTRA);assert ce['map']['sources']==[6] and ce['map']['targets']==[9] and ce['map']['frame_width']==w
  xe,extra=build(out/'source6_images.npy',we,p,t,6,a.runtime,a.workers);blocks.append((len(we),xe));labels.extend([6,j] for j in range(len(we)))
 N=len(labels);cf=out/'coefficients.npy';res=out/'residual.npy';c_tmp=cf.with_suffix('.part.npy');r_tmp=res.with_suffix('.part.npy');co=np.lib.format.open_memmap(c_tmp,mode='w+',dtype='<u4',shape=(N,n));y=np.lib.format.open_memmap(r_tmp,mode='w+',dtype='<u4',shape=(N,w));off=0
 for size,block in blocks:
  for b in range(0,size,160):
   k=min(size-b,160);cut=slice(off+b,off+b+k)
   if block is None:co[cut]=0;y[cut]=0;continue
   raw=block[b:b+k];C,Y=project(raw,U,L,panels);assert not Y[:,cols].any();np.testing.assert_array_equal((product(C,x)+Y)%P,raw);co[cut]=C;y[cut]=Y
  off+=size
 assert off==N;co.flush();y.flush();del co,y;c_tmp.replace(cf);r_tmp.replace(res);y=np.load(res,mmap_mode='r',allow_pickle=False);q=quotient_rank(y,complement,labels)
 result=dict(prime=P,target=t,strand_source=s,narrow_width=oldwidth,full_width=w,minor_sha256=sh,full_images=built,old_minor_entries_sha256=oldhash,full_frame_old_columns_sha256=narrow,columns=cols,determinant=det,panels=[d for i,j,k,T,d in panels],side_input_sha256=fsha(f),extra_source6=extra,projected_rows=N,coefficient_file_sha256=fsha(cf),residual_file_sha256=fsha(res),residual_nonzero=int(np.count_nonzero(y)),reconstructed_coordinates=N*w,quotient_minor=q)
 putjson(out/'result.json',result);print(json.dumps({k:result[k] for k in ('target','determinant','residual_nonzero','reconstructed_coordinates')}),q['order'],flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--target',required=True,type=int,choices=STRANDS);ap.add_argument('--plan',type=Path,default=maps.HERE/'classical_order1.json');ap.add_argument('--base',type=Path,required=True);ap.add_argument('--reduced',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--extra',type=Path,default=maps.HERE/'other_maps/task_7/minor.json.gz');ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--workers',type=int,default=3);a=ap.parse_args();assert 1<=a.workers<=8 and __debug__;sys.path.insert(0,str(a.runtime.resolve()));run(a)
