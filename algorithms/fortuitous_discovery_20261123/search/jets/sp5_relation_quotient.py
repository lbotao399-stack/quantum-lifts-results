"""Finite relation quotient, reduced against the actual stable Q0 image in trace blocks."""
from pfaffian_relations import *
import numpy as np
from flint import nmod_mat,fmpq_mat
from fractions import Fraction
import argparse,signal,hashlib

def pivots_from_rref(a,rank):
 ans=[];j=0
 for i in range(rank):
  while not a[i,j]:j+=1
  ans.append(j);j+=1
 return ans

def modmatrix(rows,words,p):
 ids={w:i for i,w in enumerate(words)};a=nmod_mat(len(rows),len(words),p)
 for i,row in enumerate(rows):
  for w,c in row.items():
   if w in ids:a[i,ids[w]]=int(c)%p
 return a

def row_rank(rows,p):
 words=sorted({w for row in rows for w in row});a=modmatrix(rows,words,p)
 return a.rank()

def components(images):
 parent=list(range(len(images)));owners={}
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for i,row in enumerate(images):
  for w in row:
   if w in owners:parent[find(i)]=find(owners[w])
   else:owners[w]=i
 groups={}
 for i,row in enumerate(images):
  if row:groups.setdefault(find(i),[]).append(i)
 return sorted(groups.values(),key=lambda v:(len(v),v[0]))

def reconstruct(value,p):
 if not value:return Fraction(0)
 r,s=p,int(value);t,u=0,1;bound=math.isqrt((p-1)//2)
 while abs(s)>bound:
  q=r//s;r,s=s,r-q*s;t,u=u,t-q*u
 if not u or abs(u)>bound or math.gcd(s,u)!=1:return None
 result=Fraction(s,u)
 if (result.numerator-value*result.denominator)%p:return None
 return result

def main(p,seconds,target_name):
 start=time.monotonic();work=HERE/f'quotient_sp5_{target_name.split(".")[0]}_p{p}';work.mkdir(exist_ok=True)
 formal=HERE.parent/'sp5_universal'
 source=pickle.loads((formal/'formal_339_L13.pkl').read_bytes())
 images=pickle.loads((formal/'q_339_L13.pkl').read_bytes())
 projection=pickle.loads((formal/'pivot_words_formal_339_L14_p2147483647.pkl').read_bytes())
 assert projection['rank']==17745 and projection['block_determinant_product_mod_prime']
 keep=set(projection['pivot_words'])
 images=[{w:c for w,c in row.items() if w in keep} for row in images]
 selected_global=set()
 for block_path in (formal/'rank_339_L13_p2147483647').glob('block_*.pkl'):
  block_data=pickle.loads(block_path.read_bytes());selected_global.update(block_data['selected_source_rows'])
 assert len(selected_global)==14590
 last=pickle.loads((formal/'rank_339_L13_last_projected_partial.pkl').read_bytes())
 full_relations,_=pickle.loads((HERE/target_name).read_bytes())
 relations=[{w:c for w,c in row.items() if w in keep} for row in full_relations]
 previous,_=pickle.loads((HERE/'exact_relations_sp5_L13.pkl').read_bytes())
 before,_=pickle.loads((HERE/'exact_relations_sp5_L12.pkl').read_bytes())
 nr=len(relations)
 assert (len(source),len(previous),len(before))==(37085,188,67) and 121<=nr<=183
 Rhash=hashlib.sha256(pickle.dumps(relations,protocol=4)).hexdigest()
 groups=components(images);owner={w:i for i,g in enumerate(groups) for j in g for w in images[j]}
 restricted=[[{} for r in relations] for group in groups];extra=[{} for r in relations]
 for i,row in enumerate(relations):
  for w,c in row.items():
   if w in owner:restricted[owner[w]][i][w]=c
   else:extra[i][w]=c
 print('stable blocks',len(groups),'largest',max(map(len,groups)),flush=True)
 residuals=[];coefficients=[];basis_ids=[];report={'prime':p,'relation_sha256':Rhash,'blocks':[]}
 for k,group in enumerate(groups):
  path=work/f'block_{k}.pkl'
  if path.exists():
   record,residual,C,ind=pickle.loads(path.read_bytes());assert record['relation_sha256']==Rhash
  else:
   words=sorted({w for j in group for w in images[j]});wid={w:i for i,w in enumerate(words)}
   ind=[j for j in group if j in selected_global];rank=len(ind)
   chosen=modmatrix([images[j] for j in ind],words,p)
   if group==last['source_rows']:
    pc=[wid[last['projected_words'][j]] for j in last['projected_pivot_columns']]
    assert rank==last['rank']==3507
   else:
    rr,rcheck=chosen.rref();assert rcheck==rank
    pc=pivots_from_rref(rr,rank);del rr
   T=nmod_mat([[chosen[i,j] for j in pc] for i in range(rank)],p)
   B=modmatrix(restricted[k],words,p)
   Bp=nmod_mat([[B[i,j] for j in pc] for i in range(len(relations))],p)
   Cmat=T.transpose().solve(Bp.transpose()).transpose()
   residual=np.asarray((B-Cmat*chosen).tolist(),dtype=np.int64)
   active=np.flatnonzero(np.any(residual,axis=0));residual=residual[:,active]
   C=np.asarray(Cmat.tolist(),dtype=np.int64)
   record={'block':k,'source_rows':len(group),'words':len(words),'rank':rank,'residual_columns':len(active),'relation_sha256':Rhash,'seconds':time.monotonic()-start}
   temp=path.with_suffix('.part');temp.write_bytes(pickle.dumps((record,residual,C,ind),protocol=4));temp.replace(path)
  residuals.append(residual);coefficients.append(C);basis_ids+=ind;report['blocks'].append(record)
  (work/'progress.json').write_text(json.dumps(report,indent=2))
  print('block',k,'rank',record['rank'],'residual_cols',record['residual_columns'],'seconds',time.monotonic()-start,flush=True)
 extra_words=sorted({w for row in extra for w in row});extra_matrix=np.asarray(modmatrix(extra,extra_words,p).tolist(),dtype=np.int64)
 residuals.append(extra_matrix)
 residual=np.column_stack(residuals);C=np.column_stack(coefficients)
 assert sum(r['rank'] for r in report['blocks'])==14590 and C.shape==(nr,14590)
 rr,resrank=nmod_mat(residual.T.tolist(),p).rref();pc=pivots_from_rref(rr,resrank);free=sorted(set(range(nr))-set(pc))
 kernel=[]
 for f in free:
  v=np.zeros(nr,dtype=np.int64);v[f]=1;v[pc]=[(-int(rr[i,f]))%p for i in range(resrank)];kernel.append(v)
 print('target intersection',len(kernel),'residual rank',resrank,flush=True)
 qp=HERE/'sp5_q13_identity_rows';qb=HERE/'exact_q_relations_sp5_L12.pkl'
 qprevious=[]
 for i,row in enumerate(previous):
  h,image=pickle.loads((qp/f'row_{i}.pkl').read_bytes());assert h==hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest();qprevious.append({w:c for w,c in image.items() if w in keep})
 qbefore=pickle.loads(qb.read_bytes()) if qb.exists() else [m.q(row) for row in before]
 assert len(qprevious)==len(previous) and len(qbefore)==len(before)
 nilpotence=json.loads((formal/'q2_339_L12.json').read_text())
 assert nilpotence['complete'] and nilpotence['exact_integer'] and nilpotence['checked']==38160
 words=sorted({w for row in relations for w in row}|{w for row in qprevious for w in row})
 R=modmatrix(relations,words,p);rref,rank=R.rref();assert rank==nr
 piv=pivots_from_rref(rref,rank);del rref
 T=nmod_mat([[R[i,j] for j in piv] for i in range(nr)],p)
 Q=modmatrix(qprevious,words,p);Qp=nmod_mat([[Q[i,j] for j in piv] for i in range(188)],p)
 coords=T.transpose().solve(Qp.transpose()).transpose();assert (Q-coords*R).rank()==0
 Qrr,rankQ=coords.rref();qpc=pivots_from_rref(Qrr,rankQ)
 rank_before=row_rank(qbefore,p)
 assert rank_before==67 and rankQ==121
 (work/'screen.json').write_text(json.dumps({'prime':p,'stable_Q11_rank':14590,'target_relations':nr,'stable_image_intersection':len(kernel),'rational_Q_I11_rank':rankQ,'modular_quotient_dimension':len(kernel)-rankQ,'seconds':time.monotonic()-start},indent=2))
 candidates=[];Cfield=nmod_mat(C.tolist(),p)
 for z in kernel:
  reduced=[int(x) for x in z]
  for i,c in enumerate(qpc):
   factor=reduced[c]
   if factor:reduced=[(x-factor*int(Qrr[i,j]))%p for j,x in enumerate(reduced)]
  if any(reduced):
   pre=np.asarray((nmod_mat([z.tolist()],p)*Cfield).tolist()[0],dtype=np.int64)
   candidates.append((int(np.count_nonzero(pre)),z,pre))
 assert candidates
 candidates.sort(key=lambda v:v[0]);nnz,z,pre=candidates[0]
 np.savez_compressed(work/'modular_candidate.npz',relation_coefficients=z,preimage_coefficients=pre,basis_rows=np.asarray(basis_ids,dtype=np.int32),quotient_residual=residual)
 report.update(stable_outgoing_rank=14590,stable_incoming_rank=22495,stable_centre_dimension=37085,stable_rational_cohomology=0,relation_dimension=nr,relation_intersection_dimension=len(kernel),relation_image_rank=rankQ,preceding_relation_image_rank=rank_before,partial_relation_quotient_dimension=len(kernel)-rankQ,complete_target_relation_space=nr==183,modular_candidate_source_terms=nnz,extra_relation_words=len(extra_words),seconds=time.monotonic()-start)
 rational=[reconstruct(int(c),p) for c in pre]
 failed=sum(c is None for c in rational);report['rational_reconstruction_failures']=failed
 if not failed:
  denominator=math.lcm(*(c.denominator for c in rational));poly=Counter()
  for i,c in zip(basis_ids,rational):
   if c:
    for w,n in source[i].items():poly[w]+=int(c*denominator)*n
  poly={w:c for w,c in poly.items() if c};variation=m.q(poly)
  rational_minor=fmpq_mat([[full_relations[i].get(words[j],0) for i in range(nr)] for j in piv])
  rhs=fmpq_mat([[variation.get(words[j],0)] for j in piv]);ratcoords=rational_minor.solve(rhs)
  rc=[Fraction(str(c)) for c in ratcoords.entries()];den=math.lcm(*(c.denominator for c in rc));comb=Counter()
  for c,row in zip(rc,full_relations):
   if c:
    for w,n in row.items():comb[w]+=int(c*den)*n
  expected={w:c*den for w,c in variation.items()};comb={w:c for w,c in comb.items() if c}
  exact=expected==comb;report['exact_closure']=exact
  if exact and nr==183:report['rational_cohomology_dimension']=len(kernel)-rankQ
  if exact:
   cv=[c.numerator*pow(c.denominator,-1,p)%p for c in rc]
   test=nmod_mat(coords.tolist()+[cv],p).rank();assert test==rankQ+1
   report.update(rational_fortuitous_class_lower_bound=1,cycle_terms=len(poly),cycle_variation_terms=len(variation),cycle_common_denominator=str(denominator),nonboundary_augmented_rank=test,finite_closure_relation_coefficients=[str(c) for c in rc])
   (work/'rational_cycle.pkl').write_bytes(pickle.dumps({'polynomial':poly,'original_scale_denominator':denominator,'relation_coefficients':rc},protocol=4))
 (work/'result.json').write_text(json.dumps(report,indent=2))
 print({k:v for k,v in report.items() if k not in ('blocks','finite_closure_relation_coefficients')},flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--prime',type=int,default=2147483647);a.add_argument('--seconds',type=int,default=300);a.add_argument('--target',default='exact_relations_sp5_L14.pkl');args=a.parse_args()
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('quotient reduction cap')));signal.alarm(args.seconds)
 main(args.prime,args.seconds,args.target)
