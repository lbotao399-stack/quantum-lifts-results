"""Export the reconstructed Sp(5) rational operator and an exact separating dual."""
from sp5_relation_quotient import *
import gzip

start=time.monotonic();p=2147483647
work=HERE/'quotient_sp5_target_sp5_snapshot_181_p2147483647'
cycle=pickle.loads((work/'rational_cycle.pkl').read_bytes())
poly=cycle['polynomial'];scale=cycle['original_scale_denominator'];coeffs=cycle['relation_coefficients']
relations,recipes=pickle.loads((HERE/'target_sp5_snapshot_181.pkl').read_bytes())
source,_=pickle.loads((HERE/'exact_relations_sp5_L13.pkl').read_bytes())
images=[]
for i,row in enumerate(source):
 h,image=pickle.loads((HERE/'sp5_q13_identity_rows'/f'row_{i}.pkl').read_bytes())
 assert h==hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest()
 images.append(image)
(HERE/'exact_q_relations_sp5_L13.pkl').write_bytes(pickle.dumps(images,protocol=4))
var=m.q(poly)
den=math.lcm(*(c.denominator for c in coeffs));comb=Counter()
for c,row in zip(coeffs,relations):
 for w,n in row.items():comb[w]+=int(c*den)*n
assert {w:c for w,c in comb.items() if c}=={w:c*den for w,c in var.items()}
for raising in (1,2,3):
 raised=Counter()
 for w,c in poly.items():
  for image,n in m.raise_word(w,raising).items():raised[image]+=c*n
 assert not any(raised.values())
words=sorted({w for row in relations+images for w in row})
R=modmatrix(relations,words,p);rr,rank=R.rref();assert rank==181
pc=pivots_from_rref(rr,rank)
T=nmod_mat([[R[i,j] for j in pc] for i in range(181)],p);detR=int(T.det());assert detR
Q=modmatrix(images,words,p);Qp=nmod_mat([[Q[i,j] for j in pc] for i in range(188)],p)
M=T.transpose().solve(Qp.transpose()).transpose();assert (Q-M*R).rank()==0 and M.rank()==121
c=[v.numerator*pow(v.denominator,-1,p)*pow(scale,-1,p)%p for v in coeffs]
aug=nmod_mat(M.tolist()+[c],p);rr,rank=aug.rref();assert rank==122
cols=pivots_from_rref(rr,rank)
rows=m.pivotcols(nmod_mat([[aug[i,j] for i in range(189)] for j in cols],p));assert len(rows)==122
det=int(nmod_mat([[aug[i,j] for j in cols] for i in rows],p).det());assert det

# A rational dual supported on 122 actual trace-word coefficient probes.
# Its full 188 exact pairings are checked below, independently of reconstruction.
probes=[words[j] for j in pc]
direct=nmod_mat([[row.get(w,0)%p for w in probes] for row in images]+
               [[var.get(w,0)*pow(scale,-1,p)%p for w in probes]],p)
rr,r=direct.rref();assert r==122
dcols=pivots_from_rref(rr,r)
drows=m.pivotcols(nmod_mat([[direct[i,j] for i in range(189)] for j in dcols],p))
assert len(drows)==122 and 188 in drows
minor=fmpq_mat([[images[i].get(probes[j],0) if i<188 else str(Fraction(var.get(probes[j],0),scale)) for j in dcols] for i in drows])
dual=minor.solve(fmpq_mat([[int(i==188)] for i in drows]))
dual=[Fraction(str(dual[i,0])) for i in range(122)]
dual_words=[probes[j] for j in dcols]
pairings=[sum(v*row.get(w,0) for w,v in zip(dual_words,dual)) for row in images]
pairing=sum(v*var.get(w,0) for w,v in zip(dual_words,dual))/scale
assert all(v==0 for v in pairings) and pairing==1

certificate={
 'status':'new_rational_fortuitous_class_certified',
 'group':'Sp(5)','defining_dimension':10,'native_weight':[3,3,9,0,0],'length':13,'j':36,'fermion_degree':4,
 'coefficient_convention':'Unchanged repository divided-power Taylor letters and formal highest-weight basis. Exported coefficients are the reconstructed rational coefficients.',
 'stable_character_dimension':37085,'stable_adjacent_ranks':[22495,14590],'stable_rational_cohomology':0,
 'complete_source_relation_dimension':188,'exact_rational_source_relation_image_rank':121,
 'preceding_relation_dimension_and_rank':[67,67],
 'target_identity_subspace_dimension':181,'target_identity_complete_dimension':183,
 'finite_rational_cohomology_lower_bound':1,'cycle_terms':len(poly),'variation_terms':len(var),
 'cycle_coefficient_common_denominator':str(scale),
 'formal_highest_weight_reconstruction_common_denominator':str(scale),
 'least_trace_coefficient_denominator':str(math.lcm(*(Fraction(v,scale).denominator for v in poly.values()))),
 'cycle_coefficient_common_denominator_note':'The reconstruction denominator clears the formal-HW coefficients; the unchanged trace polynomial can have a smaller least denominator.',
 'prime':p,
 'identity_basis_minor_word_indices':pc,'identity_basis_minor_determinant':detR,
 'nonboundary_minor_rows':rows,'nonboundary_minor_columns':cols,'nonboundary_minor_determinant':det,
 'nonboundary_augmented_rank':122,'dual_field':'Q','dual_boundary_pairings':[str(v) for v in pairings],
 'dual_cycle_defect_pairing':str(pairing),'dual_nonzero_coefficients':sum(bool(v) for v in dual),
 'dual_word_coefficients':[[str(v),m.encode(w)] for w,v in zip(dual_words,dual) if v],
 'defect_identity_coefficients':[str(v/scale) for v in coeffs],
 'exact_closure_checked':True,'native_highest_weights_checked':True,
 'finite_identity_proof':'Each seed is a Pfaffian antisymmetrization of twelve fundamental indices in dimension ten, with the remaining indices paired. Highest-weight projections, Taylor-profile components, and Q0 images preserve exact identities. The source identity basis saturates the exact character deficit.',
 'dependencies':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in [
  'pfaffian_relations.py','paired_relation_basis.py','path_relation_completion.py',
  'exact_relations_sp5_L12.pkl','exact_relations_sp5_L13.pkl','target_sp5_snapshot_181.pkl',
  'exact_q_relations_sp5_L12.pkl','exact_q_relations_sp5_L13.pkl','sp5_preceding_relation_rank.json',
  'sp5_relation_quotient.py','transport_target_recipes.py']},
 'seconds':time.monotonic()-start}
(HERE/'sp5_discovery_certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
operator={'group':'Sp(5)','native_weight':[3,3,9,0,0],'length':13,
          'terms':[[str(Fraction(v,scale)),m.encode(w)] for w,v in sorted(poly.items())]}
(HERE/'sp5_new_fortuitous_operator.json.gz').write_bytes(gzip.compress((json.dumps(operator,separators=(',',':'))+'\n').encode(),mtime=0))
print({k:certificate[k] for k in ['status','cycle_terms','variation_terms','nonboundary_minor_determinant','dual_cycle_defect_pairing','dual_nonzero_coefficients','seconds']},flush=True)
