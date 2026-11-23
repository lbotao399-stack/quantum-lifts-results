from relation_quotient import *
import gzip
start=time.monotonic();p=2147483647;work=HERE/'quotient_sp4_target_snapshot_181_p2147483647'
cycle=pickle.loads((work/'rational_cycle.pkl').read_bytes());poly=cycle['polynomial'];scale=cycle['original_scale_denominator'];coeffs=cycle['relation_coefficients']
relations,recipes=pickle.loads((HERE/'target_snapshot_181.pkl').read_bytes());images=pickle.loads((HERE/'exact_q_relations_sp4_L11.pkl').read_bytes())
words=sorted({w for row in relations for w in row}|{w for row in images for w in row});R=modmatrix(relations,words,p);rr,rank=R.rref();assert rank==181
pc=pivots_from_rref(rr,rank);T=nmod_mat([[R[i,j] for j in pc] for i in range(181)],p);detR=int(T.det());assert detR
Q=modmatrix(images,words,p);Qp=nmod_mat([[Q[i,j] for j in pc] for i in range(188)],p)
M=T.transpose().solve(Qp.transpose()).transpose();assert (Q-M*R).rank()==0 and M.rank()==121
c=[v.numerator*pow(v.denominator,-1,p)*pow(scale,-1,p)%p for v in coeffs]
aug=nmod_mat(M.tolist()+[c],p);rr,rank=aug.rref();assert rank==122
cols=pivots_from_rref(rr,rank);rows=m.pivotcols(nmod_mat([[aug[i,j] for i in range(189)] for j in cols],p));assert len(rows)==122
minor=nmod_mat([[aug[i,j] for j in cols] for i in rows],p);det=int(minor.det());assert det
system=nmod_mat([row+[0] for row in M.tolist()]+[c+[1]],p);reduced,rank=system.rref();assert rank==122
pivot=pivots_from_rref(reduced,rank);assert all(i<181 for i in pivot)
dual=[0]*181
for i,j in enumerate(pivot):dual[j]=int(reduced[i,181])
assert (M*nmod_mat([[x] for x in dual],p)).rank()==0
assert sum(x*y for x,y in zip(c,dual))%p==1
report=json.loads((work/'result.json').read_bytes());var=m.q(poly)
assert len(poly)==17430 and len(var)==15556
certificate={'status':'new_rational_fortuitous_class_certified','group':'Sp(4)','defining_dimension':8,'native_weight':[3,3,7,0,0],'length':11,'j':32,'fermion_degree':4,'coefficient_convention':'Unchanged repository divided-power Taylor letters and formal highest-weight basis. The coefficients below are the reconstructed rational coefficients, not a rescaled representative.','stable_character_dimension':11627,'stable_adjacent_ranks':[6940,4687],'stable_rational_cohomology':0,'complete_source_relation_dimension':188,'exact_rational_source_relation_image_rank':121,'preceding_relation_dimension_and_rank':[67,67],'target_identity_subspace_dimension':181,'target_identity_complete_dimension':183,'finite_rational_cohomology_lower_bound':1,'cycle_terms':len(poly),'variation_terms':len(var),'cycle_coefficient_common_denominator':str(scale),'prime':p,'identity_basis_minor_word_indices':pc,'identity_basis_minor_determinant':detR,'nonboundary_minor_rows':rows,'nonboundary_minor_columns':cols,'nonboundary_minor_determinant':det,'nonboundary_augmented_rank':122,'dual_on_target_identity_basis':dual,'dual_boundary_pairings':0,'dual_cycle_defect_pairing':1,'defect_identity_coefficients':[str(v/scale) for v in coeffs],'exact_closure_checked':True,'finite_identity_proof':'Each seed is a Pfaffian antisymmetrization of ten fundamental indices in dimension eight, with the remaining indices paired. Highest-weight projections, Taylor-profile components, and Q0 images preserve the exact identity. The source identity basis saturates the exact character deficit.','dependencies':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['pfaffian_relations.py','paired_relation_basis.py','path_relation_completion.py','exact_relations_sp4_L10.pkl','exact_relations_sp4_L11.pkl','target_snapshot_181.pkl','relation_rank_checks.json','formal_frontier.json']},'seconds':time.monotonic()-start}
(HERE/'sp4_discovery_certificate.json').write_text(json.dumps(certificate,indent=2))
operator={'group':'Sp(4)','native_weight':[3,3,7,0,0],'length':11,'terms':[[str(Fraction(v,scale)),m.encode(w)] for w,v in sorted(poly.items())]}
(HERE/'sp4_new_fortuitous_operator.json.gz').write_bytes(gzip.compress((json.dumps(operator,separators=(',',':'))+'\n').encode(),mtime=0))
print({k:certificate[k] for k in ['status','cycle_terms','variation_terms','nonboundary_minor_determinant','dual_cycle_defect_pairing','seconds']},flush=True)
