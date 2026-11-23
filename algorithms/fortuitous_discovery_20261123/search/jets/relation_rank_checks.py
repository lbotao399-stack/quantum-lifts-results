from relation_quotient import *

def certify(rows,p):
 words=sorted({w for row in rows for w in row});A=modmatrix(rows,words,p);rr,rank=A.rref();pc=pivots_from_rref(rr,rank)
 cols=nmod_mat([[A[i,j] for i in range(len(rows))] for j in pc],p);cr,rrank=cols.rref();assert rrank==rank
 pr=pivots_from_rref(cr,rank);minor=nmod_mat([[A[i,j] for j in pc] for i in pr],p);det=int(minor.det());assert det
 return {'rank':rank,'rows':len(rows),'word_columns':len(words),'minor_rows':pr,'minor_columns':pc,'minor_determinant':det}

start=time.monotonic();result={'group':'Sp4','weight':[3,3,7,0,0],'prime':1000003,'maps':{}}
for L in (10,11):
 rows,_=pickle.loads((HERE/f'exact_relations_sp4_L{L}.pkl').read_bytes());images=[m.q(r) for r in rows]
 assert all(not m.q(r) for r in images)
 (HERE/f'exact_q_relations_sp4_L{L}.pkl').write_bytes(pickle.dumps(images,protocol=4))
 result['maps'][L]=certify(images,1000003)
 print('Q rank',L,result['maps'][L]['rank'],'seconds',time.monotonic()-start,flush=True)
assert result['maps'][10]['rank']==67
assert result['maps'][11]['rank']==121
assert len(pickle.loads((HERE/'exact_relations_sp4_L11.pkl').read_bytes())[0])-67==121
result.update(exact_Q_squared_zero=True,rational_Q_I10_rank=67,rational_Q_I11_rank=121,seconds=time.monotonic()-start)
(HERE/'relation_rank_checks.json').write_text(json.dumps(result,indent=2))
