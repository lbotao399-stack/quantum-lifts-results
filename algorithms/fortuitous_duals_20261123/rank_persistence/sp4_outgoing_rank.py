"""Fresh Sp(4) outgoing minor in the archived Sp(3) omega grading.

A nonzero1198-minor equals the complete stable outgoing rank. Combined
with source-degree trace independence and stable exactness, this excludes
every Sp(4) cohomology class in this grading, not only a chosen lift.
"""
from pathlib import Path
import json,time,hashlib,math
from rank_runtime import HERE,POSITIVE,m,np,nmod_mat,load_stable_rows
from dual_charge_bands import coordinates
from trace_descent import dimensions

def main():
    started=time.monotonic()
    source=POSITIVE
    parents,centre,incoming,outgoing=load_stable_rows()
    stable=json.loads((source/'stable_rank_result.json').read_text())
    keep=stable['outgoing']['basis_rows']
    rows=[outgoing[i] for i in keep]
    assert len(rows)==1198 and all(m.nferm(w)==5 for row in rows for w in row)
    letters=sorted({l for row in rows for w in row for tr in w for l in tr})
    lib=m.library(HERE,8,m.P)
    chunks=[];samples=[];rank=0
    # Each independent dense exterior sample supplies six exact coordinates.
    points=math.ceil((len(rows)+18)/6)
    initial=[(min(8,points-i),812311+6151*i,list(range(min(8,points-i)*6)))
             for i in range(0,points,8)]
    for sample in initial:
        values=coordinates(rows,5,letters,[sample],8,'Sp',lib,2)
        chunks.append(values);samples.append(sample)
        print('samples',len(samples),'columns',sum(x.shape[1] for x in chunks),
              'seconds',time.monotonic()-started,flush=True)
    for extra in range(5):
        matrix=np.column_stack(chunks)
        field=nmod_mat(matrix.tolist(),m.P)
        columns=m.pivotcols(field)
        rank=len(columns)
        print('rank',rank,'formal_upper',len(rows),'seconds',time.monotonic()-started,flush=True)
        if rank==len(rows):break
        sample=(8,975221+6151*extra,list(range(48)))
        chunks.append(coordinates(rows,5,letters,[sample],8,'Sp',lib,2));samples.append(sample)
    assert rank==1198,('rank did not saturate',rank)
    minor=matrix[:,columns].astype(np.uint32)
    determinant=int(nmod_mat(minor.tolist(),m.P).det());assert determinant
    np.savez_compressed(HERE/'sp4_outgoing_minor.npz',minor=minor,
                        columns=np.asarray(columns),source_rows=np.asarray(keep))
    exact_dims={str(y):sum(dimensions((3,3,5,0,0),y,8,'Sp').values()) for y in (8,9,10)}
    assert exact_dims['8']==len(parents)==2624 and exact_dims['9']==len(centre)==2870
    result={'passed':True,'source_commit':'9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9',
            'runtime_commit':'082be7070d207e2cdecfd0a91b03d5e014fd303d',
            'weight':[3,3,5,0,0],'length':9,'group':'Sp','defining_dimension':8,
            'prime':m.P,'finite_character_dimensions':exact_dims,
            'stable_incoming_rank':1672,'stable_outgoing_rank':1198,
            'outgoing_minor_size':1198,'outgoing_minor_determinant':determinant,
            'samples':samples,'alphabet':letters,'selected_stable_outgoing_rows':keep,
            'source_rows_sha256':hashlib.sha256((source/'stable_rows.pkl').read_bytes()).hexdigest(),
            'minor_sha256':hashlib.sha256((HERE/'sp4_outgoing_minor.npz').read_bytes()).hexdigest(),
            'cohomology_upper_bound':2870-1672-1198,'seconds':time.monotonic()-started,
            'rational_bridge':'The integral stable incoming minor gives rank1672. Sp4 specialization is injective at length9 because its character dimension equals the complete stable trace dimension2870. The freshly evaluated nonzero modular outgoing minor gives rational rank>=1198, equal to the stable upper bound. Therefore rational H9 is zero.',
            'scope':'Complete fixed highest-weight grading, conditional on the archived character formula and trace-generation theorem. No statement about all other gradings or about quantum corrections.'}
    (HERE/'sp4_outgoing_rank.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['samples','alphabet','selected_stable_outgoing_rows']},indent=2),flush=True)

if __name__=='__main__':main()
