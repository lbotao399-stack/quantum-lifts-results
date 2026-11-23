"""Exact finite-rank bounds and restriction of the archived omega cycle."""
from pathlib import Path
import json,time,math,hashlib
from rank_runtime import HERE,POSITIVE,m,np,nmod_mat,load_stable_rows,load_omega
from dual_charge_bands import coordinates
from trace_descent import dimensions

def minor(values):
    a=nmod_mat(values.tolist(),m.P)
    cols=m.pivotcols(a);rows=m.pivotcols(a.transpose())
    assert len(cols)==len(rows)
    sq=values[np.ix_(rows,cols)].astype(np.uint32)
    det=int(nmod_mat(sq.tolist(),m.P).det()) if rows else 1
    assert det
    return {'rank':len(rows),'rows':rows,'columns':cols,'determinant':det},sq

def samples(dim,degree,salt):
    points=math.ceil((dim+16)/(degree+1))
    return [(min(8,points-i),salt+6151*i,list(range(min(8,points-i)*(degree+1))))
            for i in range(0,points,8)]

def main():
    started=time.monotonic();source=POSITIVE
    parents,centre,incoming,outgoing=load_stable_rows()
    stable=json.loads((source/'stable_rank_result.json').read_text())
    omega=load_omega()
    delta=m.q(omega);assert len(delta)==4502
    incoming=[incoming[i] for i in stable['incoming']['basis_rows']]
    outgoing=[outgoing[i] for i in stable['outgoing']['basis_rows']]
    result={'source_commit':'9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9',
            'weight':[3,3,5,0,0],'length':9,'prime':m.P,
            'cycle_terms':len(omega),'stable_derivative_terms':len(delta),'groups':[]}
    for size in (2,4):
        dims={y:sum(dimensions((3,3,5,0,0),y,size,'Sp').values()) for y in (8,9,10)}
        lib=m.library(HERE,size,m.P)
        samples_in=samples(dims[9],4,711233+size)
        letters_in=sorted({l for row in incoming+[omega] for w in row for tr in w for l in tr})
        before=coordinates(incoming+[omega],4,letters_in,samples_in,size,'Sp',lib,2)
        mi,sql=minor(before[:-1]);aug,sqlaug=minor(before)
        samples_out=samples(dims[10],5,912191+size)
        letters_out=sorted({l for row in outgoing for w in row for tr in w for l in tr})
        after=coordinates(outgoing,5,letters_out,samples_out,size,'Sp',lib,2)
        mo,sqr=minor(after)
        upper=dims[9]-mi['rank']-mo['rank']
        lower=aug['rank']-mi['rank']
        row={'defining_dimension':size,'character_dimensions':dims,
             'incoming':mi,'outgoing':mo,'incoming_plus_omega':aug,
             'cohomology_upper_bound':upper,'omega_nonexact_lower_bound':lower,
             'omega_evaluation_nonzeros':int(np.count_nonzero(before[-1])),
             'samples_incoming':samples_in,'samples_outgoing':samples_out,
             'incoming_alphabet':letters_in,'outgoing_alphabet':letters_out,
             'certified_H_zero':upper==0}
        assert 0<=lower<=upper
        np.savez_compressed(HERE/f'sp{size//2}_rank_minors.npz',incoming=sql,
                            outgoing=sqr,augmented=sqlaug,
                            omega_evaluation=before[-1])
        row['minors_sha256']=hashlib.sha256((HERE/f'sp{size//2}_rank_minors.npz').read_bytes()).hexdigest()
        result['groups'].append(row)
        print('group',size,'dimensions',dims,'ranks',mi['rank'],mo['rank'],
              'H_upper',upper,'omega_lower',lower,'seconds',time.monotonic()-started,flush=True)
        result['seconds']=time.monotonic()-started
        (HERE/'lower_rank_scan.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
