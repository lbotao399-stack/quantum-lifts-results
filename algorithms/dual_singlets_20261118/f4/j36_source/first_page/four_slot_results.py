"""Project native four-column probes using the published F4 interpolation minor."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np

from four_slot_obstruction import as_array, quotient_map, read_minor


def prepare(folder,out):
    data={}
    for p in (1009,1013):
        info,m=read_minor(folder,p)
        q=quotient_map(info,p)
        weights=m.solve(q)
        assert m*weights == q
        data[str(p)] = as_array(weights).astype('<u2')
    np.savez_compressed(out,**data)


def report(folder,dual,p,candidate):
    pieces = sorted(folder.glob(f'adj-{p}-*.npy'),key=lambda path:int(path.stem.split('-')[-1]))
    if pieces:
        parts=[]
        covered=0
        for file in pieces:
            assert int(file.stem.split('-')[-1])==covered
            rows=np.load(file,allow_pickle=False)
            covered+=len(rows)
            parts.append(rows)
        assert covered==703
        data=np.concatenate(parts)
    else:
        data=np.load(folder/f'adj-{p}.npy',allow_pickle=False)
    assert data.shape==(703,750) and data.dtype==np.dtype('uint16') and data.max()<p
    forms=data.T.astype(np.int64)
    with np.load(dual,allow_pickle=False) as saved:
        basis=saved[str(p)].astype(np.int64)
    assert basis.shape==(703,7)
    reduced=forms@basis%p
    _,rank=nmod_mat(reduced.tolist(),p).rref()
    matrix=nmod_mat(forms.tolist(),p)
    rref,rank0=matrix.transpose().rref()
    picked=[next(j for j in range(rref.ncols()) if rref[i,j]) for i in range(rank0)]
    test=nmod_mat(forms[picked,:].tolist(),p)
    rr,r=test.rref()
    cols=[next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(r)]
    minor=forms[np.ix_(picked,cols)]
    det=int(nmod_mat(minor.tolist(),p).det())
    assert det and r==rank0
    out=dict(prime=p, tested=len(forms), independent=rank0, selected=picked, nodes=cols,
             determinant=det, residual_rank=rank, inactive_lower_bound=696+rank,
             residual_at_selected=reduced[picked].tolist(),
             little_endian_node_table_sha256=sha256(data.astype('<u2',copy=False).tobytes()).hexdigest())
    if candidate:
        proposal=json.loads(candidate.read_bytes())
        coefficients=np.array(proposal['weights'],dtype=np.int64)
        assert proposal['prime']==p and len(coefficients)>=703 and not np.any(coefficients[703:])
        out['nonzero_candidate_pairings']=[[i,int(v)] for i,v in enumerate(forms@coefficients[:703]%p) if v]
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    a=parser.add_subparsers(dest='mode',required=True)
    make=a.add_parser('prepare');make.add_argument('minor_folder',type=Path);make.add_argument('weights',type=Path)
    read=a.add_parser('report');read.add_argument('folder',type=Path);read.add_argument('weights',type=Path)
    read.add_argument('--prime',type=int,choices=(1009,1013),required=True);read.add_argument('--candidate',type=Path)
    read.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__: parser.error('enable assertions')
    if args.mode=='prepare':prepare(args.minor_folder,args.weights)
    else:
        result=report(args.folder,args.weights,args.prime,args.candidate)
        args.out.write_text(json.dumps(result,indent=2)+'\n')
        print(result['independent'],result['residual_rank'],result['determinant'])
