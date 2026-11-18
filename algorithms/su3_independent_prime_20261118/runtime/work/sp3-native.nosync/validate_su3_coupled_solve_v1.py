"""Synthetic coupled system: a zero RHS block, four free parameters, and an inconsistent control."""
import sp_environment
from pathlib import Path
import hashlib,json,pickle,tempfile,time
import numpy as np
from flint import nmod_mat
from su3_coupled_solve_v1 import PRIME,run,homogeneous

HERE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    started=time.monotonic();rng=np.random.default_rng(2026111831)
    with tempfile.TemporaryDirectory(prefix='su3-coupled-control-',dir=HERE) as name:
        root=Path(name)
        sources=[dict(index=i,dimension=n,targets=[0,1],profile=['synthetic',i]) for i,n in enumerate((3,2,4))]
        plan=dict(source_dimension=9,sources=sources,targets=[dict(index=0,dimension=2),dict(index=1,dimension=3)])
        write(root/'plan.json',plan)
        right=rng.integers(0,PRIME,(5,9),dtype=np.int64);right[:3]=0;right[:3,:3]=np.eye(3,dtype=np.int64)
        left=rng.integers(0,PRIME,(15,5),dtype=np.int64);left[:7,3:]=0;left[-1]=0
        matrix=left@right%PRIME
        assert int(nmod_mat(matrix.tolist(),PRIME).rank())==5
        known=rng.integers(0,PRIME,9,dtype=np.int64);known[:3]=0;rhs=matrix@known%PRIME;assert not rhs[:7].any() and rhs[7:].any()
        frames={0:dict(columns=7,prime=PRIME,synthetic=True),1:dict(columns=8,prime=PRIME,synthetic=True)}
        offset=0
        for source in sources:
            i,n=source['index'],source['dimension'];basis=root/'basis'/('profile%02d'%i);basis.mkdir(parents=True)
            seeds=[('synthetic',i,j) for j in range(n)];data=pickle.dumps(seeds,protocol=4)
            (basis/'independent_seeds.pkl').write_bytes(data)
            write(basis/'result.json',dict(status='complete_rational_profile_basis_recipes',rank_lower_bound=n,dimension=n,seeds_sha256=sha(basis/'independent_seeds.pkl'),synthetic_control_only=True))
            folder=root/'maps'/('profile%02d'%i)/('00000_%05d'%n);folder.mkdir(parents=True);(folder/'source_recipes.pkl').write_bytes(data)
            write(folder/'manifest.json',dict(source=source,prime=PRIME,targets=frames,source_recipes_sha256=sha(folder/'source_recipes.pkl'),start_row=0,stop_row=n))
            for target,begin,end in ((0,0,7),(1,7,15)):
                path=folder/('t%02d_00000_%05d.npy'%(target,n));np.save(path,matrix[begin:end,offset:offset+n].T.astype(np.uint32))
                write(path.with_suffix('.json'),dict(source=i,target=target,start=0,stop=n,columns=end-begin,sha256=sha(path)))
            offset+=n
        for target,begin,end in ((0,0,7),(1,7,15)):
            folder=root/'rhs'/('profile%02d'%target);folder.mkdir(parents=True);np.save(folder/'rhs.npy',rhs[begin:end].astype(np.uint32))
            write(folder/'manifest.json',dict(frame=frames[target]))
            write(folder/'result.json',dict(status='complete_negative_Q1_evaluations',columns=end-begin,sha256=sha(folder/'rhs.npy')))
        output=root/'solve';run(root,output,root/'plan.json')
        result=json.loads((output/'result.json').read_text())
        assert result['equation_rank']==5 and result['homogeneous_dimension']==4
        solution=np.load(output/'particular.npy').astype(np.int64);free=np.load(output/'free_variables.npy').tolist()
        assert not ((matrix@solution-rhs)%PRIME).any()
        stages=[]
        for folder in sorted(output.glob('stage*')):
            record=json.loads((folder/'result.json').read_text())
            stages.append(dict(rows=np.load(folder/'rows.npy'),columns=np.load(folder/'columns.npy'),pivots=record['pivots']))
        kernel=homogeneous(stages,9,free)
        assert not (matrix@kernel.T%PRIME).any() and np.array_equal(kernel[:,free],np.eye(4,dtype=np.int64))
        folder=root/'rhs/profile01';bad=rhs[7:].copy();bad[-1]=1;np.save(folder/'rhs.npy',bad.astype(np.uint32))
        write(folder/'result.json',dict(status='complete_negative_Q1_evaluations',columns=8,sha256=sha(folder/'rhs.npy')))
        run(root,root/'inconsistent',root/'plan.json')
        assert json.loads((root/'inconsistent/result.json').read_text())['status']=='inconsistent_modular_evaluation_equations'
    report=dict(status='passed_coupled_offsets_sign_zero_rhs_and_all_homogeneous_parameters',equations=15,variables=9,rank=5,
                full_kernel_dimension_checked=4,inconsistent_control=True,seconds=time.monotonic()-started,
                source_sha256={n:sha(HERE/n) for n in ('su3_coupled_solve_v1.py','validate_su3_coupled_solve_v1.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
