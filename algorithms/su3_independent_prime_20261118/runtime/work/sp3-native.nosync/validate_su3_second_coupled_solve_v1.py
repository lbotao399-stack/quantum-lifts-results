"""Independent dense controls for coupled equations and inherited affine stages."""
import sp_environment
from pathlib import Path
import hashlib,json,pickle,tempfile,time
import numpy as np
from flint import nmod_mat
import su3_coupled_solve_v1 as first
import su3_second_coupled_solve_v1 as second

HERE=Path(__file__).resolve().parent;PRIME=first.PRIME
def write(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')
def rank(matrix):return int(nmod_mat(matrix.tolist(),PRIME).rank())

def basis(folder,item):
    folder.mkdir(parents=True,exist_ok=True)
    seeds=[('control',item['index'],i) for i in range(item['dimension'])]
    (folder/'independent_seeds.pkl').write_bytes(pickle.dumps(seeds,protocol=4))
    write(folder/'result.json',dict(status='complete_rational_profile_basis_recipes',dimension=len(seeds),
        rank_lower_bound=len(seeds),seeds_sha256=first.sha(folder/'independent_seeds.pkl'),synthetic_control_only=True))
    return seeds

def maps(root,sources,frames,matrix,basis_paths,kind):
    offset=0
    for item in sources:
        index,n=item['index'],item['dimension'];record,seeds=second.basis_record(basis_paths[index],item)
        intervals=[(0,1),(1,n)] if n>1 else [(0,1)]
        for start,stop in intervals:
            folder=root/'maps'/('profile%02d'%index)/('rows%05d_%05d'%(start,stop));folder.mkdir(parents=True)
            (folder/'source_recipes.pkl').write_bytes(pickle.dumps(seeds[start:stop],protocol=4))
            manifest=dict(source=item,prime=PRIME,targets=frames,start_row=start,stop_row=stop,
                source_recipes_sha256=first.sha(folder/'source_recipes.pkl'))
            if kind=='Q1':manifest.update(points_override=0,source_basis_sha256=record['seeds_sha256'])
            write(folder/'manifest.json',manifest)
            begin=0
            for target,frame in frames.items():
                width=frame['columns'];path=folder/('t%02d_%05d_%05d.npy'%(target,start,stop))
                np.save(path,matrix[begin:begin+width,offset+start:offset+stop].T.astype(np.uint32))
                write(path.with_suffix('.json'),dict(source=index,target=target,start=start,stop=stop,columns=width,sha256=first.sha(path)))
                begin+=width
        offset+=n

def right_hand(root,targets,frames,rhs,order):
    begin=0
    for target in targets:
        index=target['index'];frame=frames[index];width=frame['columns'];folder=root/'rhs'/('profile%02d'%index);folder.mkdir(parents=True,exist_ok=True)
        np.save(folder/'rhs.npy',rhs[begin:begin+width].astype(np.uint32))
        manifest=dict(frame=frame);record=dict(status='complete_negative_Q%d_evaluations'%order,columns=width,sha256=first.sha(folder/'rhs.npy'))
        if order==2:
            manifest.update(target=target,points_override=0,source_sha256='e'*64)
            record.update(target=index,points_override=0,source_sha256='e'*64)
        write(folder/'manifest.json',manifest);write(folder/'result.json',record);begin+=width

def main():
    started=time.monotonic();checks=[]
    with tempfile.TemporaryDirectory(prefix='second-coupled-control-',dir=HERE) as temporary:
        root=Path(temporary);first_root=root/'first';q0=root/'q0';q1=root/'q1';rhs_root=root/'second_rhs';basis_root=root/'second_basis'
        source1=[dict(index=i,dimension=n,profile=['synthetic_first',i],targets=[0,1]) for i,n in enumerate((3,2,4))]
        source2=[dict(index=i,dimension=n,profile=['synthetic_second',i],targets=[0,1]) for i,n in enumerate((3,4))]
        targets=[dict(index=0,dimension=2),dict(index=1,dimension=3)]
        plan1=dict(source_dimension=9,sources=source1,targets=targets)
        plan2=dict(source_dimension=7,sources=source2,targets=targets)
        quantum=dict(sources=source1,targets=targets)
        for name,data in (('first',plan1),('second',plan2),('quantum',quantum)):write(root/(name+'.json'),data)
        frames={i:dict(columns=n,prime=PRIME,synthetic_control_only=True) for i,n in enumerate((7,8))}
        matrix1=np.zeros((15,9),dtype=np.int64);matrix1[:5,:5]=np.eye(5,dtype=np.int64)
        matrix1[5,:5]=[1,2,3,4,5];matrix1[6,:5]=[2,0,1,0,3];matrix1[7:12,:5]=3*np.eye(5,dtype=np.int64)
        known=np.arange(3,19,dtype=np.int64);right1=matrix1@known[:9]%PRIME
        first_paths={}
        for item in source1:
            path=first_root/'basis'/('profile%02d'%item['index']);basis(path,item);first_paths[item['index']]=path
        maps(first_root,source1,frames,matrix1,first_paths,'Q0');right_hand(first_root,targets,frames,right1,1)
        first.run(first_root,root/'first_solve',root/'first.json')
        second_paths={}
        for item in source2:
            path=basis_root/'shard0'/('sources%02d'%item['index']);basis(path,item);second_paths[item['index']]=path
        for item in targets:basis(basis_root/'shard0'/('targets%02d'%item['index']),item)
        matrix2=np.zeros((15,16),dtype=np.int64)
        matrix2[0,5]=1;matrix2[1,[6,9]]=1;matrix2[2,[7,10]]=[1,2]
        matrix2[3,[8,11]]=1;matrix2[4,[12,13]]=1;matrix2[5,[0,14,15]]=1
        matrix2[6]=matrix2[0]+2*matrix2[4];matrix2[7:13]=3*matrix2[:6];matrix2[13]=matrix2[2]+matrix2[5]
        right2=matrix2@known%PRIME
        old_part=np.load(root/'first_solve/particular.npy').astype(np.int64)
        fixed_rhs=(right2-matrix2[:,:9]@old_part)%PRIME
        assert rank(np.column_stack((matrix2[:,9:],fixed_rhs)))==rank(matrix2[:,9:])+1
        checks.append('fixing_the_first_particular_solution_gives_a_false_obstruction')
        maps(q1,source1,frames,matrix2[:,:9],first_paths,'Q1');maps(q0,source2,frames,matrix2[:,9:],second_paths,'Q0')
        right_hand(rhs_root,targets,frames,right2,2)
        kwargs=dict(first_root=first_root,first_solve=root/'first_solve',q0_root=q0,q1_root=q1,rhs_root=rhs_root,
            second_basis_roots=[basis_root],first_plan_path=root/'first.json',quantum_plan_path=root/'quantum.json',second_plan_path=root/'second.json',threads=4)
        second.run(output=root/'solve',seconds=0,**kwargs)
        assert json.loads((root/'solve/result.json').read_text())['status']=='partial_second_coupled_elimination'
        second.run(output=root/'solve',**kwargs)
        result=json.loads((root/'solve/result.json').read_text());assert result['status']=='complete_modular_second_solution_with_all_first_parameters'
        matrix=np.vstack((np.column_stack((matrix1,np.zeros((15,7),dtype=np.int64))),matrix2));right=np.concatenate((right1,right2))
        solution=np.load(root/'solve/particular.npy').astype(np.int64);free=np.load(root/'solve/free_variables.npy').tolist()
        assert not ((matrix@solution-right)%PRIME).any()
        assert rank(matrix)==result['equation_rank']==11 and len(free)==result['homogeneous_dimension']==5
        stages,_=second.lift_first_stages(root/'first_solve',root/'first.json',16)
        for folder in sorted((root/'solve').glob('stage*')):
            record=json.loads((folder/'result.json').read_text())
            stages.append(dict(rows=np.load(folder/'rows.npy'),columns=np.load(folder/'columns.npy'),pivots=record['pivots']))
        kernel=first.homogeneous(stages,16,free)
        assert not (matrix@kernel.T%PRIME).any() and np.array_equal(kernel[:,free],np.eye(5,dtype=np.int64))
        checks.extend(['both_equation_residuals_and_all_five_homogeneous_parameters','nonzero_first_rhs_column_relocation','checkpoint_resume','nonzero_global_row_offsets'])
        second_bases=second.locate_second_bases([basis_root],plan2)
        assemble_args=(q0,q1,rhs_root,first_root/'basis',second_bases,plan1,quantum,plan2,targets[0])
        marker=q1/'maps/profile00/rows00001_00003/t00_00001_00003.json';saved=marker.read_bytes();marker.unlink()
        try:second.assemble_second(*assemble_args)
        except FileNotFoundError:checks.append('missing_rows_rejected')
        else:raise AssertionError('missing rows accepted')
        marker.write_bytes(saved)
        copied=marker.with_name('t00_duplicate.json');copied.write_bytes(saved);copied.with_suffix('.npy').write_bytes(marker.with_suffix('.npy').read_bytes())
        try:second.assemble_second(*assemble_args)
        except AssertionError:checks.append('overlapping_rows_rejected')
        else:raise AssertionError('overlapping rows accepted')
        copied.unlink();copied.with_suffix('.npy').unlink()
        marker=marker.parent/'manifest.json';saved=marker.read_bytes();changed=json.loads(saved);changed['targets']['0']['columns']+=1;write(marker,changed)
        try:second.assemble_second(*assemble_args)
        except AssertionError:checks.append('incompatible_evaluation_frames_rejected')
        else:raise AssertionError('mismatched frames accepted')
        marker.write_bytes(saved)
        bad=right2.copy();bad[-1]=1;right_hand(rhs_root,targets,frames,bad,2)
        assert rank(np.column_stack((matrix,np.concatenate((right1,bad)))))==rank(matrix)+1
        second.run(output=root/'inconsistent',**kwargs)
        assert json.loads((root/'inconsistent/result.json').read_text())['status']=='inconsistent_modular_second_correction_system'
        checks.append('genuine_inconsistency_agrees_with_independent_dense_rank')
    report=dict(status='passed_second_coupled_assembly_all_homogeneous_parameters_and_dense_controls',
        checks=checks,equations=30,variables=16,rank=11,homogeneous_dimension=5,
        source_sha256={name:first.sha(HERE/name) for name in ('su3_second_coupled_solve_v1.py','validate_su3_second_coupled_solve_v1.py','su3_coupled_solve_v1.py')},
        seconds=time.monotonic()-started)
    output=Path(__file__).with_suffix('.json');write(output,report);print(json.dumps(report),flush=True)

if __name__=='__main__':main()
