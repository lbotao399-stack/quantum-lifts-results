"""Couple both correction equations without fixing the first homogeneous freedom.

The first staged system keeps its pivot indices. Only its augmented RHS column
moves past the newly appended second-correction variables. Modular consistency
or inconsistency is not, by itself, a rational obstruction certificate.
"""
import sp_environment
from pathlib import Path
import argparse,json,pickle,time
import numpy as np
from flint import ctx
from bounded_flint_rref import reduce_target
from staged_sp2_first_correction import back_substitute
from exact_blas_rows import matmul_mod
from su3_coupled_solve_v1 import PRIME,HERE,sha,save,save_array,offsets,assemble as assemble_first


def lift_first_stages(folder,first_plan_path,variables):
    result=json.loads((folder/'result.json').read_text())
    assert result['status']=='complete_modular_evaluation_solution_with_all_free_parameters'
    plan=json.loads(first_plan_path.read_text());_,first_variables=offsets(plan)
    assert result['prime']==PRIME and result['variables']==first_variables<variables
    assert sha(folder/'particular.npy')==result['particular_sha256']
    assert sha(folder/'free_variables.npy')==result['free_variables_sha256']
    stages=[];pivots=set()
    for number,target in enumerate(result['processed_targets']):
        part=folder/('stage%02d_target%02d'%(number,target));marker=part/'result.json'
        digest=sha(marker);record=json.loads(marker.read_text());manifest=record['manifest']
        assert digest==result['stage_certificates'][number]
        assert manifest['previous_stages']==result['stage_certificates'][:number]
        assert manifest['plan_sha256']==sha(first_plan_path)
        assert manifest['variables']==first_variables and manifest['prime']==PRIME
        assert record['consistent'] and not pivots.intersection(record['pivots'])
        assert sha(part/'rows.npy')==record['rows_sha256'] and sha(part/'columns.npy')==record['columns_sha256']
        rows=np.load(part/'rows.npy',mmap_mode='r');columns=np.load(part/'columns.npy').copy()
        assert rows.shape==(record['rank'],len(columns)) and columns[-1]==first_variables
        assert np.all(columns[1:]>columns[:-1]) and all(0<=p<first_variables for p in record['pivots'])
        columns[-1]=variables
        stages.append(dict(rows=rows,columns=columns,pivots=record['pivots'],rank=record['rank'],
            consistent=True,certificate_sha256=digest,first_equation=True))
        pivots.update(record['pivots'])
    assert len(stages)==len(plan['targets'])==len(result['stage_certificates'])
    assert len(pivots)==result['equation_rank']
    assert np.array_equal(np.load(folder/'free_variables.npy'),sorted(set(range(first_variables))-pivots))
    assert first_variables-len(pivots)==result['homogeneous_dimension']
    return stages,result


def basis_record(folder,source):
    record=json.loads((folder/'result.json').read_text())
    assert record['status']=='complete_rational_profile_basis_recipes'
    assert record['dimension']==record['rank_lower_bound']==source['dimension']
    assert sha(folder/'independent_seeds.pkl')==record['seeds_sha256']
    seeds=pickle.loads((folder/'independent_seeds.pkl').read_bytes())
    assert len(seeds)==source['dimension']
    return record,seeds


def locate_second_bases(roots,plan):
    found={}
    for role in ('sources','targets'):
        for source in plan[role]:
            candidates=[]
            for root in roots:
                for path in sorted(root.glob('shard*/%s%02d/result.json'%(role,source['index']))):
                    record=json.loads(path.read_text())
                    if record['status']=='complete_rational_profile_basis_recipes':candidates.append(path.parent)
            if not candidates:raise FileNotFoundError('Incomplete second basis: %s %d'%(role,source['index']))
            folder=candidates[-1];record,_=basis_record(folder,source)
            found[(role,source['index'])]=dict(path=folder,seeds_sha256=record['seeds_sha256'],dimension=source['dimension'])
    return found


def assemble_second(q0_root,q1_root,rhs_root,first_basis,second_bases,first_plan,quantum_plan,second_plan,target):
    first_starts,first_variables=offsets(first_plan)
    second_starts,second_variables=offsets(second_plan);variables=first_variables+second_variables
    index=target['index'];folder=rhs_root/'rhs'/('profile%02d'%index)
    rhs=json.loads((folder/'result.json').read_text());manifest=json.loads((folder/'manifest.json').read_text())
    assert rhs['status']=='complete_negative_Q2_evaluations' and rhs['target']==index and rhs['points_override']==0
    assert sha(folder/'rhs.npy')==rhs['sha256'] and rhs['source_sha256']==manifest['source_sha256']
    width=rhs['columns'];frame=manifest['frame'];assert width==frame['columns'] and frame['prime']==PRIME
    assert manifest['target']==target and manifest['points_override']==0
    by_index={x['index']:x for x in first_plan['sources']}
    assert len(by_index)==len(quantum_plan['sources'])
    for source in quantum_plan['sources']:
        assert source['dimension']==by_index[source['index']]['dimension']
        assert source['profile']==by_index[source['index']]['profile']
    groups=[]
    for kind,root,plan,starts,offset in (('Q1',q1_root,quantum_plan,first_starts,0),('Q0',q0_root,second_plan,second_starts,first_variables)):
        for source in plan['sources']:
            if index in source['targets']:groups.append((kind,root,source,starts[source['index']]+offset))
    columns=np.asarray([i for _,_,source,start in groups for i in range(start,start+source['dimension'])]+[variables],dtype=np.int64)
    assert np.all(columns[1:]>columns[:-1])
    values=np.zeros((width,len(columns)),dtype=np.int32);cursor=0;blocks=[]
    for kind,root,source,global_start in groups:
        s=source['index'];basis=first_basis/('profile%02d'%s) if kind=='Q1' else second_bases[('sources',s)]['path']
        source_record,seeds=basis_record(basis,source)
        covered=np.zeros(source['dimension'],dtype=bool)
        source_root=root/'maps'/('profile%02d'%s)
        if not source_root.exists():raise FileNotFoundError('Missing %s source profile %d'%(kind,s))
        for part in sorted(source_root.iterdir()):
            if not part.is_dir() or not (part/'manifest.json').exists():continue
            local=json.loads((part/'manifest.json').read_text())
            assert local['source']==source and local['prime']==PRIME
            assert local['targets'][str(index)]==frame
            if kind=='Q1':
                assert local['points_override']==0 and local['source_basis_sha256']==source_record['seeds_sha256']
            assert sha(part/'source_recipes.pkl')==local['source_recipes_sha256']
            assert pickle.loads((part/'source_recipes.pkl').read_bytes())==seeds[local['start_row']:local['stop_row']]
            for marker in sorted(part.glob('t%02d_*.json'%index)):
                record=json.loads(marker.read_text());start,stop=record['start'],record['stop']
                assert record['source']==s and record['target']==index and record['columns']==width
                assert 0<=local['start_row']<=start<stop<=local['stop_row']<=source['dimension']
                assert not covered[start:stop].any(),'Overlapping map rows'
                path=marker.with_suffix('.npy');assert sha(path)==record['sha256']
                matrix=np.load(path,mmap_mode='r')
                assert matrix.shape==(stop-start,width) and np.all((matrix>=0)&(matrix<PRIME))
                values[:,cursor+start:cursor+stop]=matrix.T;covered[start:stop]=True
                blocks.append(dict(kind=kind,source=s,start=start,stop=stop,global_start=global_start+start,
                    sha256=record['sha256'],source_basis_sha256=source_record['seeds_sha256']))
        if not covered.all():raise FileNotFoundError('Target %d %s source %d: %d missing rows'%(index,kind,s,int((~covered).sum())))
        cursor+=source['dimension']
    right=np.load(folder/'rhs.npy');assert right.shape==(width,) and np.all((right>=0)&(right<PRIME))
    values[:,-1]=right
    return values,columns,dict(target=index,variables=variables,first_variables=first_variables,second_variables=second_variables,
        prime=PRIME,columns=width,blocks=blocks,rhs_sha256=rhs['sha256'],rhs_source_sha256=rhs['source_sha256'],rhs_frame=frame)


def run(first_root,first_solve,q0_root,q1_root,rhs_root,second_basis_roots,output,seconds=20000,threads=1,
        first_plan_path=None,quantum_plan_path=None,second_plan_path=None):
    assert 1<=threads<=64
    ctx.threads=threads
    started=time.monotonic();output.mkdir(parents=True,exist_ok=True)
    first_plan_path=first_plan_path or HERE/'su3_correction_profiles_v1.json'
    quantum_plan_path=quantum_plan_path or HERE/'su3_second_quantum_profiles_v1.json'
    second_plan_path=second_plan_path or HERE/'su3_second_correction_profiles_v1.json'
    first_plan=json.loads(first_plan_path.read_text());quantum_plan=json.loads(quantum_plan_path.read_text());second_plan=json.loads(second_plan_path.read_text())
    _,first_variables=offsets(first_plan);_,second_variables=offsets(second_plan);variables=first_variables+second_variables
    validation=json.loads((HERE/'bounded_flint_rref_validation.json').read_text())
    assert validation['status']=='passed_bounded_flint_rref'
    for name,digest in validation['source_sha256'].items():assert sha(HERE/name)==digest
    stages,first_result=lift_first_stages(first_solve,first_plan_path,variables)
    inherited=len(stages);processed=[]
    try:second_bases=locate_second_bases(second_basis_roots,second_plan)
    except FileNotFoundError as error:
        save(output/'result.json',dict(status='waiting_for_complete_second_bases',error=str(error),correction_solved=False));return
    for target in sorted(second_plan['targets'],key=lambda x:(x['dimension'],x['index'])):
        part=output/('stage%02d_target%02d'%(len(processed),target['index']));part.mkdir(exist_ok=True)
        try:values,columns,manifest=assemble_second(q0_root,q1_root,rhs_root,first_root/'basis',second_bases,first_plan,quantum_plan,second_plan,target)
        except FileNotFoundError as error:
            save(output/'result.json',dict(status='waiting_for_complete_second_maps',processed_targets=processed,error=str(error),correction_solved=False));return
        manifest.update(previous_stages=[s['certificate_sha256'] for s in stages],driver_sha256=sha(__file__),
            first_plan_sha256=sha(first_plan_path),second_plan_sha256=sha(second_plan_path),quantum_plan_sha256=sha(quantum_plan_path),
            first_solve_sha256=sha(first_solve/'result.json'))
        marker=part/'result.json'
        if marker.exists():
            record=json.loads(marker.read_text());assert record['manifest']==manifest
            assert sha(part/'rows.npy')==record['rows_sha256'] and sha(part/'columns.npy')==record['columns_sha256']
            stage=dict(rows=np.load(part/'rows.npy',mmap_mode='r'),columns=np.load(part/'columns.npy'),
                pivots=record['pivots'],rank=record['rank'],consistent=record['consistent'])
        else:
            stage=reduce_target(values,columns,stages,variables,PRIME)
            save_array(part/'rows.npy',stage['rows']);save_array(part/'columns.npy',stage['columns'])
            record=dict(manifest=manifest,pivots=stage['pivots'],rank=stage['rank'],consistent=stage['consistent'],
                seconds=stage['seconds'],rows_sha256=sha(part/'rows.npy'),columns_sha256=sha(part/'columns.npy'))
            save(marker,record)
        del values
        stage['certificate_sha256']=sha(marker);stages.append(stage);processed.append(target['index'])
        save(output/'progress.json',dict(processed_targets=processed,rank=sum(s['rank'] for s in stages),
            variables=variables,first_homogeneous_dimension=first_result['homogeneous_dimension'],seconds=time.monotonic()-started))
        if not stage['consistent']:
            save(output/'result.json',dict(status='inconsistent_modular_second_correction_system',processed_targets=processed,
                correction_solved=False,rational_obstruction_certified=False,
                scope='Inconsistent modulo the pinned prime with every first-correction parameter retained. Rational obstruction certification remains separate.'));return
        if time.monotonic()-started>=seconds and len(processed)<len(second_plan['targets']):
            save(output/'result.json',dict(status='partial_second_coupled_elimination',processed_targets=processed,correction_solved=False));return
    solution,free=back_substitute(stages,variables,PRIME);residuals=[]
    for target in first_plan['targets']:
        values,columns,_=assemble_first(first_root,first_plan,target)
        residual=(matmul_mod(values[:,:-1].astype(np.int64),solution[columns[:-1],None],PRIME)[:,0]-values[:,-1])%PRIME
        assert not residual.any(),('first',target['index'])
        residuals.append(dict(equation=1,target=target['index'],coordinates=len(residual),nonzero=0))
    for target in second_plan['targets']:
        values,columns,_=assemble_second(q0_root,q1_root,rhs_root,first_root/'basis',second_bases,first_plan,quantum_plan,second_plan,target)
        residual=(matmul_mod(values[:,:-1].astype(np.int64),solution[columns[:-1],None],PRIME)[:,0]-values[:,-1])%PRIME
        assert not residual.any(),('second',target['index'])
        residuals.append(dict(equation=2,target=target['index'],coordinates=len(residual),nonzero=0))
    save_array(output/'particular.npy',solution.astype(np.uint32));save_array(output/'free_variables.npy',np.asarray(free,dtype=np.int64))
    result=dict(status='complete_modular_second_solution_with_all_first_parameters',prime=PRIME,
        variables=variables,first_variables=first_variables,second_variables=second_variables,
        inherited_first_stages=inherited,first_homogeneous_dimension=first_result['homogeneous_dimension'],
        equation_rank=sum(s['rank'] for s in stages),homogeneous_dimension=len(free),processed_targets=processed,
        residuals=residuals,particular_sha256=sha(output/'particular.npy'),free_variables_sha256=sha(output/'free_variables.npy'),
        stage_certificates=[s['certificate_sha256'] for s in stages],seconds=time.monotonic()-started,flint_threads=threads,
        rational_correction_certified=False,rational_obstruction_certified=False,
        scope='Both modular correction equations are solved with the complete first homogeneous freedom. Rational certification remains separate.')
    save(output/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    for name in ('first-root','first-solve','q0-root','q1-root','rhs-root','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--second-basis-roots',type=Path,nargs='+',required=True)
    parser.add_argument('--seconds',type=int,default=20000)
    parser.add_argument('--threads',type=int,default=1)
    for name in ('first-plan-path','quantum-plan-path','second-plan-path'):parser.add_argument('--'+name,type=Path)
    run(**vars(parser.parse_args()))
