"""Solve the coupled evaluation equations, retaining the complete homogeneous freedom.

The staged reducer is the previously validated bounded FLINT implementation.
A modular solution is not a rational polynomial-identity certificate.
"""
import sp_environment
from pathlib import Path
import argparse,hashlib,json,pickle,time
import numpy as np
from bounded_flint_rref import reduce_target
from staged_sp2_first_correction import back_substitute
from exact_blas_rows import matmul_mod

HERE=Path(__file__).resolve().parent
PRIME=1000033


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,data):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)


def save_array(path,data):
    temp=path.with_suffix('.tmp.npy');np.save(temp,data);temp.replace(path)


def offsets(plan):
    result={};total=0
    for source in plan['sources']:
        result[source['index']]=total;total+=source['dimension']
    assert total==plan['source_dimension']
    return result,total


def assemble(root,plan,target):
    starts,variables=offsets(plan);index=target['index'];sources=[s for s in plan['sources'] if index in s['targets']]
    rhs_root=root/'rhs'/('profile%02d'%index)
    rhs_record=json.loads((rhs_root/'result.json').read_text());rhs_manifest=json.loads((rhs_root/'manifest.json').read_text())
    assert rhs_record['status']=='complete_negative_Q1_evaluations' and sha(rhs_root/'rhs.npy')==rhs_record['sha256']
    width=rhs_record['columns'];assert width==rhs_manifest['frame']['columns']
    columns=np.asarray([i for source in sources for i in range(starts[source['index']],starts[source['index']]+source['dimension'])]+[variables],dtype=np.int64)
    values=np.zeros((width,len(columns)),dtype=np.int32);cursor=0;blocks=[]
    for source in sources:
        s=source['index'];basis=root/'basis'/('profile%02d'%s)
        record=json.loads((basis/'result.json').read_text())
        assert record['status']=='complete_rational_profile_basis_recipes' and record['rank_lower_bound']==record['dimension']==source['dimension']
        assert sha(basis/'independent_seeds.pkl')==record['seeds_sha256']
        seeds=pickle.loads((basis/'independent_seeds.pkl').read_bytes());assert len(seeds)==source['dimension']
        covered=np.zeros(source['dimension'],dtype=bool)
        for folder in sorted((root/'maps'/('profile%02d'%s)).iterdir()):
            manifest=json.loads((folder/'manifest.json').read_text())
            assert manifest['source']==source and manifest['prime']==PRIME
            assert manifest['targets'][str(index)]==rhs_manifest['frame']
            assert sha(folder/'source_recipes.pkl')==manifest['source_recipes_sha256']
            assert pickle.loads((folder/'source_recipes.pkl').read_bytes())==seeds[manifest['start_row']:manifest['stop_row']]
            for marker in sorted(folder.glob('t%02d_*.json'%index)):
                block=json.loads(marker.read_text());first,last=block['start'],block['stop']
                assert block['source']==s and block['target']==index and block['columns']==width
                assert manifest['start_row']<=first<last<=manifest['stop_row']
                assert not covered[first:last].any()
                path=marker.with_suffix('.npy');assert sha(path)==block['sha256']
                matrix=np.load(path,mmap_mode='r');assert matrix.shape==(last-first,width)
                values[:,cursor+first:cursor+last]=matrix.T;covered[first:last]=True
                blocks.append(dict(source=s,start=first,stop=last,sha256=block['sha256']))
        if not covered.all():raise FileNotFoundError('Target %d source %d has %d missing rows'%(index,s,int((~covered).sum())))
        cursor+=source['dimension']
    # The producer already stores -Q1 x. Do not negate it a second time.
    values[:,-1]=np.load(rhs_root/'rhs.npy')
    manifest=dict(target=index,variables=variables,columns=width,prime=PRIME,blocks=blocks,
                  rhs_sha256=rhs_record['sha256'],rhs_frame=rhs_manifest['frame'])
    return values,columns,manifest


def homogeneous(stages,variables,free_columns,prime=PRIME):
    """Materialize any requested subset of the saved free-parameter basis."""
    pivots={p for stage in stages for p in stage['pivots']}
    assert len(set(free_columns))==len(free_columns) and not pivots.intersection(free_columns)
    values=np.zeros((variables,len(free_columns)),dtype=np.int64)
    values[free_columns,np.arange(len(free_columns))]=1
    for stage in reversed(stages):
        rows,columns=stage['rows'],stage['columns']
        assert columns[-1]==variables
        values[stage['pivots']]=(-matmul_mod(rows[:,:-1].astype(np.int64),values[columns[:-1]],prime))%prime
    return values.T


def run(root,output,plan_path=None,seconds=1800):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=True)
    plan_path=plan_path or HERE/'su3_correction_profiles_v1.json';plan=json.loads(plan_path.read_text())
    _,variables=offsets(plan)
    validation=json.loads((HERE/'bounded_flint_rref_validation.json').read_text())
    assert validation['status']=='passed_bounded_flint_rref'
    for name,digest in validation['source_sha256'].items():assert sha(HERE/name)==digest
    stages=[];processed=[]
    for target in sorted(plan['targets'],key=lambda x:(x['dimension'],x['index'])):
        folder=output/('stage%02d_target%02d'%(len(stages),target['index']));folder.mkdir(exist_ok=True)
        try:values,columns,manifest=assemble(root,plan,target)
        except FileNotFoundError as error:
            save(output/'result.json',dict(status='waiting_for_complete_target_blocks',processed_targets=processed,error=str(error),correction_solved=False));return
        manifest.update(previous_stages=[s['certificate_sha256'] for s in stages],driver_sha256=sha(__file__),plan_sha256=sha(plan_path))
        marker=folder/'result.json'
        if marker.exists():
            record=json.loads(marker.read_text());assert record['manifest']==manifest
            assert sha(folder/'rows.npy')==record['rows_sha256'] and sha(folder/'columns.npy')==record['columns_sha256']
            stage=dict(rows=np.load(folder/'rows.npy',mmap_mode='r'),columns=np.load(folder/'columns.npy'),pivots=record['pivots'],rank=record['rank'],consistent=record['consistent'])
        else:
            stage=reduce_target(values,columns,stages,variables,PRIME)
            save_array(folder/'rows.npy',stage['rows']);save_array(folder/'columns.npy',stage['columns'])
            record=dict(manifest=manifest,pivots=stage['pivots'],rank=stage['rank'],consistent=stage['consistent'],seconds=stage['seconds'],
                        rows_sha256=sha(folder/'rows.npy'),columns_sha256=sha(folder/'columns.npy'))
            save(marker,record)
        del values
        stage['certificate_sha256']=sha(marker);stages.append(stage);processed.append(target['index'])
        save(output/'progress.json',dict(processed_targets=processed,rank=sum(s['rank'] for s in stages),variables=variables,seconds=time.monotonic()-started))
        if not stage['consistent']:
            save(output/'result.json',dict(status='inconsistent_modular_evaluation_equations',processed_targets=processed,correction_solved=False,not_a_rational_lift_certificate=True));return
        if time.monotonic()-started>=seconds and len(stages)<len(plan['targets']):
            save(output/'result.json',dict(status='partial_coupled_elimination',processed_targets=processed,correction_solved=False));return
    solution,free=back_substitute(stages,variables,PRIME)
    residuals=[]
    for target in plan['targets']:
        values,columns,_=assemble(root,plan,target)
        residual=(matmul_mod(values[:,:-1].astype(np.int64),solution[columns[:-1],None],PRIME)[:,0]-values[:,-1])%PRIME
        assert not residual.any(),target['index']
        residuals.append(dict(target=target['index'],coordinates=len(residual),nonzero=0))
    save_array(output/'particular.npy',solution.astype(np.uint32));save_array(output/'free_variables.npy',np.asarray(free,dtype=np.int64))
    result=dict(status='complete_modular_evaluation_solution_with_all_free_parameters',prime=PRIME,variables=variables,
                equation_rank=sum(s['rank'] for s in stages),homogeneous_dimension=len(free),processed_targets=processed,
                residuals=residuals,particular_sha256=sha(output/'particular.npy'),free_variables_sha256=sha(output/'free_variables.npy'),
                stage_certificates=[s['certificate_sha256'] for s in stages],seconds=time.monotonic()-started,
                rational_correction_certified=False,target_frame_injectivity_certified=False,
                scope='All stored coupled evaluation equations are solved. The triangular stages parameterize every homogeneous solution; rational lifting and polynomial identities remain separate.')
    save(output/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--plan-path',type=Path);parser.add_argument('--seconds',type=int,default=1800)
    run(**vars(parser.parse_args()))
