"""Checkpointed elimination of complete target equations, preserving every free variable.

Each saved stage expresses new pivot variables through later pivots and final
free variables. Earlier pivot rows are never discarded or falsely frozen.
All arithmetic is exact modulo the pinned prime. Rational lifting is separate.
"""
import argparse
from datetime import datetime,timezone
import gc
import json
from pathlib import Path
import tempfile
import time

import sp_environment
import numpy as np
from flint import nmod_mat
from exact_blas_rows import matmul_mod
from integer_j30_boundary_system import HERE,PRIME,digest
from projected_checkpoint import atomic,json_bytes


def save_array(path,values):
    temporary=path.with_suffix('.tmp')
    with temporary.open('wb') as stream:np.save(stream,values,allow_pickle=False)
    temporary.replace(path)


def compact(values,columns,variables):
    keep=np.any(values!=0,axis=0)|(columns==variables)
    return values[:,keep],columns[keep]


def reduce_target(values,columns,previous,variables,prime):
    start=time.monotonic()
    assert values.dtype==np.int32 and columns[-1]==variables
    values,columns=compact(values,columns,variables)
    for number,stage in enumerate(previous):
        lookup={int(c):i for i,c in enumerate(columns)}
        common=[(j,lookup[p]) for j,p in enumerate(stage['pivots']) if p in lookup]
        if not common:continue
        coefficients=values[:,[j for i,j in common]]
        active=np.any(coefficients!=0,axis=0)
        if not np.any(active):continue
        indices=np.flatnonzero(active)
        coefficients=coefficients[:,indices]
        rows=stage['rows'][[common[i][0] for i in indices]]
        union=np.union1d(columns,stage['columns'])
        expanded=np.zeros((len(values),len(union)),dtype=np.int32)
        expanded[:,np.searchsorted(union,columns)]=values
        target_columns=np.searchsorted(union,stage['columns'])
        # Bounded row batches cap binary64 temporaries. The helper proves every
        # integer dot-product subtotal is strictly below 2**53.
        for begin in range(0,len(values),256):
            end=min(begin+256,len(values))
            product=matmul_mod(coefficients[begin:end],rows,prime)
            expanded[begin:end,target_columns]=(expanded[begin:end,target_columns]-product)%prime
        del values,coefficients,rows
        values,columns=compact(expanded,union,variables)
        del expanded
        remaining=set(map(int,columns))
        assert not remaining.intersection(stage['pivots'])
        print(json.dumps(dict(stage='previous_pivots_eliminated',previous_stage=number,
                              equations=len(values),columns=len(columns),seconds=time.monotonic()-start)),flush=True)
    values=values[np.any(values!=0,axis=1)]
    if not len(values):
        return dict(rows=np.zeros((0,len(columns)),dtype=np.int32),columns=columns,pivots=[],
                    consistent=True,rank=0,seconds=time.monotonic()-start)
    exact=nmod_mat(values.tolist(),prime)
    del values
    gc.collect()
    reduced,rank=exact.rref()
    del exact
    rank=int(rank)
    rows=np.asarray(reduced.tolist(),dtype=np.int32)[:rank].copy()
    del reduced
    gc.collect()
    local_pivots=[int(np.flatnonzero(row)[0]) for row in rows]
    assert np.array_equal(rows[:,local_pivots],np.eye(rank,dtype=np.int32))
    pivots=[int(columns[i]) for i in local_pivots]
    assert not set(pivots).intersection(p for s in previous for p in s['pivots'])
    return dict(rows=rows,columns=columns,pivots=pivots,
                consistent=not pivots or pivots[-1]!=variables,rank=rank,seconds=time.monotonic()-start)


def back_substitute(stages,variables,prime,seed=None):
    pivots={p for stage in stages for p in stage['pivots']}
    assert variables not in pivots
    free=sorted(set(range(variables))-pivots)
    values=np.zeros(variables,dtype=np.int64)
    if seed is not None:
        values[free]=np.random.default_rng(seed).integers(0,prime,size=len(free))
    for stage in reversed(stages):
        rows=stage['rows'];columns=stage['columns']
        assert columns[-1]==variables
        assert len(columns)*(prime-1)**2<2**63
        for begin in range(0,len(rows),256):
            end=min(begin+256,len(rows))
            rhs=rows[begin:end,-1].astype(np.int64) if seed is None else 0
            product=rows[begin:end,:-1] @ values[columns[:-1]]
            values[stage['pivots'][begin:end]]=(rhs-product)%prime
    return values,free


def controls():
    start=time.monotonic();cases=0;coordinates=0
    for seed in range(12):
        random=np.random.default_rng(seed+7751)
        variables=17;equations=27
        if seed==0:
            matrix=np.zeros((equations,variables),dtype=np.int32)
        else:
            left=random.integers(0,PRIME,size=(equations,8),dtype=np.int64)
            right=random.integers(0,PRIME,size=(8,variables),dtype=np.int64)
            matrix=((left@right)%PRIME).astype(np.int32)
            matrix[:,::5]=0
        rhs=(matrix.astype(np.int64)@random.integers(0,PRIME,size=variables))%PRIME
        augmented=np.column_stack((matrix,rhs)).astype(np.int32)
        if seed%3==2:
            augmented[-1,:]=0;augmented[-1,-1]=1
        order=random.permutation(equations)
        stages=[]
        for begin in range(0,equations,6):
            stage=reduce_target(augmented[order[begin:begin+6]].copy(),np.arange(variables+1),stages,variables,PRIME)
            stages.append(stage)
            if not stage['consistent']:break
        consistent=all(s['consistent'] for s in stages)
        map_rank=int(nmod_mat(augmented[:,:-1].tolist(),PRIME).rank())
        augmented_rank=int(nmod_mat(augmented.tolist(),PRIME).rank())
        assert consistent==(map_rank==augmented_rank)
        if consistent:
            assert sum(s['rank'] for s in stages)==map_rank
            solution,free=back_substitute(stages,variables,PRIME)
            assert len(free)==variables-map_rank
            assert not np.any((matrix.astype(np.int64)@solution-rhs)%PRIME)
            coordinates+=equations
            for trial in (137,619,2027):
                solution,free=back_substitute(stages,variables,PRIME,trial)
                assert not np.any(matrix.astype(np.int64)@solution%PRIME)
                coordinates+=equations
        cases+=1
    # Later stages can introduce lower-numbered variables and must inherit
    # columns created by substitution, even if absent from their raw equations.
    variables=17
    matrix=np.zeros((6,variables),dtype=np.int32)
    matrix[0,[8,14]]=[1,2]
    matrix[1,[2,8]]=[1,3]
    matrix[2,[0,2,11]]=[2,1,4]
    matrix[3,[11,14]]=[1,1]
    matrix[4]=matrix[1]+matrix[3]
    rhs=matrix.astype(np.int64)@np.arange(1,variables+1,dtype=np.int64)%PRIME
    augmented=np.column_stack((matrix,rhs)).astype(np.int32)
    stages=[]
    for row in augmented:
        stages.append(reduce_target(row[None,:].copy(),np.arange(variables+1),stages,variables,PRIME))
    assert [p for s in stages for p in s['pivots']]==[8,2,0,11]
    for seed in (None,137,619,2027):
        solution,free=back_substitute(stages,variables,PRIME,seed)
        expected=rhs if seed is None else 0
        assert not np.any((matrix.astype(np.int64)@solution-expected)%PRIME)
        coordinates+=len(matrix)
    cases+=1
    record=dict(status='passed',independent_dense_flint_cases=cases,solution_residual_coordinates=coordinates,
                source_sha256=digest(Path(__file__)),blas_source_sha256=digest(HERE/'exact_blas_rows.py'),
                seconds=time.monotonic()-start)
    atomic(HERE/'staged_solver_controls.json',json_bytes(record))
    print(json.dumps(record),flush=True)
    return record


def target_input(audit,target):
    key=[target['degree'],target['number']]
    edges=[e for e in audit['edges'] if e['target']==key]
    assert edges and all(e['complete'] for e in edges) and target['complete']
    sources=[s for s in audit['sources'] if any(e['source']==[s['degree'],s['number']] for e in edges)]
    assert all(s['certificate'] for s in sources)
    manifest=dict(target=key,variables=audit['full_source_dimension'],prime=PRIME,
                  columns=target['columns'],rhs_sha256=target['vector_sha256'],
                  sources=[dict(source=[s['degree'],s['number']],offset=s['offset'],dimension=s['dimension'],
                                rational_source_sha256=s['certificate']['rational_basis_sha256']) for s in sources],
                  blocks=[{k:e[k] for k in ('source','target','folder','source_dimension','matrix_sha256','manifest_sha256')} for e in edges])
    return manifest,sources,edges


def assemble_target(manifest,sources,edges):
    variables=manifest['variables']
    columns=np.asarray([i for s in sources for i in range(s['offset'],s['offset']+s['dimension'])]+[variables],dtype=np.int64)
    values=np.zeros((manifest['columns'],len(columns)),dtype=np.int32)
    offset=0
    for source in sources:
        edge=next(e for e in edges if e['source']==[source['degree'],source['number']])
        path=HERE/'sp2_j30_first_correction_blocks'/edge['folder']/'matrix.npy'
        assert digest(path)==edge['matrix_sha256']
        block=np.load(path,mmap_mode='r',allow_pickle=False)
        values[:,offset:offset+source['dimension']]=block.T
        offset+=source['dimension']
    rhs=HERE/'sp2_j30_first_correction_rhs'/('target%d_%d'%tuple(manifest['target']))/'vector.npy'
    assert digest(rhs)==manifest['rhs_sha256']
    values[:,-1]=(-np.load(rhs,allow_pickle=False))%PRIME
    return values,columns


def load_stage(folder,manifest):
    record=json.loads((folder/'result.json').read_text())
    assert record['status']=='passed' and record['input']==manifest
    for name in ('rows','columns'):
        assert digest(folder/(name+'.npy'))==record[name+'_sha256']
    return dict(rows=np.load(folder/'rows.npy',mmap_mode='r',allow_pickle=False),
                columns=np.load(folder/'columns.npy',allow_pickle=False),pivots=record['pivots'],
                consistent=record['consistent'],rank=record['rank'],record=record,
                certificate_sha256=digest(folder/'result.json'))


def full_residual(audit,coefficients):
    sources={(s['degree'],s['number']):s for s in audit['sources']}
    records=[]
    for target in audit['targets']:
        key=[target['degree'],target['number']]
        values=np.zeros((len(coefficients),target['columns']),dtype=np.int64)
        for edge in (e for e in audit['edges'] if e['target']==key):
            source=sources[tuple(edge['source'])]
            block=np.load(HERE/'sp2_j30_first_correction_blocks'/edge['folder']/'matrix.npy',mmap_mode='r',allow_pickle=False)
            assert source['dimension']*(PRIME-1)**2<2**63
            part=coefficients[:,source['offset']:source['offset']+source['dimension']]
            for begin in range(0,target['columns'],512):
                end=min(begin+512,target['columns'])
                values[:,begin:end]=(values[:,begin:end]+part@block[:,begin:end])%PRIME
        rhs=np.load(HERE/'sp2_j30_first_correction_rhs'/('target%d_%d'%tuple(key))/'vector.npy',allow_pickle=False)
        values[0]=(values[0]+rhs)%PRIME
        nonzero=int(np.count_nonzero(values))
        assert nonzero==0,(key,nonzero)
        records.append(dict(target=key,coordinates=values.size,nonzero=nonzero))
    return records


def run(seconds,max_stages,watch):
    start=time.monotonic()
    from assemble_sp2_first_correction import main as assemble
    audit=assemble()
    out=HERE/'sp2_j30_staged_correction_modular'
    out.mkdir(exist_ok=True)
    variables=audit['full_source_dimension']
    control=json.loads((HERE/'staged_solver_controls.json').read_text())
    assert control['status']=='passed' and control['source_sha256']==digest(Path(__file__))
    order=sorted(audit['targets'],key=lambda t:t['columns'])
    stages=[];processed=[]
    for target in order:
        while not all(e['complete'] for e in audit['edges'] if e['target']==[target['degree'],target['number']]):
            if not watch or time.monotonic()-start>=seconds:return checkpoint(out,stages,processed,variables,start)
            checkpoint(out,stages,processed,variables,start)
            time.sleep(15)
            audit=assemble()
        manifest,sources,edges=target_input(audit,target)
        manifest.update(driver_sha256=digest(Path(__file__)),blas_sha256=digest(HERE/'exact_blas_rows.py'),
                        previous_certificates=[s['certificate_sha256'] for s in stages])
        folder=out/('stage%02d_target%d_%d'%(len(stages),*manifest['target']))
        folder.mkdir(exist_ok=True)
        if not (folder/'result.json').exists():
            values,columns=assemble_target(manifest,sources,edges)
            reduced=reduce_target(values,columns,stages,variables,PRIME)
            save_array(folder/'rows.npy',reduced['rows'])
            save_array(folder/'columns.npy',reduced['columns'])
            record=dict(status='passed',input=manifest,pivots=reduced['pivots'],rank=reduced['rank'],
                        consistent=reduced['consistent'],rows_sha256=digest(folder/'rows.npy'),
                        columns_sha256=digest(folder/'columns.npy'),seconds=reduced['seconds'])
            atomic(folder/'result.json',json_bytes(record))
            del values,columns,reduced
            gc.collect()
        stage=load_stage(folder,manifest)
        stages.append(stage);processed.append(manifest['target'])
        progress=checkpoint(out,stages,processed,variables,start)
        if not stage['consistent'] or (max_stages and len(stages)>=max_stages):return progress
        if time.monotonic()-start>=seconds:return progress
    assert len(stages)==len(audit['targets'])==13
    solution,free=back_substitute(stages,variables,PRIME)
    homogeneous=[back_substitute(stages,variables,PRIME,seed)[0] for seed in (137,619,2027)]
    residuals=full_residual(audit,np.asarray([solution]+homogeneous))
    save_array(out/'solution.npy',solution)
    result=checkpoint(out,stages,processed,variables,start)
    result.update(status='complete_modular_solution',equation_consistent=True,free_columns=free,
                  all_homogeneous_freedoms_retained=True,solution_sha256=digest(out/'solution.npy'),
                  residual_checks=residuals,checked_solution_coordinates=sum(t['columns'] for t in audit['targets']),
                  checked_homogeneous_coordinates=3*sum(t['columns'] for t in audit['targets']),
                  homogeneous_check_seeds=[137,619,2027])
    atomic(out/'result.json',json_bytes(result))
    print(json.dumps({k:v for k,v in result.items() if k not in ('free_columns','stage_certificates')}),flush=True)
    return result


def checkpoint(out,stages,processed,variables,start):
    consistent=all(s['consistent'] for s in stages)
    augmented=sum(s['rank'] for s in stages)
    rank=augmented-int(not consistent)
    result=dict(status='checkpoint' if consistent else 'modular_inconsistent',prime=PRIME,
                processed_targets=processed,complete_targets=len(processed),required_targets=13,
                variables=variables,processed_coefficient_rank=rank,processed_augmented_rank=augmented,
                remaining_free_variables=variables-rank,processed_equations_consistent=consistent,
                stage_certificates=[s['certificate_sha256'] for s in stages],
                rational_obstruction_decided=False,new_lift=False,
                updated_at=datetime.now(timezone.utc).isoformat(),seconds=time.monotonic()-start)
    atomic(out/'progress.json',json_bytes(result))
    print(json.dumps({k:v for k,v in result.items() if k!='stage_certificates'}),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--controls',action='store_true')
    parser.add_argument('--seconds',type=int,default=2400)
    parser.add_argument('--max-stages',type=int,default=0)
    parser.add_argument('--watch',action='store_true')
    args=parser.parse_args()
    controls() if args.controls else run(args.seconds,args.max_stages,args.watch)
