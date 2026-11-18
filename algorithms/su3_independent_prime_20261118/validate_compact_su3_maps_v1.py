"""Compare original and packed coupled matrices, and reject damaged inputs."""
from pathlib import Path
import argparse,json,shutil,sys,tempfile,time
import numpy as np

parser=argparse.ArgumentParser();parser.add_argument('--runtime',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();sys.path.insert(0,str(args.runtime))
import validate_su3_second_coupled_solve_v1 as fixture
import su3_coupled_solve_v1 as first
import su3_second_coupled_solve_v1 as second
import compact_su3_maps_v1 as compact

started=time.monotonic();prime=first.PRIME
with tempfile.TemporaryDirectory(prefix='compact-SU3-control-') as temporary:
    root=Path(temporary);old=root/'old';new=root/'new'
    sources1=[dict(index=i,dimension=n,profile=['first',i],targets=[0,1]) for i,n in enumerate((3,2))]
    sources2=[dict(index=i,dimension=n,profile=['second',i],targets=[0,1]) for i,n in enumerate((2,4))]
    targets=[dict(index=0,dimension=2),dict(index=1,dimension=3)]
    plan1=dict(source_dimension=5,sources=sources1,targets=targets)
    plan2=dict(source_dimension=6,sources=sources2,targets=targets)
    quantum=dict(sources=sources1,targets=targets)
    frames={i:dict(columns=n,prime=prime,synthetic_control_only=True) for i,n in enumerate((3,4))}
    basis1={};basis2={}
    for item in sources1:
        p=old/'first/basis'/('profile%02d'%item['index']);fixture.basis(p,item);basis1[item['index']]=p
    for item in sources2:
        p=old/'basis/shard0'/('sources%02d'%item['index']);fixture.basis(p,item);basis2[item['index']]=p
    for item in targets:fixture.basis(old/'basis/shard0'/('targets%02d'%item['index']),item)
    rng=np.random.default_rng(20261118)
    matrix1=rng.integers(0,prime,size=(7,5),dtype=np.int64)
    matrixQ0=rng.integers(0,prime,size=(7,6),dtype=np.int64)
    matrixQ1=rng.integers(0,prime,size=(7,5),dtype=np.int64)
    fixture.maps(old/'first',sources1,frames,matrix1,basis1,'Q0')
    fixture.maps(old/'q0',sources2,frames,matrixQ0,basis2,'Q0')
    fixture.maps(old/'q1',sources1,frames,matrixQ1,basis1,'Q1')
    fixture.right_hand(old/'first',targets,frames,np.arange(7,dtype=np.int64),1)
    fixture.right_hand(old/'rhs',targets,frames,np.arange(7,dtype=np.int64)+7,2)
    records=[]
    for old_name,new_name,basis,plan,kind in [('first','first',old/'first/basis',plan1,'first_Q0'),('q0','q0',old/'basis/shard0',plan2,'second_Q0'),('q1','q1',old/'first/basis',quantum,'second_Q1')]:
        records.append(compact.run(old/old_name,basis,new/new_name,plan,kind,prime,workers=2))
    (new/'first/basis').symlink_to(old/'first/basis',target_is_directory=True)
    (new/'first/rhs').symlink_to(old/'first/rhs',target_is_directory=True)
    bases=second.locate_second_bases([old/'basis'],plan2);coordinates=0
    for target in targets:
        before=first.assemble(old/'first',plan1,target);after=first.assemble(new/'first',plan1,target)
        assert np.array_equal(before[0],after[0]) and np.array_equal(before[1],after[1]);coordinates+=before[0].size
        before=second.assemble_second(old/'q0',old/'q1',old/'rhs',old/'first/basis',bases,plan1,quantum,plan2,target)
        after=second.assemble_second(new/'q0',new/'q1',old/'rhs',old/'first/basis',bases,plan1,quantum,plan2,target)
        assert np.array_equal(before[0],after[0]) and np.array_equal(before[1],after[1]);coordinates+=before[0].size
    marker=next((old/'q1/maps/profile00').glob('*/t00_*.json'));data=marker.read_bytes();marker.unlink()
    try:compact.run(old/'q1',old/'first/basis',root/'missing',quantum,'second_Q1',prime,workers=1)
    except AssertionError as e:assert str(e)=='Missing source rows'
    else:raise AssertionError('Missing rows were accepted')
    marker.write_bytes(data)
    values=np.load(marker.with_suffix('.npy'));values[0,0]=(values[0,0]+1)%prime;np.save(marker.with_suffix('.npy'),values)
    try:compact.run(old/'q1',old/'first/basis',root/'damaged',quantum,'second_Q1',prime,workers=1)
    except AssertionError:pass
    else:raise AssertionError('Corrupted payload was accepted')
result=dict(status='passed_lossless_compaction_full_assembly_and_input_rejection',prime=prime,
    matrix_coordinates_compared=coordinates,input_blocks=sum(r['input_blocks'] for r in records),
    output_blocks=sum(r['output_blocks'] for r in records),checks=['first_and_second_assembled_matrices_identical',
    'source_order_and_global_columns_identical','missing_rows_rejected','corrupted_payload_rejected'],
    driver_sha256=compact.sha(compact.__file__),seconds=time.monotonic()-started)
args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
