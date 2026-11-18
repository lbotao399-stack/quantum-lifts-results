"""Verify checkpoint blocks once and pack each source-target matrix contiguously.

The original row order, recipes, evaluation frames, and residues are preserved.
No rank reduction or choice of correction is performed.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,hashlib,json,pickle,time
import numpy as np

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda:stream.read(4*1024**2),b''):digest.update(data)
    return digest.hexdigest()

def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def compact_source(root,basis,output,source,kind,prime):
    started=time.monotonic();index=source['index'];dimension=source['dimension']
    folder=output/'maps'/('profile%02d'%index)/('rows%05d_%05d'%(0,dimension))
    folder.mkdir(parents=True,exist_ok=True)
    declared=json.loads((basis/'result.json').read_text());seed_path=basis/'independent_seeds.pkl'
    assert declared['dimension']==declared['rank_lower_bound']==dimension
    assert sha(seed_path)==declared['seeds_sha256'];seeds=pickle.loads(seed_path.read_bytes());assert len(seeds)==dimension
    parts=[];frames=None;template=None;input_manifests=[]
    for part in sorted((root/'maps'/('profile%02d'%index)).iterdir()):
        if not part.is_dir() or not (part/'manifest.json').exists():continue
        manifest=json.loads((part/'manifest.json').read_text())
        assert manifest['source']==source and manifest['prime']==prime
        assert 0<=manifest['start_row']<manifest['stop_row']<=dimension
        if frames is None:frames=manifest['targets'];template=manifest
        else:assert manifest['targets']==frames,'Incompatible evaluation frames'
        assert set(map(int,frames))==set(source['targets'])
        recipe_path=part/'source_recipes.pkl';assert sha(recipe_path)==manifest['source_recipes_sha256']
        assert pickle.loads(recipe_path.read_bytes())==seeds[manifest['start_row']:manifest['stop_row']]
        if kind=='Q1':assert manifest['source_basis_sha256']==declared['seeds_sha256'] and manifest['points_override']==0
        parts.append((part,manifest))
        input_manifests.append(dict(path=str((part/'manifest.json').relative_to(root)),sha256=sha(part/'manifest.json')))
    assert parts,'Missing source chunks'
    matrices={};covered={};payloads=[]
    for target in source['targets']:
        width=frames[str(target)]['columns'];assert frames[str(target)]['prime']==prime
        path=folder/('t%02d_%05d_%05d.tmp.npy'%(target,0,dimension))
        matrices[target]=np.lib.format.open_memmap(path,mode='w+',dtype=np.uint32,shape=(dimension,width))
        covered[target]=np.zeros(dimension,dtype=bool)
    for part,manifest in parts:
        for marker in sorted(part.glob('t*.json')):
            record=json.loads(marker.read_text());target=record['target'];first,last=record['start'],record['stop']
            assert record['source']==index and target in matrices
            assert manifest['start_row']<=first<last<=manifest['stop_row']
            assert not covered[target][first:last].any(),'Overlapping source rows'
            width=frames[str(target)]['columns'];assert record['columns']==width
            path=marker.with_suffix('.npy');assert sha(path)==record['sha256']
            values=np.load(path,mmap_mode='r');assert values.shape==(last-first,width)
            assert np.all((values>=0)&(values<prime))
            matrices[target][first:last]=values;covered[target][first:last]=True
            payloads.append(dict(path=str(path.relative_to(root)),sha256=record['sha256'],target=target,start=first,stop=last))
    assert all(rows.all() for rows in covered.values()),'Missing source rows'
    records=[]
    for target in source['targets']:
        matrices[target].flush();del matrices[target]
        name='t%02d_%05d_%05d'%(target,0,dimension);path=folder/(name+'.npy')
        (folder/(name+'.tmp.npy')).replace(path)
        record=dict(source=index,target=target,start=0,stop=dimension,columns=frames[str(target)]['columns'],sha256=sha(path))
        save(folder/(name+'.json'),record);records.append(record)
    recipe_path=folder/'source_recipes.pkl';recipe_path.write_bytes(seed_path.read_bytes())
    manifest=dict(template,start_row=0,stop_row=dimension,rows=dimension,batch=dimension,
        source_recipes_sha256=sha(recipe_path),production_driver_sha256=template.get('driver_sha256'),
        driver_sha256=sha(__file__),compaction_only=True)
    save(folder/'manifest.json',manifest)
    save(folder/'input_payloads.json',dict(manifests=input_manifests,payloads=payloads))
    result=dict(status='complete_profile_compaction',kind=kind,source=index,prime=prime,rows=dimension,
        input_blocks=len(payloads),output_blocks=len(records),records=records,
        input_payloads_sha256=sha(folder/'input_payloads.json'),seconds=time.monotonic()-started,
        scope='Every stored residue and its row order is preserved; no algebraic conclusion is added.')
    save(folder/'result.json',result);return result

def run(root,basis,output,plan,kind,prime,workers=4):
    started=time.monotonic();root=Path(root);basis=Path(basis);output=Path(output)
    assert output.resolve()!=root.resolve();output.mkdir(parents=True,exist_ok=True)
    finished=[]
    def task(source):
        name=('sources%02d' if kind=='second_Q0' else 'profile%02d')%source['index']
        return compact_source(root,basis/name,output,source,'Q1' if kind=='second_Q1' else 'Q0',prime)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(task,s) for s in plan['sources']]):
            finished.append(future.result());save(output/'progress.json',dict(completed=len(finished),required=len(plan['sources'])))
    result=dict(status='complete_lossless_verified_map_compaction',prime=prime,kind=kind,
        input_blocks=sum(r['input_blocks'] for r in finished),output_blocks=sum(r['output_blocks'] for r in finished),
        rows=sum(r['rows'] for r in finished),profiles=len(finished),seconds=time.monotonic()-started,
        driver_sha256=sha(__file__),profiles_records=finished)
    save(output/'compaction.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('root','basis','output','plan'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--kind',choices=['first_Q0','second_Q0','second_Q1'],required=True)
    p.add_argument('--prime',type=int,required=True);p.add_argument('--workers',type=int,default=4)
    args=vars(p.parse_args());args['plan']=json.loads(args['plan'].read_text())
    result=run(**args);print(json.dumps({k:v for k,v in result.items() if k!='profiles_records'}))
