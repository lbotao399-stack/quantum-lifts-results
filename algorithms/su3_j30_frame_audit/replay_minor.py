"""Rebuild a modest selected minor directly from the original saved rows; use FLINT det."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def replay(root,folder):
    from flint import nmod_mat
    load=lambda p:json.loads(p.read_bytes());r=load(folder/'result.json');m=load(folder/'manifest.json');assert r['manifest_sha256']==sha(folder/'manifest.json')
    w=r['witness'];rank=w['rank'];assert 0<rank<=1600 and len(set(w['columns']))==len(w['columns'])==len(set(map(tuple,w['rows'])))==rank
    assert min(w['columns'])>=0 and max(w['columns'])<r['width'];seeds={};images={};p=r['prime'];assert p==1000033 and m['prime']==p
    for g in r['chunks']:
        first,last=g['start'],g['stop'];arr=folder/('rows%05d_%05d.npy'%(first,last));assert sha(arr)==g['sha256'];v=np.load(arr,allow_pickle=False);assert v.shape==(last-first,r['width']);seeds.update((first+i,row) for i,row in enumerate(v))
    plan=load(root/'runtime/work/sp3-native.nosync/su3_correction_profiles_v1.json')
    for g in r['packed_Q0_sources']:
        i=g['source'];n=plan['sources'][i]['dimension'];folder1=root/'output/compacted/first/maps'/f'profile{i:02d}'/f'rows00000_{n:05d}';a=folder1/f"t{r['target']:02d}_00000_{n:05d}.npy";assert sha(a)==g['sha256'];images[i]=np.load(a,allow_pickle=False,mmap_mode='r')
    rows=[]
    for tag,*ind in w['rows']:
        if tag=='seed':row=seeds[ind[0]]
        else:assert tag=='image' and len(ind)==2;row=images[ind[0]][ind[1]]
        rows.append(row[w['columns']].tolist())
    det=int(nmod_mat(rows,p).det());assert det!=0 and det==w['determinant_in_selected_order']
    return dict(equation=r['equation'],target=r['target'],width=r['width'],dimension=r['dimension'],rank=rank,FLINT_determinant=det,certificate_sha256=sha(folder/'result.json'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('folder',type=Path);a=ap.parse_args();sys.path.insert(0,str(a.root/'runtime/work/sp3-native.nosync'));import sp_environment
    print(json.dumps(replay(a.root,a.folder)))
