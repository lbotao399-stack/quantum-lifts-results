"""Compare the ordered sources, not target evaluation frames, at three primes."""
import argparse, hashlib, json, pickle, zipfile
from pathlib import Path

sha = lambda b: hashlib.sha256(b).hexdigest()
ROOT = Path(__file__).resolve().parents[2]
RUNTIME = 'runtime/work/sp3-native.nosync'
PINS = ('20ef1af3ba7b8d1cacc2bc66479b3bf25180df28630d439c41dca742b1cb7aeb',
        '95d55d6045d0a52a3f6e3755ab40e61f873d3b0d0712fde6a5c67227262f22af')
PACKAGES = ('3cd48eea41430f6cd08449f651f729db8125417ce04b1d6b8dbd6f4b77867149',
            'f155245a54cea2677007f73067865a39d8cc3be11df89e040420c2d71e786ada')

def compare(archives, trees):
    assert __debug__
    manifests = []
    for tree,pin in zip(trees,PACKAGES,strict=True):
        raw=(tree/'SOURCE_MANIFEST.json').read_bytes(); assert sha(raw)==pin
        manifests.append(json.loads(raw))
    result=[]; offset=0
    for stage,archive in enumerate(archives):
        assert sha(archive.read_bytes())==PINS[stage]
        tag=('su3_correction_profiles_v1.json','su3_second_correction_profiles_v1.json')[stage]
        plan_path=RUNTIME+'/'+tag
        plans=[(tree/plan_path).read_bytes() for tree in trees]
        assert plans[0]==plans[1] and all(sha(raw)==m[plan_path] for raw,m in zip(plans,manifests))
        plan=json.loads(plans[0]); count=(17,57)[stage];assert len(plan['sources'])==count
        with zipfile.ZipFile(archive) as z:
            for i,source in enumerate(plan['sources']):
                rel=('first_basis/profile%02d','second_basis/shard0/sources%02d')[stage]%i
                original=('basis/profile%02d','second_basis/sources%02d')[stage]%i
                for name in ('manifest.json','result.json','independent_seeds.pkl'):
                    raw=z.read(original+'/'+name); key=rel+'/'+name
                    for tree,m in zip(trees,manifests):
                        assert (tree/key).read_bytes()==raw and sha(raw)==m[key],key
                    if name=='manifest.json':
                        a=json.loads(raw);assert a['task']==source and a['mu']==plan['mu'] and a['parity']==-1 and a['length']==9-stage and a['plan_sha256']==sha(plans[0])
                    elif name=='result.json':
                        r=json.loads(raw); assert r['profile_index']==i==source['index'] and r['dimension']==r['rank_lower_bound']==source['dimension'] and r['fermions']==source['fermions']
                    else:
                        # Only unpickle after validation against the pinned archive and both package manifests.
                        assert sha(raw)==r['seeds_sha256'] and len(pickle.loads(raw))==source['dimension']
                result.append(dict(equation=stage+1,index=i,profile=source['profile'],length=a['length'],dimension=source['dimension'],start=offset,stop=offset+source['dimension'],seeds_sha256=r['seeds_sha256']))
                offset+=source['dimension']
        assert offset==(28311,77639)[stage]
    return dict(primes=[1000003,1000033,1000037],original_archives=list(PINS),production_manifests=list(PACKAGES),counts=[17,57],block_sizes=[28311,49328],variables=offset,comparison='Byte equality of seed lists, basis manifests and basis results; ordered global offsets in coupled solver convention, stop exclusive; does not compare target evaluations.',profiles=result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('first_zip',type=Path);p.add_argument('second_zip',type=Path);p.add_argument('third_tree',type=Path);p.add_argument('--second-tree',type=Path,default=ROOT/'algorithms/su3_independent_prime_20261118');p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=compare([a.first_zip,a.second_zip],[a.second_tree,a.third_tree]);a.output.write_text(json.dumps(r,indent=2)+'\n');print('Identical source recipes:',r['counts'],r['block_sizes'],r['variables'])
