"""Make the Q0 and diagonal restriction plan on the very same length-nine recipes."""
import argparse, hashlib, json, sys
from pathlib import Path
from collections import Counter

P=1000033
BASE='3cd48eea41430f6cd08449f651f729db8125417ce04b1d6b8dbd6f4b77867149'
PREV='b8c669fa244e851f38a25b91186c3091e98a939935242805b0f71a9f8c00d609'
def sha(b):return hashlib.sha256(b).hexdigest()
def frozen(x):return tuple((tuple(k),n) for k,n in x)

def make(root,cartan):
    rt=root/'runtime/work/sp3-native.nosync';sys.path.insert(0,str(rt))
    from su3_correction_profiles_v1 import profile_dimensions,contents,profile,Q_letter
    raw=(root/'SOURCE_MANIFEST.json').read_bytes();assert sha(raw)==BASE;manifest=json.loads(raw)
    raw=(rt/'su3_second_correction_profiles_v1.json').read_bytes();assert sha(raw)==PREV;prev=json.loads(raw);mu=tuple(prev['mu'])
    assert mu==(2,2,5,2,2) and prev['parity']==-1 and prev['target_length']==9
    source=prev['targets']; dims=profile_dimensions(mu,10); targets=[dict(index=i,fermions=k[0],profile=k[1],dimension=v) for i,(k,v) in enumerate(sorted(dims.items()))]
    assert len(source)==37 and sum(t['dimension'] for t in source)==63983 and sum(t['dimension'] for t in targets)==46876
    si={frozen(s['profile']):s['index'] for s in source};ti={t['profile']:t['index'] for t in targets};edges={i:set() for i in range(len(source))}
    # Q0 is a derivation: its profile replacement is independent of trace order and of basis projection.
    from functools import lru_cache
    @lru_cache(None)
    def splits(letter):return {(len(left[2]),sum(left[:2]),len(right[2]),sum(right[:2])) for _,left,right in Q_letter(letter)}
    for content in contents(mu,9):
        p=profile(content)
        if p not in si:continue
        for letter in set(content):
            kind=(len(letter[2]),sum(letter[:2]));counts=Counter(dict(p));counts[kind]-=1
            for n,a,m,b in splits(letter):
                cc=counts.copy();cc[n,a]+=1;cc[m,b]+=1;q=tuple(sorted((k,v) for k,v in cc.items() if v))
                if q in ti:edges[si[p]].add(ti[q])
    cart=json.loads(cartan.read_bytes());assert cart['mu']==list(mu) and cart['length']==9 and cart['dimension']==937
    cd={frozen(t['profile']):t['dimension'] for t in cart['profiles']};assert sum(cd.values())==937 and set(cd)<=set(si)
    for s in source:
        i=s['index'];d=f'second_basis/shard0/targets{i:02d}/';fr={}
        for f in ('independent_seeds.pkl','manifest.json','result.json'):
            b=(root/d/f).read_bytes();h=sha(b);assert h==manifest[d+f];fr[f]=h
            if f=='manifest.json':assert json.loads(b)['task']==s and json.loads(b)['length']==9
            if f=='result.json':assert json.loads(b)['seeds_sha256']==fr['independent_seeds.pkl'] and json.loads(b)['dimension']==s['dimension']
        s.update(targets=sorted(edges[i]),cartan_dimension=cd.get(frozen(s['profile']),0),recipe_sha256=fr['independent_seeds.pkl'])
    return dict(prime=P,mu=mu,parity=-1,source_length=9,target_length=10,source_dimension=63983,target_dimension=46876,cartan_dimension=937,source_manifest=BASE,previous_plan_sha256=PREV,cartan_dimensions_sha256=sha(cartan.read_bytes()),sources=source,targets=targets)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('root',type=Path);ap.add_argument('cartan_dimensions',type=Path);ap.add_argument('out',type=Path);a=ap.parse_args();assert __debug__
    x=make(a.root.resolve(),a.cartan_dimensions);a.out.write_text(json.dumps(x,indent=2)+'\n');print(len(x['targets']),sum(len(s['targets']) for s in x['sources']),sum(bool(s['cartan_dimension']) for s in x['sources']),sha(a.out.read_bytes()))
