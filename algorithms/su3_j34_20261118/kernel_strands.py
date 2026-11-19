"""Identify the four homogeneous blocks and seven coupled correction families at SU(3) j=34."""
import argparse, hashlib, json
from pathlib import Path

HERE=Path(__file__).resolve().parent
STRANDS=((17,13,7),(18,14,8),(19,15,9),(22,18,12))
FAMILIES=(16,17,18,19,20,21,22)

def load(root):
    paths=[root/'classical_length9.json',root/'classical_order1.json']
    hashes=[hashlib.sha256(x.read_bytes()).hexdigest() for x in paths]
    assert hashes==['47e60de0fbcf1f2f08c8c61be4c53a168c5a97a952771348701c3e18a0c9726f','1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c']
    pre=json.loads(paths[0].read_bytes());plan=json.loads(paths[1].read_bytes())['incoming']
    assert pre['mu']==plan['mu'] and pre['parity']==plan['parity']==1 and pre['target_length']==plan['source_length']==10
    assert pre['targets']==[{k:v for k,v in x.items() if k!='targets'} for x in plan['sources']]
    return pre,plan

def describe(root):
    pre, plan=load(root); source=plan['sources'];target=plan['targets']; families={}
    for i in FAMILIES:
        pr=pre['sources'][i]; kept=[s for s in pr['targets'] if 12<=s<=18]; omitted=[s for s in pr['targets'] if s not in kept]
        left={t for s in kept for t in source[s]['targets']};right={t for s in omitted for t in source[s]['targets']}
        allowed=sorted(left&right); assert set(allowed)<={13,14,15}
        families[str(i)]=dict(predecessor_dimension=pr['dimension'],kept=kept,complement=omitted,shared_final_targets=allowed,only_kept_reaches=sorted(left-right))
    strands=[]
    for i,s,t in STRANDS:
        assert [j for j in range(12,19) if t in source[j]['targets']]==[s]
        def difference(a,b):
            x={tuple(k):n for k,n in a['profile']};y={tuple(k):n for k,n in b['profile']}
            return {k:y.get(k,0)-x.get(k,0) for k in x.keys()|y.keys() if y.get(k,0)!=x.get(k,0)}
        assert difference(pre['sources'][i],source[s])==difference(source[s],target[t])=={(1,0):2,(2,0):-1}
        strands.append(dict(predecessor=i,source=s,target=t,dimensions=[pre['sources'][i]['dimension'],source[s]['dimension'],target[t]['dimension']]))
    assert sum(source[i]['dimension'] for i in range(12,19))==95313
    return dict(strands=strands, families=families, restricted_source_dimension=95313)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=HERE);p.add_argument('--wrapped-plan',type=Path);a=p.parse_args()
    if a.wrapped_plan:
        pre,_=load(a.root);a.wrapped_plan.write_text(json.dumps({'order':1,'incoming':pre},indent=2)+'\n')
    print(json.dumps(describe(a.root),indent=2))
