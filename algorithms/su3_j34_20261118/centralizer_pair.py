"""Two incident parents vanish identically per target. Evaluate the residual over Z."""
import argparse
from fractions import Fraction
from pathlib import Path
from random import Random
import commuting_parent_test as old

s,ft=old.s,old.ft

def block_assignment(t,seed):
    rng=Random(seed); ids=bytes.fromhex(t); n=sum(ft.ODD[k] for k in ids); M={}
    for k in sorted(set(ids)):
        for h in (range(n) if ft.ODD[k] else [-1]):
            A=ft.traceless(rng,3)
            if k==0: A=((1,0,0),(0,1,0),(0,0,-2))
            elif k in (4,6): A=tuple(tuple(a if (i==2)==(j==2) else 0 for j,a in enumerate(row)) for i,row in enumerate(A))
            assert sum(row[i] for i,row in enumerate(A))==0
            M[k,h]=A
    A,B,C=[M[k,-1] for k in (0,4,6)]
    assert ft.mm(A,B)==ft.mm(B,A) and ft.mm(A,C)==ft.mm(C,A) and ft.mm(B,C)!=ft.mm(C,B)
    return M,n

def run(rhs,shared,frontier):
    z=s.read(frontier); vh,rows=old.residuals(rhs,shared,old.TARGETS); assert z['vector_sha']==vh and z['rhs']==s.INPUT
    records={}
    for t in old.TARGETS:
        p=old.split_pairs(t); item=z['target_parents'][t]; allowed=item['confined_to_failures']+item['touch_one_success']; assert set(item['all'])==set(p) and sorted((p[g][1:3] for g in allowed))==[(0,4),(0,6)]
        trials=[]
        for seed in (817234,419071,112903):
            M,n=block_assignment(t,seed); value=sum((c*ft.exterior_value(w,M,n) for w,c in rows[t].items()),Fraction());assert value and value.denominator==1
            trials.append(dict(seed=seed,exterior_degree=n,matrices={f'{i:02x}:{h}':m for (i,h),m in M.items()},integer_value=int(value)))
        records[t]=dict(splits={g:[*p[g]] for g in allowed},residual_words=len(rows[t]),assignments=trials)
    a,b=(records[t]['assignments'] for t in old.TARGETS);assert all(r['integer_value']==-v['integer_value'] for r,v in zip(a,b))
    return dict(rhs=s.INPUT,vector_sha=vh,targets=records,sign='Q0 O1 + Q1 O0',scope='All corrections supported on parents touching at most one previously successful content are excluded on each of these targets; no unrestricted obstruction follows.')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--shared',type=Path,default=s.HERE/'joint_70/vector_shared.json.gz');a.add_argument('--frontier',type=Path,default=s.HERE/'joint_70/parent_frontier.json');a.add_argument('--out',type=Path,required=True);v=a.parse_args();assert __debug__;z=run(v.rhs,v.shared,v.frontier);s.write(v.out,z);print({t:[a['integer_value'] for a in r['assignments']] for t,r in z['targets'].items()})
