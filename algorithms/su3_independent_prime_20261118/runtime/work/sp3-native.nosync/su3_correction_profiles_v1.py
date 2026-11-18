"""Exact C-odd cochain dimensions and all Q0 profile couplings for Q0 y=-Q1 x."""
import sp_environment
from su3_charge_dimensions_v1 import parity_parts, dimensions
from g2_full_classical_screen_v1 import contents
from sp_full_profiles import profile
from fortuity.full.letters import is_ferm, Q_letter
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path
import hashlib
import json
import time

MU=(1,1,6,3,3)


def profile_dimensions(mu,length):
    answer=defaultdict(int);rho=(2,1,0)
    for spin_shift,spin_sign in ((0,1),(1,-1)):
        for order in permutations(range(3)):
            sign=(-1)**sum(order[i]>order[j] for i in range(3) for j in range(i+1,3))
            shifted=(mu[0]+spin_shift,mu[1]-spin_shift)+tuple(mu[i+2]+rho[i]-rho[order[i]] for i in range(3))
            if min(shifted)<0:continue
            for content in contents(shifted,length):
                counts=Counter(content)
                bosons=tuple(sorted((n for L,n in counts.items() if not is_ferm(L)),reverse=True))
                fermions=tuple(sorted((n for L,n in counts.items() if is_ferm(L)),reverse=True))
                answer[sum(fermions),profile(content)]+=spin_sign*sign*parity_parts(bosons,fermions)[1]
    assert all(value>=0 for value in answer.values())
    answer={key:value for key,value in answer.items() if value}
    totals=defaultdict(int)
    for (degree,_),value in answer.items():totals[degree]+=value
    assert dict(totals)==dimensions(mu,length,-1)
    return answer


def main():
    start=time.monotonic()
    source=profile_dimensions(MU,9);target=profile_dimensions(MU,10)
    source_by_profile={p:(i,d,n) for i,((d,p),n) in enumerate(sorted(source.items()))}
    target_by_profile={p:(i,d,n) for i,((d,p),n) in enumerate(sorted(target.items()))}
    edges=defaultdict(set)
    for content in contents(MU,9):
        p=profile(content)
        if p not in source_by_profile:continue
        i=source_by_profile[p][0]
        for position,L in enumerate(content):
            for _,left,right in Q_letter(L):
                q=profile(content[:position]+content[position+1:]+(left,right))
                if q in target_by_profile:edges[i].add(target_by_profile[q][0])
    graph=defaultdict(set)
    for i,js in edges.items():
        for j in js:graph['s%d'%i].add('t%d'%j);graph['t%d'%j].add('s%d'%i)
    remaining=set(graph);components=[]
    while remaining:
        todo=[min(remaining)];visited=set()
        while todo:
            vertex=todo.pop()
            if vertex in visited:continue
            visited.add(vertex);todo.extend(graph[vertex]-visited)
        remaining-=visited
        ss=sorted(int(x[1:]) for x in visited if x[0]=='s')
        tt=sorted(int(x[1:]) for x in visited if x[0]=='t')
        components.append(dict(source_profiles=ss,target_profiles=tt,
            source_dimension=sum(n for _,(i,d,n) in source_by_profile.items() if i in ss),
            target_dimension=sum(n for _,(i,d,n) in target_by_profile.items() if i in tt)))
    result=dict(status='exact_integer_profiles_and_structural_couplings',mu=MU,parity=-1,
        source_dimension=sum(source.values()),target_dimension=sum(target.values()),
        sources=[dict(index=i,fermions=d,dimension=n,profile=p,targets=sorted(edges[i]))
                 for p,(i,d,n) in sorted(source_by_profile.items(),key=lambda x:x[1])],
        targets=[dict(index=i,fermions=d,dimension=n,profile=p)
                 for p,(i,d,n) in sorted(target_by_profile.items(),key=lambda x:x[1])],
        connected_components=components,seconds=time.monotonic()-start,
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='All structural Q0 couplings are retained. No matrix rank, correction solution or lift is inferred.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(source_dimension=result['source_dimension'],target_dimension=result['target_dimension'],
        source_profiles=len(source),target_profiles=len(target),components=components,seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
