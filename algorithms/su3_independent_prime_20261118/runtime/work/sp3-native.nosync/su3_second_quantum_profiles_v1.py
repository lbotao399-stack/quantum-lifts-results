"""Structural Q1 couplings from every first-correction profile."""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import json,time

from su3_second_correction_map_v1 import HERE,sha,contents
from sp_full_profiles import profile
from su_quantum_integer_v1 import species_rules


def main():
    started=time.monotonic()
    first=json.loads((HERE/'su3_correction_profiles_v1.json').read_text())
    second=json.loads((HERE/'su3_second_correction_profiles_v1.json').read_text())
    source_by_profile={tuple((tuple(key),n) for key,n in row['profile']):row['index'] for row in first['sources']}
    target_by_profile={tuple((tuple(key),n) for key,n in row['profile']):row['index'] for row in second['targets']}
    edges=defaultdict(set);cache={};content_count=0
    for content in contents(tuple(first['mu']),9):
        current=profile(content)
        if current not in source_by_profile:continue
        content_count+=1;index=source_by_profile[current]
        for i,j in combinations(range(len(content)),2):
            for key in ((content[i],content[j]),(content[j],content[i])):
                if key not in cache:cache[key]=species_rules(key,1)
                for left,right in cache[key]:
                    following=profile(tuple(letter for position,letter in enumerate(content) if position not in (i,j))+(left,right))
                    if following in target_by_profile:edges[index].add(target_by_profile[following])
    result=dict(status='all_first_correction_profile_structural_Q1_couplings',
        source_mu=first['mu'],target_mu=second['mu'],source_length=9,target_length=9,
        source_dimension=first['source_dimension'],target_dimension=second['target_dimension'],
        sources=[dict(row,targets=sorted(edges[row['index']])) for row in first['sources']],
        targets=second['targets'],structural_edges=sum(map(len,edges.values())),
        enumerated_contents=content_count,species_pairs=len(cache),seconds=time.monotonic()-started,
        sources_sha256={name:sha(HERE/name) for name in
            ('su3_second_quantum_profiles_v1.py','su3_correction_profiles_v1.json',
             'su3_second_correction_profiles_v1.json','su_quantum_integer_v1.py')},
        scope='Structural support only. All 28311 first-correction coordinates remain; the first equation restricts their homogeneous freedom to 10308 dimensions. No rank or lift is inferred.')
    (HERE/'su3_second_quantum_profiles_v1.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('sources','targets','sources_sha256')}),flush=True)


if __name__=='__main__':main()
