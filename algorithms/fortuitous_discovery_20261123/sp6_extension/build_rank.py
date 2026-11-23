"""Re-expand stored Sp(5) identity recipes at a chosen rank.

Appending scalar pairs supplies genuine higher-rank identities.  Linear
independence is checked; no cohomology claim is inferred from transport.
"""
from fast_paths import *
import argparse
import hashlib
import json
import math
import pickle
import time

A=(0,0,(0,))
JETS=BASE_PACKAGE/'symplectic/sp5' if BASE_PACKAGE else HERE.parent/'fortuity_discovery_20261123/jets'


def load_baseline_target():
    original=HERE.parent/'fortuity_rank_extension_20261123/complete_target_sp5.pkl'
    if BASE_PACKAGE is None:return pickle.loads(original.read_bytes())
    rows,recipes=pickle.loads((JETS/'target_sp5_snapshot_181.pkl').read_bytes())
    def tuples(x):return tuple(tuples(y) for y in x) if isinstance(x,list) else x
    extra=json.loads((BASE_PACKAGE/'extension/additional_identities_sp5.json').read_text())
    for record in extra['identities']:
        polynomial=Counter();recipe=[]
        for coefficient,graph in record['terms']:
            graph=tuples(graph);recipe.append((graph,coefficient))
            if record['kind']=='paths':image=expand_paths(graph,5)
            else:
                paths,closed=graph;image=multiply(expand_paths(paths,5),{(closed,):1})
            for w,c in image.items():polynomial[w]+=coefficient*c
        row=pf.primitive(clean(polynomial))
        digest=hashlib.sha256(json.dumps([[str(row[w]),pf.m.encode(w)] for w in sorted(row)],separators=(',',':')).encode()).hexdigest()
        assert digest==record['polynomial_sha256']
        rows.append(row);recipes.append((record['kind'],recipe))
    return rows,recipes


def main():
    cache_native_canonicalization()
    ap=argparse.ArgumentParser()
    ap.add_argument('--rank',type=int,default=6)
    ap.add_argument('--seconds',type=int,default=600)
    args=ap.parse_args()
    assert args.rank>=5
    rank=args.rank;extra=rank-5;delta=2*extra
    output=HERE/f'frames_sp{rank}';output.mkdir(exist_ok=True)
    started=time.monotonic();all_rows={};qcache={};records=[]
    def qrow(length,index):
        key=(length,index)
        if key not in qcache:
            folder=output/f'Q{length}';folder.mkdir(exist_ok=True)
            path=folder/f'row_{index}.pkl'
            if path.exists():
                digest,row=pickle.loads(path.read_bytes())
                assert digest==hashlib.sha256(pickle.dumps(all_rows[length][index],protocol=4)).hexdigest()
            else:
                row=pf.m.q(all_rows[length][index])
                path.write_bytes(pickle.dumps((hashlib.sha256(pickle.dumps(all_rows[length][index],protocol=4)).hexdigest(),row),protocol=4))
            qcache[key]=row
        return qcache[key]

    expand=lru_cache(maxsize=256)(expand_paths)
    def paths_plus(paths):
        return tuple(paths)+((A,A),)*extra
    def graph_expand(paths,cycles=()):
        return multiply(expand(paths_plus(paths),rank),{tuple(cycles):1})
    def rebuild(recipe,length,index):
        if not isinstance(recipe[0],str):
            letters,remaining=recipe
            paths,cycles,_,_=graph_components(len(letters),5,remaining)
            raw=graph_expand(tuple(tuple(letters[i] for i in p) for p in paths),
                             tuple(tuple(letters[i] for i in c) for c in cycles))
            return pf.singlet(raw,3)
        kind=recipe[0]
        if kind in ('paired','paths','paths_cycle'):
            out=Counter()
            for graph,c in recipe[1]:
                if kind=='paths_cycle':paths,cycle=graph;row=graph_expand(paths,(cycle,))
                else:row=graph_expand(graph)
                for w,n in row.items():out[w]+=c*n
        elif kind=='Q':
            _,previous,j=recipe;out=qrow(previous+delta,j)
        elif kind=='Qcomponent':
            _,previous,j,profile=recipe
            profile=tuple(sorted(tuple(profile)+((0,1),)*delta))
            out=pf.profile_components(qrow(previous+delta,j)).get(profile,{})
        elif kind=='component':
            _,j,profile=recipe
            profile=tuple(sorted(tuple(profile)+((0,1),)*delta))
            out=pf.profile_components(all_rows[length][j]).get(profile,{})
        else:raise ValueError(kind)
        return pf.primitive(clean(out))

    for old_length,expected in [(12,67),(13,188),(14,183)]:
        length=old_length+delta
        old_rows,recipes=(pickle.loads((JETS/f'exact_relations_sp5_L{old_length}.pkl').read_bytes())
                         if old_length<14 else load_baseline_target())
        assert len(old_rows)==len(recipes)==expected
        rows=[];all_rows[length]=rows;echelon=pf.Echelon()
        folder=output/f'L{length}';folder.mkdir(exist_ok=True)
        for i,recipe in enumerate(recipes):
            path=folder/f'row_{i}.pkl'
            signature=hashlib.sha256(pickle.dumps((rank,recipe),protocol=4)).hexdigest()
            if path.exists():
                oldsignature,row=pickle.loads(path.read_bytes());assert signature==oldsignature
            else:
                row=rebuild(recipe,length,i)
                assert row and not pf.spin_action(row,True),('highest weight',length,i)
                if rank==5:assert row==old_rows[i],('baseline mismatch',length,i)
                path.write_bytes(pickle.dumps((signature,row),protocol=4))
            assert echelon.add(row),('dependent',rank,length,i)
            rows.append(row)
            if i%10==0 or i+1==expected:
                print(json.dumps(dict(rank=rank,length=length,rows=i+1,total=expected,terms=len(row),
                                      seconds=time.monotonic()-started)),flush=True)
            if time.monotonic()-started>args.seconds:
                print('Checkpoint time cap; rerun to continue.',flush=True);return
        target=output/f'relations_L{length}.pkl'
        target.write_bytes(pickle.dumps((rows,recipes),protocol=4))
        records.append(dict(length=length,dimension=expected,independent_modulo=1000003,
                            total_terms=sum(map(len,rows)),sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                            exact_highest_weight=True))
        qcache.clear();expand.cache_clear();cycle.cache_clear()
    result=dict(group=f'Sp({rank})',native_weight=[3,3,2*rank-1,0,0],frames=records,
                genuine_pfaffian_recipes=True,cohomology_not_yet_determined=True,
                baseline_exact_replay=rank==5,seconds=time.monotonic()-started)
    (output/'construction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
