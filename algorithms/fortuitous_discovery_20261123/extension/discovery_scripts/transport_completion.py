"""Re-expand the two missing graph recipes at a larger symplectic rank.

This transports identity recipes only. No transport of cohomology is
assumed. Every resulting polynomial is generated at the target rank.
"""
from pathlib import Path
from collections import Counter
import argparse
import json
import pickle
import time
import complete_target as ct

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,default=5)
    args=ap.parse_args();rank=args.rank;assert rank==5
    start=time.monotonic()
    source_rows,source_recipes=pickle.loads((HERE/'complete_target_sp4.pkl').read_bytes())
    rows,recipes=pickle.loads((ct.JETS/'target_sp5_snapshot_181.pkl').read_bytes())
    echelon=ct.pf.Echelon()
    for row in rows:assert echelon.add(row)
    scalar=(0,0,(0,));pair=(scalar,scalar)
    records=[]
    for source_index in (181,182):
        kind,old_recipe=source_recipes[source_index]
        polynomial=Counter();recipe=[]
        for graph,coefficient in old_recipe:
            if kind=='paths':
                sign,new_graph=ct.pc.canonical_paths(graph+(pair,))
                expanded=ct.pc.expand_paths(new_graph,rank)
            else:
                assert kind=='paths_cycle'
                paths,cycle=graph
                sign,new_paths=ct.pc.canonical_paths(paths+(pair,))
                new_graph=(new_paths,cycle)
                expanded=ct.expand_cycle(new_graph,rank)
            assert sign
            recipe.append((new_graph,coefficient*sign))
            for word,c in expanded.items():polynomial[word]+=coefficient*sign*c
        row=ct.pf.primitive({w:c for w,c in polynomial.items() if c})
        assert not ct.pf.spin_action(row,True)
        assert echelon.add(row)
        rows.append(row);recipes.append((kind,tuple(recipe)))
        record=dict(source_index=source_index,new_index=len(rows)-1,trace_terms=len(row),recipe_terms=len(recipe))
        records.append(record);print(json.dumps(record),flush=True)
    assert len(rows)==183
    output=HERE/'complete_target_sp5.pkl'
    output.write_bytes(pickle.dumps((rows,recipes),protocol=4))
    result=dict(rank=rank,source_rank=4,complete=True,relation_dimension=len(rows),
                additions=records,elapsed_seconds=time.monotonic()-start,
                cohomology_transport_assumed=False)
    (HERE/'complete_target_sp5_transport.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
