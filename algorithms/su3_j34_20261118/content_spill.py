"""Which of all 500 first-order equations can an internal 12-content repair change?

Use all forward classical letter variations and separately reverse the admissible
merges. The resulting incidence is a support bound; predicted nonzero edges are
not assumed to survive trace identities. Counts outside that support are exact
for every correction to the archived even vector using only the named parents.
"""
import argparse
from collections import Counter
from pathlib import Path
import joint_content_correction as jc

s,l=jc.s,jc.l

def run(graph,probe):
    z=s.read(graph);p=s.read(probe);assert z['base']==p['vector_sha'] and z['rhs']==p['rhs_sha']
    domain=set(l.contents(11,0,l.CHARGES));base=set(l.contents(10,0,l.CHARGES));assert len(domain)==500 and len(base)==581
    F={g:{bytes(sorted(g[:i]+g[i+1:]+pair)) for i,k in enumerate(g) for _,pair in l.f.SPLITS[k]} for g in base}
    back={t:set(l.tri.out.incident(t)) for t in domain};assert all(t in domain for vals in F.values() for t in vals);assert {(g,t) for t,gg in back.items() for g in gg}=={(g,t) for g,ts in F.items() for t in ts}
    internal=set(map(bytes.fromhex,z['internal']));can=set.union(*(F[g] for g in internal));fixeds={bytes.fromhex(x) for x in p['exact_zero']};fixedf={bytes.fromhex(x) for x in p['mathematically_nonzero']};assert fixeds|fixedf==domain and not fixeds&fixedf
    x12=set(map(bytes.fromhex,z['targets']));x70=jc.side.replay.domain('joint');assert can&x70==x12
    spill=can-x70;locked_failed=fixedf-can;locked_ok=fixeds-can
    # Components with no nonzero RHS cannot require a cross-component repair.
    todo=set(domain);components=[]
    while todo:
        t=min(todo);ts={t};parents=set();front={t}
        while front:
            ins=set.union(*(back[v] for v in front))-parents;parents.update(ins);front=({t for g in ins for t in F[g]}-ts);ts.update(front)
        todo-=ts;components.append({'targets':len(ts),'parents':len(parents),'rhs':len(ts&fixedf),'min_target':min(ts).hex()})
    enc=lambda x:[t.hex() for t in sorted(x)]
    return dict(input_rhs=z['rhs'],base=z['base'],previous_exact_zero=len(fixeds),previous_nonzero=len(fixedf),parents=len(base),targets=len(domain),edges=sum(map(len,F.values())),all_components=components,internal_parents=len(internal),direct_targets=enc(x12),additional_possible_targets=enc(spill),locked_zero=enc(locked_ok),locked_failed=enc(locked_failed),old_failures_in_spill=enc(spill&fixedf),old_passes_in_spill=enc(spill&fixeds),conditional=dict(condition='Any internal 12-content repair which solves all 12',min_pass=len(locked_ok|x12),max_pass=len(domain-locked_failed),max_unresolved=len(spill),remaining_definite_failures=len(locked_failed)),scope='Upper support only; no edge is assumed to have nonzero physical image. Unresolved bounds concern this correction space, not arbitrary first corrections.')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--graph',type=Path,required=True);a.add_argument('--probe',type=Path,required=True);a.add_argument('--out',type=Path,required=True);v=a.parse_args();assert __debug__;r=run(v.graph,v.probe);s.write(v.out,r);print('edges, components, additional, fixed failures, fixed successes, conditional counts:',r['edges'],len(r['all_components']),len(r['additional_possible_targets']),len(r['locked_failed']),len(r['locked_zero']),r['conditional'])
