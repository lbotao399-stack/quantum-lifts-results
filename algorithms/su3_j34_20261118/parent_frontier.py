"""Determine which parent contents can affect each failed equation without touching a success.

Forward splits and reverse merges are checked independently. These are upper
bounds on support, never assertions that a column has nonzero physical image.
"""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import joint_content_correction as j

s,l=j.s,j.l

def incidence():
    parents=l.contents(10,0,l.CHARGES); targets=set(l.contents(11,0,l.CHARGES))
    forward={}; back={t:{} for t in targets}
    for g in parents:
        possible=defaultdict(set)
        for i,k in enumerate(g):
            for c,pair in l.f.SPLITS[k]:
                t=bytes(sorted(g[:i]+g[i+1:]+pair)); possible[t].add((k,c,*pair))
        forward[g]=set(possible)
        for t,terms in possible.items():
            (k,c,a,b),(k2,v,b2,a2)=sorted(terms);assert k==k2 and (a,b)==(a2,b2) and a!=b
            assert v==-c*(-1)**(int(l.f.ODD[a])*int(l.f.ODD[b]))
            back[t][g]=(a,b)
    assert (len(parents),len(targets))==(581,500)
    assert all(set(back[t])==l.tri.out.incident(t) for t in targets)
    return forward,back

def frontier(probe,old,spill):
    F,B=incidence(); p=s.read(probe);o=s.read(old);w=s.read(spill)
    assert not p['unresolved'] and not o['unresolved'] and w['input_rhs']==p['rhs_sha']==o['rhs_sha'] and w['base']==o['vector_sha']
    dec=lambda seq:{bytes.fromhex(t) for t in seq}; enc=lambda seq:[g.hex() for g in sorted(seq)]
    ok=dec(p['exact_zero']);bad=dec(p['mathematically_nonzero']);before=dec(o['exact_zero']);assert ok.isdisjoint(bad) and ok|bad==set(B)
    locked=dec(w['locked_failed']);before_locked=dec(w['locked_zero']);assert locked<=bad and before_locked<=ok
    safe={g for g,v in F.items() if v and not v&ok}; reachable={t for g in safe for t in F[g]}; other={g for g in F if len(F[g]&ok)==1}
    changes={t.hex():{'all':enc(B[t]),'confined_to_failures':enc(set(B[t])&safe),'touch_one_success':enc(set(B[t])&other)} for t in sorted(bad)}
    rows={g.hex():enc(F[g]) for g in sorted(safe)}
    return dict(vector_sha=p['vector_sha'],rhs=p['rhs_sha'],compared_with=o['vector_sha'],unchanged_locked_failures=enc(locked),unchanged_locked_successes=enc(before_locked),gained=enc(ok-before),lost=enc(before-ok),successes=len(ok),failures=len(bad),letter_edges=sum(map(len,F.values())),even_parent_dimension=sum(j.E.parity_frame_dims(g)[0] for g in F),even_target_dimension=sum(j.E.parity_frame_dims(t)[0] for t in B),safe_parents=rows,failed_targets_without_safe_parent=enc(bad-reachable),locked_targets_without_safe_parent=enc(locked-reachable),reachable_locked=len(locked&reachable),target_parents=changes,failed_by_diagonal_count=dict(sorted(Counter(t.count(0) for t in bad).items())),scope='Changes by a sum of parents individually confined to failing targets preserve every current success. Other sums may also preserve successes by cancellation. Dimensions are fixed weight before the flavour highest-weight quotient.')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--probe',type=Path,required=True);a.add_argument('--base-probe',type=Path,required=True);a.add_argument('--spill',type=Path,required=True);a.add_argument('--out',type=Path,required=True);v=a.parse_args();assert __debug__;z=frontier(v.probe,v.base_probe,v.spill);s.write(v.out,z);print('single-parent supports:',len(z['safe_parents']), 'locked reachable',z['reachable_locked'], 'all failures without one', len(z['failed_targets_without_safe_parent']))
