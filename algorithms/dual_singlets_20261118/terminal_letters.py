"""Check the complete terminal and predecessor alphabets against the full-letter runtime."""
import argparse
import json
from itertools import combinations_with_replacement
from pathlib import Path
from charge_maps import RUNTIME
from g2_full_classical_screen_v1 import contents
from fortuity.full.letters import Q_letter, is_ferm, letters_up_to
from terminal_commutators import weight_terms


def degree(l): return l[0]+l[1]+len(l[2])


def compute():
    alphabet = letters_up_to(2, 2)
    simple = {x for x in alphabet if degree(x)==1}
    double = {x for x in alphabet if degree(x)==2}
    assert len(simple)==5 and len(double)==12 and all(not Q_letter(x) for x in simple)
    expected = {}
    for x,y in combinations_with_replacement(sorted(simple),2):
        if set(x[2]) & set(y[2]): continue
        l=x[0]+y[0],x[1]+y[1],tuple(sorted(x[2]+y[2])); assert l not in expected
        expected[l]=x,y
    assert set(expected)==double
    rules=[]
    for l,(x,y) in sorted(expected.items()):
        action=Q_letter(l); a={(v,w):c for c,v,w in action}; assert len(a)==len(action)
        if x==y: assert is_ferm(x) and a=={(x,x):1}
        else: assert set(a)=={(x,y),(y,x)} and abs(a[x,y])==1 and a[y,x]==-(-1)**(is_ferm(x)*is_ferm(y))*a[x,y]
        rules.append(dict(letter=l,action=action))
    rows=[]
    for h in range(5):
        L=12-h
        for w,sign in weight_terms(h):
            if min(w)<0:continue
            assert sum(w)==L
            top=contents(w,L); pre=contents(w,L-1)
            assert len(top)==1 and all(x in simple for x in top[0])
            assert all(len(c)==L-1 and sum(x in double for x in c)==1 and all(x in simple|double for x in c) for c in pre)
            feasible={l for l in double if all(v<=q for v,q in zip((l[0],l[1])+tuple(int(i in l[2]) for i in range(3)),w))}
            assert len(pre)==len(feasible) and feasible=={x for c in pre for x in c if degree(x)==2}
            rows.append(dict(bidegree=h,weight=w,sign=sign,terminal_length=L,predecessor_contents=len(pre)))
    return dict(weight_one_letters=sorted(simple),all_weight_two_rules=rules,dominant_and_shifted_content_checks=rows)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);a=ap.parse_args();assert __debug__;r=compute()
    if a.output:a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(len(r['all_weight_two_rules']),len(r['dominant_and_shifted_content_checks']),sum(e['predecessor_contents'] for e in r['dominant_and_shifted_content_checks']))
