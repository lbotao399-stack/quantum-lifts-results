"""Evaluate all degree-three anticommutator terms on one exact exterior coefficient.

For Q1 after Q2 and Q2 after Q1, reverse the target letter counts before
generating a trace. Each contraction deletes respectively two / three letters
and inserts two.  Keep the full three-colour reconnection of every candidate.
The saved Q1 input is NEGATIVE; change its sign before applying Q2.
"""
import argparse,gzip,hashlib,json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations,permutations
from pathlib import Path
from time import monotonic
from replay_third_action import grouped,basis_values,load
from third_action import HERE,fw,read_source,Q3d_multi
from verify import add_scaled
import su_quantum_integer_v1 as qi
def letter(x):return x[0],x[1],tuple(x[2])
def shape(w):return tuple(sorted(L for t in w for L in t))
def weigh(d,selection):
    n=sum(1 for _ in d);assert not n or set(map(shape,d))=={selection}
    if not d:return '0'
    return basis_values({selection:d},346712)['contents'][0]['coefficient']
def candidates(g,target,order):
    n=order+1;a=Counter(g);t=Counter(target);a_minus=tuple(sorted((a-t).elements()));t_minus=tuple(sorted((t-a).elements()))
    # Count differences determine how many contracted species can also reappear.
    if len(a_minus)>n or len(t_minus)>2:return False,0
    hits=0
    for removed in set(combinations(g,n)):
        unchanged=a-Counter(removed)
        if unchanged-t:continue
        inserted=tuple(sorted((t-unchanged).elements()));assert len(inserted)==2
        for seq in set(permutations(removed)):
            if any(tuple(sorted(z))==inserted for z in rules(seq,order)):hits+=1
    return bool(hits),hits
@lru_cache(None)
def rules(seq,order):return qi.species_rules(seq,order)
def calc(p,target,order,sign,tally_only=False):
    docs={};keep={};hd=hashlib.sha256();n=0;since=monotonic()
    with gzip.open(p,'rb') as st:
        for line in st:
            hd.update(line);c,q=json.loads(line);w=tuple(tuple(map(letter,tr)) for tr in q);g=shape(w);n+=1
            if g not in docs:docs[g]=candidates(g,target,order)
            if docs[g][0]:keep[w]=sign*Fraction(c)
    print(p.name,'order',order,'candidate contents',sum(d[0] for d in docs.values()),'words',len(keep),'of',n,flush=True)
    record=dict(input_sha256=hd.hexdigest(),input_words=n,content_count=len(docs),retained_content_count=sum(d[0] for d in docs.values()),retained_words=len(keep))
    if tally_only:return record
    result={};ordered=sorted(keep)
    for pos in range(0,len(ordered),64):
        selected=ordered[pos:pos+64];sub={w:keep[w] for w in selected};z=qi.act(sub,order,3);add_scaled(result,{w:c for w,c in z.items() if shape(w)==target},1)
        if pos%8192==0:print(order,pos,'of',len(keep),'elapsed',int(monotonic()-since),flush=True)
    record.update(selected_output_words=len(result),coefficient=weigh(result,target));return record
def run(replay,out,tally_only=False):
    rec=json.loads(replay.read_bytes());test=rec['Q0_of_direct_Q3_three_colour_test'];u=test['contents'][0];assert test['seed']==346712 and Fraction(u['coefficient']);target=tuple(map(letter,u['content']));a=calc(HERE/'Q2.jsonl.gz',target,1,1,tally_only);assert a['input_sha256']==json.loads((HERE/'second_action.json').read_bytes())['output_jsonl_sha256']
    b=calc(HERE/'negative_Q1.jsonl.gz',target,2,-1,tally_only);assert b['input_sha256']==json.loads((HERE/'first_action.json').read_bytes())['rhs_uncompressed_sha256']
    record=dict(bare_Q0_Q3=u['coefficient'],target=target,Q1_Q2=a,Q2_Q1=b,seed=346712)
    if not tally_only:
        # For Q3 Q0, the output letters carry total derivative (3,3).
        # Pure BMN Q3 always inserts jets (k,3-k) and (3-k,k), and leaves all
        # other derivatives zero. The selected target instead has three jet letters.
        assert sum(bool(i or j) for i,j,S in target)==3;record['Q3_Q0']="0";parts=[Fraction(u['coefficient']),Fraction(a['coefficient']),Fraction(b['coefficient'])]
        assert sum(parts)==0;record.update(exact_sum=str(sum(parts)),Q3_Q0_selection='Underived input to Q3 has exactly two derivative-carrying output letters. This target has three.');print([str(p) for p in parts],str(sum(parts)),flush=True)
    out.write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);h=HERE/'direct_third';p.add_argument('--replay',type=Path,default=h/'independent.json');p.add_argument('--out',type=Path,required=True);p.add_argument('--tally-only',action='store_true');assert __debug__;run(**vars(p.parse_args()))
