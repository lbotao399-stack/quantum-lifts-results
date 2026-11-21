"""Test exact order-two nilpotency on a nonzero three-colour target coefficient."""
import argparse,json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from fortuity.full.letters import Q_letter,is_ferm
from replay_third_action import basis_values,load
from third_consistency import calc,shape
from third_action import fw,HERE
from verify import add_scaled

def possible(g):
    found=set()
    for t in set(g):
        rest=list(g);rest.remove(t)
        for c,L,R in Q_letter(t):
            target=tuple(sorted(rest+[L,R]))
            if sum(bool(a or b) for a,b,ss in target)==3 and sum(map(is_ferm,target))<=5:found.add(target)
    return found

def run(out, level32=False):
    if level32:
        location = HERE.parent / "dual_singlets_20261118/taylor_hooks/su3_level32"
        first_file, first_sign = location / "Q1.jsonl.gz", 1
        first_sha, second_sha = (json.loads((location / f"order{i}.json").read_bytes())["output_sha256"]
                                 for i in (1,2))
    else:
        location = HERE
        first_file, first_sign = HERE / "negative_Q1.jsonl.gz", -1
        first_sha = json.loads((HERE / "first_action.json").read_bytes())["rhs_uncompressed_sha256"]
        second_sha = json.loads((HERE / "second_action.json").read_bytes())["output_jsonl_sha256"]
    d,digest=load(location/'Q2.jsonl.gz');assert digest==second_sha;counts=Counter(map(shape,d));forward={g:possible(g) for g in counts};cost=Counter()
    for g,n in counts.items():
        for t in forward[g]:cost[t]+=n
    tested=[];found=None
    for target,n in sorted(cost.items(),key=lambda item:(item[1],item[0])):
        src={g for g in counts if target in forward[g]};row={};retained=0
        for w,c in d.items():
            if shape(w) not in src:continue
            retained+=1;add_scaled(row,{v:k for v,k in fw.Q_multi(w).items() if shape(v)==target},c)
        assert retained==n;ev=basis_values({target:row},346712) if row else None;val=Fraction(ev['contents'][0]['coefficient']) if ev else Fraction(0);tested.append(dict(target=target,incoming_Q2_image_words=retained,formal_terms=len(row),coefficient=str(val)))
        if val:found=(target,retained,row,ev,val);break
    assert found is not None;target,retained,row,ev,val=found
    # No field in Q0 O0 carries a derivative. Q2 inserts a complementary
    # pair of derivative jets on its two outputs. No term can have three.
    assert sum(bool(a or b) for a,b,S in target)==3
    other=calc(first_file,target,1,first_sign);assert other['input_sha256']==first_sha
    assert val+Fraction(other['coefficient'])==0
    doc=dict(Q2_input_sha256=digest,target=target,number_of_targets_examined=len(tested),examined=tested,Q0_Q2_input_words=retained,Q0_Q2_formal_words=len(row),Q0_Q2=str(val),Q1_Q1=other,Q2_Q0='0',sum='0',specialization=ev,scope='One exact finite-colour coefficient. Direct Q2 O0 is not Q0-closed. This neither settles the first correction nor establishes a nonzero d2.')
    out.write_text(json.dumps(doc,indent=2)+'\n');print('Q0 Q2, Q1 Q1, Q2 Q0:',val,other['coefficient'],0,'candidates tested',len(tested),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--out',type=Path,required=True);a.add_argument('--level32',action='store_true');assert __debug__;args=a.parse_args();run(args.out,args.level32)
