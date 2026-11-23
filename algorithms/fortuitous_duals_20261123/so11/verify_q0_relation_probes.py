#!/usr/bin/env python3
"""Replay full high-factor Q0 identities, then the two requested coordinates."""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import hashlib,json,subprocess,time
import independent_trace_ops as tr
from probe_q0_relations import source,q_inputs,run

HERE=Path(__file__).resolve().parent

def parse(raw):
    p=defaultdict(int)
    for line in raw.splitlines():
        w,c=line.split()
        if w=='!':continue
        s,k=tr.canonical(tuple(w.rstrip('.').split('.')))
        if s:p[k]+=s*int(c)
    return {w:c for w,c in p.items() if c}

def differentiated(p):
    out=defaultdict(int)
    for words,c in p.items():
        assert sum(t.count('D') for t in words)==1
        for i,t in enumerate(words):
            if 'D' not in t:continue
            for repl,sign in [('BC',1),('CB',-1)]:
                w=words[:i]+(t.replace('D',repl),)+words[i+1:]
                s,k=tr.canonical(w)
                if s:out[k]+=c*s*sign
    return {w:c for w,c in out.items() if c}

def main():
    start=time.monotonic();records=[]
    for index,expected in [(1,['4096','0']),(5,['4608','-1152'])]:
        spec=source(index)
        inp='6 0\n'+' '.join(str('ABCx'.index(c)) for c in spec['species'])+'\n'
        inp+=' '.join(['0']*13+['1'])+'\n'+' '.join(map(str,spec['endpoints']))+'\n'
        raw=subprocess.check_output([str(HERE/'probe_trials/bmngram14')],input=inp,text=True)
        pred=parse(raw);image=differentiated(pred)
        split=defaultdict(int)
        for sign,qinp in q_inputs(spec,6):
            qraw=subprocess.check_output([str(HERE/'probe_trials/bmngram')],input=qinp,text=True)
            for w,c in parse(qraw).items():split[w]+=sign*c
        split={w:c for w,c in split.items() if c}
        assert image==split,(index,len(image),len(split))
        r=run(spec);assert r['query_values']==expected
        r.update(source_canonical_terms=len(pred),derivative_canonical_terms=len(image),
                 all_retained_Q0_coordinates_match=True,
                 source_raw_sha256=hashlib.sha256(raw.encode()).hexdigest(),source_input=inp)
        records.append(r)
    det=Fraction(4096)*(-1152)-Fraction(4608)*0;assert det==-4718592
    assert [28575*int(records[0]['query_values'][i])-8000*int(records[1]['query_values'][i]) for i in range(2)]==[80179200,9216000]
    result={'passed':True,'records':records,'two_probe_image_rank':2,'minor_determinant':str(det),
            'closed_combination_probe_preimage':[28575,-8000],
            'matched_closed_combination_probes':[80179200,9216000],
            'scope':'Actual full Gram identities define predecessors in I. All retained high-factor coordinates of their Q0 images are independently replayed. The two chosen coordinate functionals have no nonzero linear combination annihilating full Q0(I). The two-probe preimage is not a preimage of the full closed-combination derivative.',
            'seconds':time.monotonic()-start}
    (HERE/'verified_q0_relation_probes.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
