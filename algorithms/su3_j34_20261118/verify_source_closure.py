"""Expand the published SU(3) source's classical differential at finite colour."""
import argparse,hashlib,json,subprocess
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from itertools import product
from math import lcm
from pathlib import Path
from second_action import read_source,HERE,BMN_LETTERS,SOURCE_SHA
from verify import add_scaled
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm
def prepare(source):
    src=read_source(source);raw={};inputs=defaultdict(dict);labels={L:i for i,L in enumerate(sorted(BMN_LETTERS))};odd=bytes(i for L,i in labels.items() if is_ferm(L)).hex()
    for w,c in src.items():add_scaled(raw,fw.Q_multi(w),c)
    for w,c in raw.items():
        word=tuple(bytes(labels[L] for L in tr) for tr in w);g=bytes(sorted(b''.join(word)));inputs[g][word]=c
    assert len(raw)==25299 and len(src)==6270 and len(inputs)==15
    return labels,odd,inputs
def write_words(path,group,scalar,odd):
    D=lcm(*(c.denominator for c in group.values()))
    with path.open('w') as st:
        st.write(f'{scalar} {odd} {len(group)}\n')
        for w,c in sorted(group.items()):st.write(f'{int(c*D)} {len(w)} '+' '.join(t.hex() for t in w)+'\n')
    return D
def run(source,exe,work,out,workers,limit=None):
    work.mkdir(exist_ok=True,parents=True);labels,odd,inputs=prepare(source);evens={i for L,i in labels.items() if not is_ferm(L)}
    def one(t):
        row=inputs[t];c=Counter(t);sc=min(evens,key=lambda i:(-c[i],i));g=t.hex();path=work/(g+'.txt');target=work/(g+'.json');D=write_words(path,row,sc,odd);subprocess.run([str(exe.resolve()),str(path),str(target)],check=True);v=json.loads(target.read_bytes());assert v['formal_words']==len(row) and v['nonzero_monomials']==0
        # A small positive control from the same content prevents vacuous zero output.
        control=work/(g+'.control.txt');ans=work/(g+'.control.json');checked=0
        for w,coef in sorted(row.items()):
            write_words(control,{w:coef},sc,odd);subprocess.run([str(exe.resolve()),str(control),str(ans)],check=True);a=json.loads(ans.read_bytes());checked+=1
            if a['nonzero_monomials']:break
        assert a['nonzero_monomials']>0;v.update(denominator=str(D),diagonal_letter=sc,diagonal_multiplicity=c[sc],control_words_tried=checked,control_word=[z.hex() for z in w],control_nonzero_monomials=a['nonzero_monomials'],control_first_integer_coefficient=a['first_integer_coefficient']);print(g,len(row),v['distinct_monomials'],flush=True);return g,v
    ordered=sorted(inputs,key=lambda g:(len(inputs[g]),g));chosen=ordered[:limit] if limit is not None else ordered
    with ThreadPoolExecutor(workers) as pool:records=dict(pool.map(one,chosen))
    result=dict(source_sha256=SOURCE_SHA,letter_table=sorted(labels,key=labels.get),formal_derivative_words=25299,all_contents=len(inputs),tested_contents=len(records),tested_words=sum(z['formal_words'] for z in records.values()),monomials=sum(z['distinct_monomials'] for z in records.values()),nonzero_monomials=sum(z['nonzero_monomials'] for z in records.values()),complete=set(records)==set(t.hex() for t in inputs),contents=records)
    out.write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=HERE/'O34_C_even.json.gz');p.add_argument('--exe',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=1);p.add_argument('--limit',type=int);assert __debug__;print(run(**vars(p.parse_args()))['monomials'])
