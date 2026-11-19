"""Evaluate one exterior coefficient of the rational RHS exactly and verify every allowed classical replacement.

Only Python integers enter the trace calculation; the saved matrices use the first three odd generators.
"""
import argparse,gzip,hashlib,json,math,sys
from fractions import Fraction
from functools import lru_cache
from collections import Counter
from pathlib import Path
import numpy as np
import projected_q0 as maps
import edge_annihilator as witness
from coupled_first_replay import RHS_SHA,PRIME

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));import sp_environment
    from fortuity.full.letters import is_ferm,Q_letter
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents
    witness.Frame,witness.Q_letter,witness.is_ferm=Frame,Q_letter,is_ferm
    plan=json.loads((a.root/'classical_order1.json').read_bytes())['incoming'];mu=tuple(plan['mu']);t=plan['targets'][6];pat=maps.pattern(t);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});allrec=json.loads((a.root/'joint/edge_annihilator.json').read_bytes());rec=next(r for r in allrec['records'] if r['target']==6 and r['frame']=='diagonal_scalars');
    expected=[3,4,7,8,9]; point_index=0
    if a.mask is not None:
        from edge_masks import make
        records=json.loads((a.root/'joint/edge_masks.json').read_bytes());assert records['rhs_stream_sha256']==RHS_SHA
        rec=next(x for x in records['records'] if x['mask']==a.mask);expected=rec['possible'];idx=next(i for i,x in enumerate(rec['values']) if x);point_index,component=divmod(idx,10);assert component==0
        rec={**rec,'rhs_first_nonzero':[0,rec['values'][idx]],'matrices_sha256':rec['matrix_sha256']}
    if a.matrix:
        raw=json.loads(a.matrix.read_bytes());v={(x,y,tuple(z)):np.asarray(m,dtype=object) for (x,y,z),m in raw['entries']};assert set(v)==set(alpha)
    else:
        _,points=(witness.int_frame(alpha,3,rec['frame'],rec['points'],rec['seed']) if a.mask is None else make(alpha,a.mask,rec['points'],rec['seed']));assert hashlib.sha256(b''.join(points[l].astype('<i8').tobytes() for l in alpha)).hexdigest()==rec['matrices_sha256']
        v={l:points[l][point_index,:(3 if is_ferm(l) else 1)].astype(object) for l in alpha};raw=dict(entries=[[l,v[l].tolist()] for l in alpha],exterior_coefficient=[0,1,2]);a.out.mkdir(parents=True,exist_ok=True);(a.out/'point.json').write_text(json.dumps(raw,separators=(',',':'))+'\n')
    assert raw['exterior_coefficient']==[0,1,2];assert all(np.trace(m)==0 for val in v.values() for m in val)
    def sign(x,y):
        if x&y:return 0
        return (-1)**sum(1 for i in range(3) for j in range(i) if x>>i&1 and y>>j&1)
    blocks={l:({1<<i:m for i,m in enumerate(ms)} if is_ferm(l) else {0:ms[0]}) for l,ms in v.items()}
    def multiply(x,y,matrix):
        result={}
        for i,one in x.items():
            for j,two in y.items():
                sg=sign(i,j)
                if sg:result[i|j]=result.get(i|j,0)+sg*((one@two) if matrix else one*two)
        return result
    rules={};allowed=[]
    for s in plan['sources']:
        if 6 not in s['targets']:continue
        arr=[];groups=witness.matching(alpha,pat,maps.pattern(s));assert groups
        for l,rr in groups:
            expr={}
            for c,x,y in rr:
                for k,b in multiply(blocks[x],blocks[y],True).items():expr[k]=expr.get(k,0)+c*b
            nz=sum(np.count_nonzero(z) for z in expr.values());arr.append([l,int(nz)])
        rules[str(s['index'])]=arr
        if any(n for _,n in arr):allowed.append(s['index'])
    assert set(allowed)<=set(expected)
    @lru_cache(None)
    def trace(w):
        value=blocks[w[0]]
        for l in w[1:]:value=multiply(value,blocks[l],True)
        return {k:int(np.trace(b)) for k,b in value.items() if np.trace(b)}
    h=hashlib.sha256();numer=terms=total=0
    with gzip.open(a.root/'negative_Q1.jsonl.gz','rb') as stream:
        for line in stream:
            h.update(line);total+=1;c,ww=json.loads(line);w=maps.decode(ww)
            if witness.pattern(l for tr in w for l in tr)!=pat:continue
            terms+=1;c=Fraction(c)*120;assert c.denominator==1
            val={0:1}
            for tr in w:
                val=multiply(val,trace(tr),False)
                if not val:break
            numer+=int(c)*val.get(7,0)
    assert total==517468 and terms==31720 and h.hexdigest()==RHS_SHA;value=Fraction(numer,120);mod=value.numerator*pow(value.denominator,-1,PRIME)%PRIME;assert value and [0,mod]==rec['rhs_first_nonzero']
    result=dict(target=6,point_sha256=hashlib.sha256((json.dumps(raw,separators=(',',':'))+'\n').encode()).hexdigest(),rhs_sha256=h.hexdigest(),target_terms=terms,trace_cache=trace.cache_info().currsize,exterior_coefficient=[0,1,2],rhs_numerator_at_common_denominator_120=str(numer),exact_pairing=str(value),modulus=PRIME,residue=mod,possible_source_profiles=allowed,rule_images_nonzero_entries=rules);a.out.mkdir(parents=True,exist_ok=True);(a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rule_images_nonzero_entries'}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=maps.HERE);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--out',type=Path,required=True);p.add_argument('--matrix',type=Path);p.add_argument('--mask',type=int,choices=range(32));run(p.parse_args())
