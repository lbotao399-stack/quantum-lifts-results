"""Directly verify a fixed rational correction from the original -Q1 input and
canonical Q_multi, without using the solver's projected matrices or pivots.

The polynomial expander uses integer trace-word columns. All rational residual
coefficients are collected with unbounded Python integers. A nonzero integer
monomial certifies failure of this particular correction on an optional holdout.
"""
import argparse,hashlib,json,subprocess
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from functools import lru_cache
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import replay_formal_blocks as s
from fortuity.full.words import Q_multi
from fortuity.full.letters import is_ferm
import scalar_filtration as lay

def content(w):return bytes(sorted(b''.join(w)))

def domain(mode):
    a=set(lay.contents(11,0,lay.CHARGES));assert len(a)==500
    high={c for c in a if c.count(0)>=6};assert len(high)==67
    return (set(lay.tri.CONTENTS) if mode=='three' else high if mode=='top' else high|set(lay.tri.CONTENTS))

def replay(rhs,co,scope,exe,work,jobs,holdouts=(),selected_contents=None):
    selected=domain(scope);
    if selected_contents is not None:
        declared=set(selected_contents);assert declared and declared<=selected;selected=declared
    extras=set(holdouts);assert not extras&selected;test=extras|selected
    assert all(len(c)==11 and c.count(0)>=4 and c in lay.contents(11,0,lay.CHARGES) for c in extras);hd=hashlib.sha256();blocks={k:{} for k in test};n=0;rhs_count=Counter()
    for w,coef in s.stream_rows(rhs,hd):
        n+=1;k=content(w)
        if k in test:assert w not in blocks[k];blocks[k][w]=-Fraction(coef);rhs_count[k]+=1
    assert n==517468 and hd.hexdigest()==s.INPUT;hits=Counter();raw=s.read(co);src=raw['solution'] if isinstance(raw,dict) else raw;seen=set()
    for r,coef in src:
        w=s.unpack_hex(r);coef=Fraction(coef);assert w not in seen and coef and len(b''.join(w))==10;seen.add(w)
        for u,m in Q_multi(s.expanded(w)).items():
            v=s.packed(u);k=content(v)
            if k not in test:continue
            a=blocks[k];hits[k]+=1;a[v]=a.get(v,0)+coef*m
            if not a[v]:del a[v]
    odd=bytes(i for i,t in enumerate(s.ALPH) if is_ferm(t)).hex();rows=sorted((k,v) for k,v in blocks.items() if v);exe=str(exe.resolve());work.mkdir(parents=True,exist_ok=True)
    def expand(row):
        k,r=row;words=sorted(r);den=lcm(*(v.denominator for v in r.values()));num=[int(r[w]*den) for w in words];
        with TemporaryDirectory(dir=work) as d:
            d=Path(d);inp=d/'words';out=d/'polys';inp.write_text(str(len(words))+'\n'+''.join(f"{len(w)} {' '.join(s.hexes(w))}\n" for w in words));subprocess.run([exe,str(inp),odd,str(out)],check=True)
            z=Counter();count=0
            with out.open() as lines:
                for count,line in enumerate(lines,1):
                    a=list(map(int,line.split()));assert 1+2*a[0]==len(a);c=num[count-1]
                    for j in range(1,len(a),2):z[a[j]]+=c*a[j+1]
            assert count==len(words)
        nz={j:v for j,v in z.items() if v};first=min(nz) if nz else None
        return dict(content=k.hex(),formal_residual_terms=len(r),monomials_collected=len(z),nonzero_monomials=len(nz),first_key=first,first_coefficient=str(Fraction(nz[first],den)) if nz else None)
    with ThreadPoolExecutor(jobs) as pool:detail={bytes.fromhex(v['content']):v for v in pool.map(expand,rows)}
    def fields(k):
        row=detail.get(k,dict(content=k.hex(),formal_residual_terms=0,monomials_collected=0,nonzero_monomials=0,first_key=None,first_coefficient=None));return dict(rhs_terms=rhs_count[k],nonzero_parent_word_incidences=hits[k],**row)
    result=dict(rhs_uncompressed_sha256=s.INPUT,ordered_source_sha256=hashlib.sha256(json.dumps(src,separators=(',',':')).encode()).hexdigest(),source_terms=len(src),scope=scope if selected_contents is None or selected==domain(scope) else scope+'-subset',complete_contents=len(selected),total_rhs_terms=sum(rhs_count[k] for k in selected),all_input_rhs_terms=n,contents=[fields(k) for k in sorted(selected)],holdouts=[fields(k) for k in sorted(extras)])
    assert not any(r['nonzero_monomials'] for r in result['contents']);return result
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--scope',choices=('three','top','joint'),required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--vector',type=Path,required=True);a.add_argument('--exe',type=Path,required=True);a.add_argument('--work',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--jobs',type=int,default=2);a.add_argument('--select',type=Path,help='Optional JSON with the exact enforced content list for a subset replay.');a.add_argument('--holdouts',nargs='*',default=[]);b=a.parse_args();assert __debug__;res=replay(b.rhs,b.vector,b.scope,b.exe,b.work,b.jobs,map(bytes.fromhex,b.holdouts),None if not b.select else map(bytes.fromhex,s.read(b.select)['enforced']));s.write(b.out,res);print('Source, contents, input RHS terms; formal residual terms, collected monomials, nonzero monomials:',res['source_terms'],res['complete_contents'],res['total_rhs_terms'],*(sum(x[k] for x in res['contents']) for k in ('formal_residual_terms','monomials_collected','nonzero_monomials')),flush=True);print('holdout',[(x['content'],x['nonzero_monomials'],x['first_coefficient']) for x in res['holdouts']],flush=True)
