"""Identify the finite-colour cycles using nonzero minors of the stored boundary rows.

Only the 225 previously certified boundary directions are read; no singlets are
regenerated. The F row is evaluated as published (denominator cleared by 252).
The determinant of a full-row-rank prefix also certifies every preceding prefix.
"""
import argparse
from collections import Counter
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import mixed_incoming as m
from mixed33_exact import apply,read,sha


def inputs():
    p=m.HERE/'mixed33';ip=m.HERE/'mixed33_incoming.json.gz';op=p/'outgoing.json.gz';fp=m.HERE.parent/'so7_cartan_target_20261118/representatives.json.gz';cp=p/'formal_certificate.json.gz'
    inc,out,cert,F=map(read,(ip,op,cp,fp));assert sha(ip)==cert['incoming_plan'] and sha(op)==cert['outgoing_plan']
    iw=list(map(m.decode,inc['target_dictionary']));ow=list(map(m.decode,out['source_dictionary']));src=list(map(dict,out['sources']));qr=list(map(dict,out['images']))
    boundary=[{iw[i]:c for i,c in inc['images'][j]} for j in cert['incoming_basis_rows']];assert len(boundary)==225
    cycles=[]
    for recipe in cert['formal_class_recipes']:
        assert not apply(dict(recipe),qr);cycles.append({ow[i]:c for i,c in apply(dict(recipe),src).items()})
    assert len(cycles)==2
    letters=[(i,j,tuple(t)) for i,j,t in F['letters']];den=lcm(*(d for _,d,_ in F['F']));z=Counter()
    for num,d,word in F['F']:
        w=tuple(tuple(letters[l] for l in t) for t in word);s,v=m.so_canon(w);z[v]+=s*num*(den//d)
    g={w:c for w,c in z.items() if c};assert den==252 and len(g)==2969 and len(m.q(g))==1279
    for k in (1,2,3):
        v=Counter();
        for w,c in g.items():
            for t,n in m.raise_word(w,k).items():v[t]+=c*n
        assert not any(v.values())
    return boundary,cycles,g,dict(formal_certificate_sha256=sha(cp),incoming_sha256=sha(ip),outgoing_sha256=sha(op),F_sha256=sha(fp),F_denominator=den,class_recipes=cert['formal_class_recipes'],cycle_words=list(map(len,cycles)),F_terms=len(g),formal_Q0F_terms=1279)


def run(group,p,workers):
    boundary,cycles,F,meta=inputs();states=boundary+cycles+([F] if group=='SO_odd' else [])
    ws,rows=m.flatten(states);plan=dict(target_dictionary=[m.encode(w) for w in ws],images=rows,source_words=len(states))
    with TemporaryDirectory() as d:r=m.run(plan,group,p,Path(d),workers,[6,8],8,5658813 if p==1000003 else 6351915,12)
    n=len(states);assert r['rank']==n and r['rows']==list(range(n))
    return dict(**meta,group=group,prime=p,rank_prefixes=[[225,225],[227,227]]+([[228,228]] if n==228 else []),matrix_shape=r['matrix_shape'],trace_words=len(ws),evaluation_sha256=r['evaluation_sha256'],columns=r['columns'],determinant=r['determinant'],minor_sha256=r['minor_sha256'],frames=r['frames'],native_source_sha256=m.digest((m.RT/'trace_frame_bmn_rank_threadsafe_v1.cpp').read_bytes()))

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--group',choices=('Sp','SO_odd'),required=True);a.add_argument('--prime',type=int,choices=(1000003,1000033),required=True);a.add_argument('--workers',type=int,default=4);a.add_argument('--output',required=True,type=Path);x=a.parse_args();assert __debug__
    r=run(x.group,x.prime,x.workers);m.write(x.output,r);print(x.group,x.prime,r['rank_prefixes'],r['determinant'])
