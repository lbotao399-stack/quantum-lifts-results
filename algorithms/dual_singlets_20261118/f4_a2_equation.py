"""The scalar correction in an A2 factor of the maximal A2+A2 subalgebra."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

import mixed_incoming as m
from f4_subgroup_equation import integral_rows
from fortuity.full.letters import is_ferm
from g2_full_classical_screen_v1 import dimensions
from mixed33_exact import apply, exact
from rank_four_sixth import G


def canonical(traces):
    sign=1
    cycles=[]
    for tr in traces:
        s,t=m.fw.canon_word(tr)
        if not s:
            return 0,None
        sign*=s
        cycles.append(t)
    s,t=m.fw.canon_multi(tuple(cycles))
    return sign*s,t


def raise_word(word,axis):
    ans=Counter()
    for i,tr in enumerate(word):
        for j,(a,b,theta) in enumerate(tr):
            if axis==3:
                if not b:
                    continue
                new=(a+1,b-1,theta)
                coef=a+1
            else:
                if axis not in theta or axis-1 in theta:
                    continue
                new=(a,b,tuple(sorted(axis-1 if t==axis else t for t in theta)))
                coef=1
            s,v=canonical(word[:i]+(tr[:j]+(new,)+tr[j+1:],)+word[i+1:])
            if s:
                ans[v]+=s*coef
    return {w:c for w,c in ans.items() if c}


def rows(weight,length,degree):
    blocks=defaultdict(set)
    for c in m.contents(weight,length):
        if sum(map(is_ferm,c))!=degree:
            continue
        for w in m.multitraces_with_content(c):
            blocks[m.shape(w,"spin")].add(w)
    out=[]
    for _,block in sorted(blocks.items()):
        words=sorted(block)
        equations=[{(a,v):c for a in (1,2,3) for v,c in raise_word(w,a).items()} for w in words]
        ids={x:i for i,x in enumerate(sorted({k for row in equations for k in row}))}
        matrix=[{ids[k]:n for k,n in row.items()} for row in equations]
        _,_,kernel=exact(matrix)
        for k in kernel:
            assert not apply(k,matrix)
            out.append({words[i]:c for i,c in k.items()})
    return out


def q(expr):
    result=Counter()
    for w,c in expr.items():
        for v,n in m.fw.Q_multi(w).items():
            result[v]+=c*n
    return {v:c for v,c in result.items() if c}


def run(workers):
    w=(2,2,1,1,1)
    src=rows(w,6,3); end=rows(w,7,4)
    ds,dt=dimensions(w,6,"SU3")[3],dimensions(w,7,"SU3")[4]
    image=[q(r) for r in src]
    a=((0,0,(0,)),(0,0,(1,)),(0,0,(2,)))
    rhs=Counter()
    for tr,weight in ((a,1),((a[0],a[2],a[1]),-1)):
        s,v=canonical((G,G,tr))
        rhs[v]+=s*weight
    assert not q(rhs)
    sp,sf=integral_rows(src,"SU3",3,ds,540909,workers)
    mat,tf=integral_rows(end+image+[rhs],"SU3",4,dt,921873,workers)
    x=mat.tolist()
    bd=fmpz_mat(x[len(end):-1]);aug=fmpz_mat(x[len(end):])
    ranks=[sp.rank(),fmpz_mat(x[:len(end)]).rank(),bd.rank(),aug.rank()]
    assert ranks[:2]==[ds,dt]
    result=dict(group="SU3",power=2,weight=list(w),source_character=ds,target_character=dt,
                source_frame=sf,target_frame=tf,ranks=ranks,
                rhs=[[m.encode(w),c] for w,c in sorted(rhs.items())])
    if ranks[2]<ranks[3]:
        from f4_pure_target import covector
        result["separating_covector"]=covector(aug)
        result["membership"]="excluded over Q"
    else:
        red,rank=fmpq_mat(aug).transpose().rref()
        co={}
        for i in range(rank):
            column=next(j for j in range(red.ncols()) if red[i,j])
            assert column<len(src)
            if red[i,len(src)]:
                co[column]=red[i,len(src)]
        assert (fmpq_mat([[co.get(j,0) for j in range(len(src))]])*fmpq_mat(bd)).tolist()[0]==x[-1]
        result["primitive_in_source_rows"]=[[k,str(c)] for k,c in co.items()]
        result["membership"]="included over Q"
    print("A2",ds,dt,ranks,result["membership"],flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    assert __debug__
    m.write(a.out,run(a.workers))
