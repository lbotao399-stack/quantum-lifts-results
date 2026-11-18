"""Exact rational highest-weight differential at derivative bidegree (4,4)."""
import argparse
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
import json, sys
import numpy as np
from flint import fmpq_mat
from charge_maps import RUNTIME, dimensions
from g2_full_classical_screen_v1 import contents, sector
from sp_environment import canonical, fw
from full_extremal import project
from fortuity.full.wordenum import multitraces_with_content
from fortuity.hwcoh.evalg import table

MU=(4,4,0,0,0)
WEIGHT=(0,0,8,8,8)

def mat(rows, width):
    return fmpq_mat(len(rows),width,[str(x) for row in rows for x in row])

def pivots(matrix):
    r,rank=matrix.rref()
    return [next(j for j in range(r.ncols()) if r[i,j]) for i in range(rank)]

def array(expressions,words):
    return mat([[e.get(w,0) for w in words] for e in expressions],len(words))

def q(expr):
    return canonical(fw.Q_expr(expr))

def independent(exprs):
    words=sorted({w for e in exprs for w in e})
    idx=pivots(array(exprs,words).transpose())
    return [exprs[i] for i in idx]

def formal(length):
    seeds=sorted({v for letters in contents(MU,length) for w in multitraces_with_content(letters)
                  for v in canonical({w:1})})
    exprs=independent([project({w:1},sector(MU,length)) for w in seeds])
    return exprs

class IntegerEvaluation:
    """Each exterior coefficient is an int64; rational sums of traces use FLINT."""
    def __init__(self,group,length,letters,points):
        self.k=length
        n=7 if group=='SO_odd' else 6
        # Each individual trace product has at most k! n^k signed terms per coordinate.
        assert factorial(length)*n**length < 2**63
        rng=np.random.default_rng(2026111815+length)
        self.values={}
        for letter in sorted(letters):
            upper=np.triu(rng.integers(-1,2,(points,length,n,n),dtype=np.int64))
            if group=='SO_odd':
                m=upper-upper.swapaxes(-1,-2)
            else:
                sym=upper+upper.swapaxes(-1,-2)
                ix=np.arange(n); sym[...,ix,ix]=upper[...,ix,ix]
                m=np.concatenate((-sym[...,3:,:],sym[...,:3,:]),axis=-2)
            assert not letter[2] and np.abs(m).max()<=1
            self.values[letter]=m
    @lru_cache(None)
    def prefix(self,word):
        if len(word)==1: return self.values[word[0]]
        a=self.prefix(word[:-1]);b=self.values[word[-1]]
        i,j,s,st=table(self.k,len(word)-1,1)
        return np.add.reduceat((a[:,i]@b[:,j])*s[None,:,None,None],st,axis=1)
    @lru_cache(None)
    def trace(self,w):
        return self.prefix(w).trace(axis1=-2,axis2=-1)
    def word(self,w):
        value=self.trace(w[0]);deg=len(w[0])
        for trace in w[1:]:
            i,j,s,st=table(self.k,deg,len(trace))
            value=np.add.reduceat(value[:,i]*self.trace(trace)[:,j]*s[None,:],st,axis=1)
            deg+=len(trace)
        assert deg==self.k and value.shape[1]==1
        return value[:,0].tolist()
    def evaluate(self,exprs,width):
        words=sorted({w for e in exprs for w in e})
        return array(exprs,words)*mat([self.word(w) for w in words],width)
    def clear(self):
        self.trace.cache_clear();self.prefix.cache_clear()

def solve(group):
    rows=[]; prior=[]; maps=[]; S=(((0,1,()),(1,0,())),)
    for length in range(2,9):
        d=sum(dimensions(group,MU,length).values())
        points=d+4
        e=formal(length); im=[q(x) for x in prior]
        assert all(not q(x) for x in im)
        letters={l for x in e+im for w in x for t in w for l in t}
        frame=IntegerEvaluation(group,length,letters,points)
        ev=frame.evaluate(e,points)
        idx=pivots(ev.transpose()); assert len(idx)==d
        basis=[e[i] for i in idx]; base=frame.evaluate(basis,points)
        col=pivots(base); square=mat([[base[i,j] for j in col] for i in range(d)],d)
        det=square.det()
        assert det
        def coordinates(exprs):
            v=frame.evaluate(exprs,points)
            u=mat([[v[i,j] for j in col] for i in range(len(exprs))],d)*square.inv()
            assert u*base==v
            return u
        diff=coordinates(im)
        if maps: assert not any((maps[-1]*diff).entries())
        maps.append(diff)
        row=dict(length=length,dimension=d,formal_trace_dimension=len(e),incoming_rank=diff.rank(),
            eval_columns=col,basis_eval_determinant=str(det),differential=[[str(x) for x in line] for line in diff.tolist()])
        if length==8:
            s=canonical({S*4:1})
            assert not q(s)
            x=coordinates([s]); row['S4_class_rank_increment']=mat(diff.tolist()+x.tolist(),d).rank()-diff.rank()
            row['S4_coordinates']=[str(a) for a in x.entries()]
        frame.clear()
        prior=basis; rows.append(row)
    for i,r in enumerate(rows):
        r['outgoing_rank']=rows[i+1]['incoming_rank'] if i+1<len(rows) else 0
        r['cohomology']=r['dimension']-r['incoming_rank']-r['outgoing_rank']
        assert r['cohomology']>=0
    assert not contents(MU,9)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    groups={g:solve(g) for g in ('SO_odd','Sp')}
    for g,rows in groups.items():
        print(g,[(r['length'],r['dimension'],r['incoming_rank'],r['outgoing_rank'],r['cohomology']) for r in rows],rows[-1]['S4_class_rank_increment'],flush=True)
    a.output.write_text(json.dumps(dict(mu=MU,target_charges=WEIGHT,field='Q',groups=groups),indent=2)+'\n')
