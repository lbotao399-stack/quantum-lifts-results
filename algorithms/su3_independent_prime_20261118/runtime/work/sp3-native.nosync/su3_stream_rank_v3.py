"""Checkpointed random trace sampling for SU(3) classical rank lower bounds.

Native cyclic projection and direct matrix differentiation, without exhaustive formal word enumeration. Every saved independent row has its
seed trace word, so it can be reconstructed with the rational projector.
The terminal condition is an exact-character cohomology upper bound, never
an assertion that a stalled sample exhausted a cochain space.
"""
import sp_environment
import su_trace_adapter_v1 as algebra
import su3_charge_dimensions_v1 as character
import native_su_projection_v1 as native_projection
import native_su_derivation_v1 as native_derivation
from g2_full_classical_screen_v1 import contents,sector
from sp_full_projected_basis import evaluate
from fortuity.full.evalg import Frame,nferm
from fortuity.full.letters import Q_letter
from flint import nmod_mat
from collections import defaultdict
from fractions import Fraction
from math import ceil,comb
from pathlib import Path
import argparse
import hashlib
import json
import os
import pickle
import random
import signal
import time
import numpy as np


def modular(expr,prime):
    return {word:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime
            for word,c in expr.items() if c}


class RowSpace:
    def __init__(self,columns,prime):
        self.columns,self.prime=columns,prime
        self.rows=nmod_mat(0,columns,prime)
        self.pivots=[]

    def extend(self,values):
        if not len(values):
            return []
        matrix=nmod_mat(values.tolist(),self.prime)
        if self.pivots:
            coefficients=nmod_mat([[int(matrix[i,j]) for j in self.pivots] for i in range(matrix.nrows())],self.prime)
            matrix=matrix-coefficients*self.rows
        transposed,rank=matrix.transpose().rref()
        if not rank:
            return []
        selected=[next(j for j in range(transposed.ncols()) if transposed[i,j]) for i in range(rank)]
        reduced,_=matrix.rref()
        newrows=nmod_mat([[int(reduced[i,j]) for j in range(self.columns)] for i in range(rank)],self.prime)
        newpivots=[next(j for j in range(self.columns) if newrows[i,j]) for i in range(rank)]
        if self.pivots:
            coefficients=nmod_mat([[int(self.rows[i,j]) for j in newpivots] for i in range(self.rows.nrows())],self.prime)
            self.rows=self.rows-coefficients*newrows
        pairs=list(zip(self.pivots,self.rows.tolist()))+list(zip(newpivots,newrows.tolist()))
        pairs.sort()
        self.pivots=[pivot for pivot,_ in pairs]
        self.rows=nmod_mat([row for _,row in pairs],self.prime)
        return selected

    def dump(self):
        return dict(columns=self.columns,prime=self.prime,pivots=self.pivots,
                    rows=np.asarray(self.rows.tolist(),dtype=np.uint32).reshape(len(self.pivots),self.columns))

    @classmethod
    def restore(cls,state):
        out=cls(state['columns'],state['prime'])
        out.pivots=state['pivots']
        if out.pivots:
            out.rows=nmod_mat(state['rows'].tolist(),out.prime)
        return out


def random_word(content,rng,max_trace_length=6):
    # Short traces avoid expensive redundant long words. A length cap only
    # changes sampling; no completeness follows from the cap or a plateau.
    letters=list(content)
    rng.shuffle(letters)
    traces=[]
    while letters:
        possible=[length for length in range(2,min(len(letters),max_trace_length)+1)
                  if len(letters)-length!=1]
        size=rng.choice(possible)
        traces.append(tuple(letters[:size]))
        del letters[:size]
    sign,word=algebra.canonical_multi(tuple(traces))
    return word if sign else None


def run(mu,length,parity,prime,output,seconds=300,batch=24,stop_rank=None,resume_from=None,max_trace_length=6):
    output.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    dims=character.dimensions(mu,length+1,parity)
    source_dims=character.dimensions(mu,length,parity)
    all_contents=contents(mu,length)
    by_degree=defaultdict(list)
    for content in all_contents:
        degree=sum(algebra.is_ferm(letter) for letter in content)
        if degree in source_dims:
            by_degree[degree].append(content)
    target_alphabet=sorted({letter for content in contents(mu,length+1) for letter in content})
    required=set(target_alphabet)
    required.update(letter for content in all_contents for letter in content)
    while True:
        enlarged=required.union(replacement for letter in required
                               for _,left,right in Q_letter(letter) for replacement in (left,right))
        if enlarged==required:break
        required=enlarged
    alphabet=sorted(required)
    frames={}
    for degree,dimension in dims.items():
        k=degree+2 if degree else 0
        frames[degree]=Frame(3,max(1,ceil((dimension+24)/comb(k,degree))),k,[degree],alphabet,
                             p=prime,seed=671953+7919*degree+length*31)
    columns=sum(frame.ncols for frame in frames.values())
    sources={Path(path).name:hashlib.sha256(Path(path).read_bytes()).hexdigest()
             for path in (__file__,algebra.__file__,character.__file__,algebra.projection.__file__,native_projection.__file__,native_derivation.__file__)}
    assert max_trace_length>=3
    manifest=dict(mu=mu,length=length,parity=parity,prime=prime,source_dims=source_dims,max_trace_length=max_trace_length,
                  target_dims=dims,columns=columns,sources=sources)
    # JSON conversion makes tuples and integer keys stable across resumes.
    manifest=json.loads(json.dumps(manifest))
    marker=output/'manifest.json'
    if marker.exists():
        assert json.loads(marker.read_text())==manifest
    else:
        marker.write_text(json.dumps(manifest,indent=2)+'\n')
    rng=random.Random(883781+length*97+prime)
    space=RowSpace(columns,prime)
    seeds=[]
    seen=set()
    attempts=projections=0
    elapsed_before=0
    state_path=output/'checkpoint.pkl'
    loaded_path=state_path
    if not loaded_path.exists() and resume_from is not None:
        prior=json.loads((resume_from/'manifest.json').read_text())
        assert all(prior[key]==manifest[key] for key in ('mu','length','parity','prime','source_dims','target_dims','columns'))
        assert alphabet==target_alphabet,'Old frames need an explicit alphabet bridge before resuming'
        loaded_path=resume_from/'checkpoint.pkl'
    if loaded_path.exists():
        state=pickle.loads(loaded_path.read_bytes())
        space=RowSpace.restore(state['space'])
        seeds,seen=state['seeds'],state['seen']
        attempts,projections=state['attempts'],state['projections']
        rng.setstate(state['rng'])
        elapsed_before=state['seconds']
    stop_rank=sum(source_dims.values()) if stop_rank is None else stop_rank
    pending,pending_seeds=[],[]
    degrees=list(by_degree)
    weights=[source_dims[degree] for degree in degrees]
    last_save=start
    timings=dict(projection=0.0,differentiation=0.0,row_reduction=0.0)
    peak_terms=0

    def checkpoint(final=False):
        state=dict(space=space.dump(),seeds=seeds,seen=seen,attempts=attempts,
                   projections=projections,rng=rng.getstate(),seconds=elapsed_before+time.monotonic()-start)
        temporary=state_path.with_suffix('.tmp')
        temporary.write_bytes(pickle.dumps(state,protocol=4))
        temporary.replace(state_path)
        (output/'independent_seeds.pkl').write_bytes(pickle.dumps(seeds,protocol=4))
        result=dict(status='sampled_classical_rank_lower_bound',mu=mu,length=length,parity=parity,
                    prime=prime,source_dimension=sum(source_dims.values()),target_dimension=sum(dims.values()),
                    rank_lower_bound=len(space.pivots),attempts=attempts,nonzero_projections=projections,
                    retained_seeds=len(seeds),peak_projected_source_terms=peak_terms,seconds=state['seconds'],
                    current_run_seconds=time.monotonic()-start,timings=timings,max_trace_length=max_trace_length,
                    stopped_by='rank_target' if len(space.pivots)>=stop_rank else ('time_cap' if final else 'running'),
                    independent_seeds_sha256=hashlib.sha256((output/'independent_seeds.pkl').read_bytes()).hexdigest(),
                    checkpoint_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),
                    source_sha256=sources,
                    scope='Rank lower bound over Q via good reduction. Sampling completeness is not asserted.')
        temporary=output/'result.json.tmp'
        temporary.write_text(json.dumps(result,indent=2)+'\n')
        temporary.replace(output/'result.json')
        print(json.dumps({key:value for key,value in result.items() if key not in ('source_sha256',)}),flush=True)

    def flush():
        if not pending:
            return
        before=time.monotonic()
        parts=[]
        for degree,frame in frames.items():
            parts.append(native_derivation.evaluate(pending,frame))
        values=np.concatenate(parts,axis=1) if parts else np.zeros((len(pending),0),dtype=np.int64)
        timings['differentiation']+=time.monotonic()-before
        before=time.monotonic()
        selected=space.extend(values)
        timings['row_reduction']+=time.monotonic()-before
        seeds.extend(pending_seeds[index] for index in selected)
        assert len(seeds)==len(space.pivots)<=min(sum(source_dims.values()),sum(dims.values()))
        pending.clear()
        pending_seeds.clear()

    while time.monotonic()-start<seconds and len(space.pivots)<stop_rank:
        attempts+=1
        degree=rng.choices(degrees,weights=weights,k=1)[0]
        word=random_word(rng.choice(by_degree[degree]),rng,max_trace_length)
        if word is None or word in seen:
            continue
        seen.add(word)
        expr=algebra.parity_project({word:1},parity)
        if not expr:
            continue
        before=time.monotonic()
        expr=native_projection.project(expr,sector(mu,length),prime)
        timings['projection']+=time.monotonic()-before
        if not expr:
            continue
        assert modular(algebra.conjugate(expr),prime)=={word:(parity*c)%prime for word,c in expr.items()}
        projections+=1
        peak_terms=max(peak_terms,len(expr))
        pending.append(expr)
        pending_seeds.append(word)
        if len(pending)>=batch:
            flush()
        if time.monotonic()-last_save>30:
            flush()
            checkpoint()
            last_save=time.monotonic()
    flush()
    checkpoint(True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--mu',required=True)
    parser.add_argument('--length',required=True,type=int)
    parser.add_argument('--parity',type=int,choices=(-1,1),default=-1)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--seconds',type=int,default=300)
    parser.add_argument('--batch',type=int,default=24)
    parser.add_argument('--stop-rank',type=int)
    parser.add_argument('--resume-from',type=Path)
    parser.add_argument('--max-trace-length',type=int,default=6)
    args=parser.parse_args()
    run(tuple(map(int,args.mu.split(','))),args.length,args.parity,args.prime,args.output,args.seconds,args.batch,args.stop_rank,args.resume_from,args.max_trace_length)
