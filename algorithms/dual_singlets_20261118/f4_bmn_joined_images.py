"""F4 BMN boundaries from adjoint and Jordan contractions, including a spectator trace."""
import argparse
import hashlib
import random
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_bmn_primary import LETTERS, SIGNS, WEIGHT, bmn_convention
from f4_bmn_trace_images import baseline
from f4_dressed_jordan import short_arm
from f4_joined_boundaries import graph, native, vary
from f4_scalar_correction import evaluate_rows
from f4_singlet_characters import singlets
from mixed33_exact import read


@lru_cache(maxsize=100000)
def canonical(words, channel):
    if channel == "adjoint":
        return graph(words,channel)
    (s,first),(t,second) = (short_arm(w) for w in words[:2])
    if not s*t:
        return 0,words
    u, result = graph((first,second)+words[2:],channel)
    return s*t*u,result


@lru_cache(None)
def arm_sizes(length,channel,spectator):
    minimum = 4 if channel=="jordan" else 2
    ends = range(2,length-2*minimum+1) if spectator else [0]
    return [(a,length-end-a,end) for end in ends for a in range(minimum,(length-end)//2+1)]


def recipe(length, channel, spectator, seed):
    rng = random.Random(seed)
    kind = "p" if length==11 else rng.choice(("f","pp","pp"))
    permutation = list(range(length))
    rng.shuffle(permutation)
    a,b,end = rng.choice(arm_sizes(length,channel,spectator))
    result = Counter()
    for letters, coeff in scalar.assignments(kind):
        ordered = tuple(LETTERS[letters[i]] for i in permutation)
        for i in letters:
            coeff *= SIGNS[i]
        words = (ordered[:a],ordered[a:a+b])
        if end:
            words += (ordered[a+b:],)
        s,g = canonical(words,channel)
        if s:
            result[g] += coeff*s
    return {g:c for g,c in result.items() if c}


def reduce_graphs(expr,channel):
    result = Counter()
    for graph_word, coeff in expr.items():
        s,g = canonical(graph_word,channel)
        result[g] += s*coeff
    return {g:c for g,c in result.items() if c}


def run(previous,work,workers,draws):
    work.mkdir(exist_ok=True,parents=True)
    old = read(previous)
    p = old["prime"]
    bmn_convention()
    lib = {channel:native(work,p,channel) for channel in ("jordan","adjoint")}
    closed = {27:m.library(work,27,p)}
    records = []
    for rec in old["records"]:
        length = rec["length"]
        dims = singlets("F4", WEIGHT, length+1)
        degree = next(iter(dims))
        alphabet = set(LETTERS[:3+(3 if degree else 0)])
        mat, labels = baseline(rec,dims,p,closed,lib["jordan"],work,workers,715180+length*73)
        limit = 61 if length==11 else 141-records[0]["minor"]["rank"]
        history = []
        for channel,spectator in (("jordan",True),("adjoint",False),("adjoint",True)):
            for offset in range(0,draws,56):
                seeds = range(3700001+length*1123+offset,3700001+length*1123+min(offset+56,draws))
                labels_batch, images = [], []
                for seed in seeds:
                    row = recipe(length,channel,spectator,seed)
                    if not row:
                        continue
                    assert all(not reduce_graphs(vary(row,axis,channel),channel) for axis in (1,2))
                    image = reduce_graphs(vary(row,0,channel),channel)
                    if image:
                        assert not reduce_graphs(vary(image,0,channel),channel)
                        images.append(image)
                        labels_batch.append((channel,spectator,seed))
                if not images:
                    continue
                block,_ = evaluate_rows(images,alphabet,rec["frames"],p,work,workers,
                                        lib[channel],spectator,channel)
                mat = np.vstack((mat,block))
                labels += labels_batch
                minor = m.minor(mat,p)
                piv = minor["rows"]
                mat = mat[piv]
                labels = [labels[i] for i in piv]
                assert len(piv)<=limit
                history.append([channel,spectator,seeds.stop-1,len(piv)])
                print(length,channel,spectator,len(labels_batch),len(piv),flush=True)
                if len(piv)==limit:
                    break
            if len(labels)==limit:
                break
        records.append(dict(length=length,rows=labels,minor=m.minor(mat,p),history=history))
        m.write(work/"joined_images.json",dict(input_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
                prime=p,weight=WEIGHT,dimensions=old["dimensions"],records=records,
                classical_cohomology_upper=141-sum(r["minor"]["rank"] for r in records) if len(records)==2 else None))
    return [r["minor"]["rank"] for r in records]


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previous",type=Path,default=Path(__file__).with_name("f4")/"primary.json")
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=12)
    p.add_argument("--draws",type=int,default=392)
    args=p.parse_args()
    assert __debug__ and min(args.workers,args.draws)>0
    print(run(**vars(args)),flush=True)
