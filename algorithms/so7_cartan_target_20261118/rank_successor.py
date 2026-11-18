"""Certify a minor of Q0 from length eleven to twelve at weight (3,3,2,2,2)."""
import argparse
import gzip
import hashlib
import json
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor
from math import ceil, comb
from pathlib import Path

import numpy as np
from flint import nmod_mat
from reduce_incoming import Echelon, P

HERE = Path(__file__).resolve().parent
MU = (3,3,2,2,2)


def setup(runtime, native, plan):
    global choices, frame, project, canonical, classical, evaluate, charge
    sys.path[:0] = [str(native), str(runtime)]
    from sp_environment import canonical, fw
    from native_full_projection_sparse_alphabet_v1 import project
    from native_so7_terminal_trace_v1 import evaluate
    from g2_full_classical_screen_v1 import contents, sector
    from sp_full_profiles import profile
    from fortuity.full.evalg import Frame, nferm
    classical = fw.Q_expr
    charge = sector(MU, 11)
    profiles = json.loads(Path(plan).read_text())["targets"]
    choices = []
    for row in profiles:
        pattern = tuple((tuple(k),v) for k,v in row["profile"])
        choices.append([w for w in contents(MU,11) if profile(w) == pattern])
    target = contents(MU,12)
    assert all(nferm((w,)) == 6 for w in target)
    letters = sorted({l for w in target for l in w})
    frame = Frame(3,ceil((2308+96)/comb(9,6)),9,[6],letters,seed=81118811,p=P,group="SO_odd")


def action(word):
    expression = project({word:1},charge,P)
    return {w:c%P for w,c in canonical(classical(expression)).items() if c%P}


def sample(serial):
    rng = random.Random(81118819 + 84991*serial)
    words, images = [], []
    for _ in range(8):
        while True:
            letters = list(rng.choice(choices[rng.randrange(3)]))
            rng.shuffle(letters)
            traces = []
            while letters:
                candidates = [n for n in range(2,len(letters)+1) if len(letters)-n != 1]
                n = rng.choice(candidates)
                traces.append(tuple(letters[:n]))
                letters = letters[n:]
            word = next(iter(canonical({tuple(traces):1})),None)
            if word:
                image = action(word)
                if image:
                    words.append(word)
                    images.append(image)
                    break
    return words, evaluate(images,frame).astype(np.uint32)


def decode(word):
    return tuple(tuple((i,j,tuple(s)) for i,j,s in t) for t in word)


def controls():
    import full_extremal
    from sp_full_projected_basis import evaluate as reference_evaluate
    words, actual = sample(0)
    for w, value in zip(words[:2], actual):
        exact = full_extremal.project({w:1}, charge)
        assert not canonical(classical(canonical(classical(exact))))
        def residue(c):
            from fractions import Fraction
            c = Fraction(c)
            return c.numerator * pow(c.denominator,-1,P) % P
        expected = {k:residue(c) for k,c in canonical(classical(exact)).items() if residue(c)}
        assert expected == action(w)
        # Independent Python evaluator at the first Cartan-unrestricted matrix point.
        from copy import copy
        point = copy(frame)
        point.P, point.ncols = 1, frame.cpp
        point.L = {l:(d,arr[:1]) for l,(d,arr) in frame.L.items()}
        row = reference_evaluate([expected],point)[0]
        assert np.array_equal(row, value[:frame.cpp])


def compute(args):
    init = tuple(str(p.resolve()) for p in (args.runtime,args.native,args.plan))
    setup(*init)
    controls()
    if args.control:
        print("Independent classical actions and point evaluations passed.",flush=True)
        return
    args.output.mkdir(parents=True,exist_ok=True)
    if args.replay:
        saved = json.loads(gzip.decompress(args.replay.read_bytes()))
        columns = saved["columns"]
        words = [decode(w) for w in saved["words"]]
        rows = []
        for i in range(0,len(words),8):
            rows.extend(evaluate([action(w) for w in words[i:i+8]],frame)[:,columns].tolist())
        det = int(nmod_mat(rows,P).det())
        assert len(words)==len(columns)==len(set(columns)) and det==saved["determinant"] and det
        print(json.dumps(dict(order=len(words),determinant=det,prime=P)),flush=True)
        return
    space = Echelon(frame.ncols,True)
    words, values = [], []
    with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=init) as pool:
        for start in range(0,args.batches,2*args.workers):
            stop = min(start+2*args.workers,args.batches)
            for recipes,rows in pool.map(sample,range(start,stop)):
                words.extend(recipes)
                values.extend(rows)
            n=8*start
            block=np.column_stack((np.asarray(values[n:]),np.zeros(len(words)-n,dtype=np.uint32)))
            space.add(block,list(range(n,len(words))))
            print(json.dumps(dict(rows=len(words),rank=space.rank)),flush=True)
            if space.rank>=2307:
                break
    ids = [i for panel in space.panels for i in panel["ids"]]
    columns = [int(i) for panel in space.panels for i in panel["pivots"]]
    assert len(set(ids))==len(set(columns))==space.rank
    minor = np.asarray([values[i][columns] for i in ids],dtype=np.uint32)
    det = int(nmod_mat(minor.tolist(),P).det())
    assert det
    record = dict(prime=P,weight=MU,length=11,dimension=14832,next_dimension=2308,
                  rows_generated=len(words),minor_order=space.rank,determinant=det,
                  previous_rank_upper=14832-space.rank,
                  minor_sha256=hashlib.sha256(minor.astype("<u4").tobytes()).hexdigest())
    archive = gzip.compress(json.dumps(dict(columns=columns,words=[words[i] for i in ids],
                                     determinant=det),separators=(",",":")).encode(),mtime=0)
    (args.output/"minor.json.gz").write_bytes(archive)
    record["recipes_sha256"]=hashlib.sha256(archive).hexdigest()
    (args.output/"result.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(record),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runtime",type=Path,required=True)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--plan",type=Path,default=HERE/"source_profiles.json")
    p.add_argument("--output",type=Path,default=Path("output"))
    p.add_argument("--workers",type=int,default=1)
    p.add_argument("--batches",type=int,default=2000)
    p.add_argument("--control",action="store_true")
    p.add_argument("--replay",type=Path)
    compute(p.parse_args())
