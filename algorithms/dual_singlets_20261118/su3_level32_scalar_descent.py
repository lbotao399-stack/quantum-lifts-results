"""Cancel the five-scalar residual by four-scalar trace parents within each trace component."""
import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
import hashlib
from pathlib import Path

from su3_level32_layers import content, DATA, first_targets, formal, joint, residual
from su3_level32_levi import extensions
from su3_level32_centralizer import general_colour_test
from su3_level32_content import physical
from verify import add_scaled


def serialized(row):
    return [[joint.s.hexes(w),str(c)] for w,c in sorted(row.items())]


def restored(rows):
    return {joint.s.unpack_hex(w):Fraction(c) for w,c in rows}


def prepare(data, correction, folder):
    terms = restored(joint.s.read(correction))
    assert all(content(w).count(0)==5 for w in terms)
    rows = residual(data,terms,lambda t:t.count(0)==5)
    erase = {m:k for k,(m,_) in extensions().items()}
    groups = defaultdict(dict)
    for t,row in rows.items():
        bare = bytes(sorted(erase.get(i,i) for i in t if i))
        groups[bare.hex()][t.hex()] = serialized(row)
    folder.mkdir(exist_ok=True,parents=True)
    entries = {}
    for bare,values in sorted(groups.items()):
        path = folder/(bare+".json.gz")
        joint.s.write(path,values)
        entries[bare] = dict(terms=sum(map(len,values.values())),contents=len(values),
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    assert len(groups)==14 and sum(x["contents"] for x in entries.values())==119
    joint.s.write(folder/"groups.json",dict(correction_sha256=hashlib.sha256(correction.read_bytes()).hexdigest(),
                                            groups=entries))


@lru_cache(maxsize=10000)
def scalar_splits(w):
    return {v:c for v,c in formal.q(w).items() if sum(t.count(0) for t in v)==5}


def solve_component(seed):
    targets={seed}
    queue=deque([seed])
    sources={}
    while queue:
        w=queue.popleft()
        for u in formal.parents(w):
            if u in sources or sum(t.count(0) for t in u)!=4:
                continue
            image=scalar_splits(u)
            if w not in image:
                continue
            sources[u]=image
            for t in image:
                if t not in targets:
                    targets.add(t)
                    queue.append(t)
    return targets,sources


def reduce_component(rhs, sources):
    pivots={}
    for u,image in sorted(sources.items()):
        row=image.copy()
        operations={}
        while row:
            k=min(row)
            if k not in pivots:
                pivots[k]=row,u,operations
                break
            earlier=pivots[k][0]
            scale=Fraction(row[k],earlier[k])
            add_scaled(row,earlier,-scale)
            operations[k]=scale
    remaining={}
    active=rhs.copy()
    solution={}
    while active:
        k=min(active)
        if k not in pivots:
            remaining[k]=active.pop(k)
            continue
        row=pivots[k][0]
        scale=Fraction(active[k],row[k])
        add_scaled(active,row,-scale)
        solution[k]=scale
    delta={}
    for k in sorted(pivots,reverse=True):
        if solution.get(k):
            _,u,operations=pivots[k]
            delta[u]=-solution[k]
            add_scaled(solution,operations,-solution[k])
    return delta,remaining


def descend(folder, bare, engine, work, full_solve):
    entry=joint.s.read(folder/"groups.json")["groups"][bare]
    path=folder/(bare+".json.gz")
    assert hashlib.sha256(path.read_bytes()).hexdigest()==entry["sha256"]
    values=joint.s.read(path)
    original={w:c for v in values.values() for w,c in restored(v).items()}
    assert len(original)==entry["terms"]
    delta={}
    remaining={}
    seen_components=0
    if full_solve:
        pending=set(original)
        while pending:
            seed=min(pending)
            targets,sources=solve_component(seed)
            hits=targets & pending
            pending.difference_update(hits)
            part,rest=reduce_component({w:original[w] for w in hits},sources)
            assert delta.keys().isdisjoint(part) and remaining.keys().isdisjoint(rest)
            delta.update(part)
            remaining.update(rest)
            seen_components+=1
    else:
        remaining=original.copy()
    checked=original.copy()
    for w,c in delta.items():
        add_scaled(checked,scalar_splits(w),c)
    assert checked==remaining
    grouped=defaultdict(dict)
    for w,c in remaining.items():
        grouped[content(w).hex()][w]=c
    assert set(grouped)<=set(values)
    work.mkdir(exist_ok=True,parents=True)
    checks={}
    for t in sorted(values):
        row=grouped[t]
        if not row:
            checks[t]={"kind":"formal_zero"}
            continue
        for seed in (894611,620003):
            value,_=general_colour_test(bytes.fromhex(t),row,seed)
            if Fraction(value):
                checks[t]=dict(kind="nonzero",seed=seed,coefficient=value)
                break
        else:
            exact=physical.polynomials([row],engine.resolve(),work)[0]
            checks[t]=dict(kind="nonzero" if exact["nonzero_monomials"] else "polynomial_zero",exact=exact)
    result=dict(bare=bare,input_sha256=entry["sha256"],reduced=full_solve,components=seen_components,
                original_words=len(original),residual_words=len(remaining),correction=serialized(delta),checks=checks)
    joint.s.write(work/(bare+"_result.json.gz"),result)
    print(bare,len(original),len(delta),len(remaining),dict(Counter(r["kind"] for r in checks.values())),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("prepare","solve","test"))
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--folder",type=Path,required=True)
    p.add_argument("--correction",type=Path)
    p.add_argument("--bare")
    p.add_argument("--engine",type=Path)
    p.add_argument("--work",type=Path)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    if args.mode=="prepare" and args.correction:
        prepare(args.data,args.correction,args.folder)
    elif args.mode!="prepare" and args.bare and args.engine and args.work:
        descend(args.folder,args.bare,args.engine,args.work,args.mode=="solve")
    else:
        p.error("prepare requires a correction and solve/test require a bare content, engine and work directory")
