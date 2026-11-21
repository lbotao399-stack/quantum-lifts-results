"""Solve selected letter contents of the level-32 first equation at three colours."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from math import lcm
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from su3_level32_trace_parents import DATA, formal, incident
import inhomogeneous_block as physical
import joint_content_correction as joint
from su_trace_adapter_v1 import parity_project
from verify import add_scaled


def layer_component(target, scalar_count):
    """All target contents connected through parents with one fewer copy of the first scalar."""
    from g2_full_classical_screen_v1 import contents
    universe = {bytes(sorted(formal.ID[k] for k in raw))
                for raw in contents((1,1,7,3,3),11) if raw.count(formal.ALPH[0]) == scalar_count}
    if target not in universe:
        raise ValueError("selected content is outside the required weight or scalar layer")
    nearby = defaultdict(set)
    for t in universe:
        for parent in incident(t):
            if parent.count(0) == scalar_count-1:
                nearby[parent].add(t)
    component = {target}
    while True:
        neighbours = set().union(*(ts for ts in nearby.values() if component.intersection(ts)))
        if neighbours <= component:
            break
        component.update(neighbours)
    groups = {k.hex():sorted(t.hex() for t in ts) for k,ts in nearby.items() if component.intersection(ts)}
    assert all(set(ts) <= {t.hex() for t in component} for ts in groups.values())
    return sorted(component),groups


def prepare(data, work, upper_band=False, component=None):
    from g2_full_classical_screen_v1 import contents
    targets=[bytes.fromhex("0000000000000505050710")]
    partner=bytes(sorted(formal.ID[b,a,fl] for k in targets[0] for a,b,fl in (formal.ALPH[k],)))
    targets.append(partner)
    assert len(set(targets))==2
    if component is not None:
        seed = bytes.fromhex(component)
        targets, groups = layer_component(seed, seed.count(0))
    direct=defaultdict(dict)
    digest=hashlib.sha256()
    with gzip.open(data/"Q1.jsonl.gz","rb") as stream:
        for line in stream:
            digest.update(line)
            c, raw=json.loads(line)
            w=tuple(bytes(formal.ID[a,b,tuple(f)] for a,b,f in tr) for tr in raw)
            t=bytes(sorted(b"".join(w)))
            if t in targets or (upper_band and t.count(0)>=6):
                direct[t][w]=-Fraction(c)
    assert digest.hexdigest()==json.loads((data/"order1.json").read_bytes())["output_sha256"]
    if upper_band:
        top={bytes(sorted(formal.ID[L] for L in raw)) for raw in contents((1,1,7,3,3),11)
             if raw.count(formal.ALPH[0])==7}
        assert len(top)==13 and not top.intersection(direct)
        targets=sorted(set(direct)|top)
        assert len(targets)==19
    rhs={t:direct[t] for t in targets}
    den=lcm(*(c.denominator for r in rhs.values() for c in r.values()))
    if component is None:
        groups=defaultdict(list)
        for t in targets:
            for g in incident(t):
                groups[g.hex()].append(t.hex())
    dims={key:joint.E.parity_frame_dims(bytes.fromhex(key))[0] for key in list(groups)+[t.hex() for t in targets]}
    work.mkdir(parents=True,exist_ok=True)
    for t,r in rhs.items():
        joint.s.write(work/(t.hex()+".rhs.json.gz"),[[joint.s.hexes(w),str(int(c*den))] for w,c in sorted(r.items())])
    graph=dict(targets=[t.hex() for t in targets],groups=groups,internal=sorted(groups),dimensions=dims,
               denominator=den,rhs=digest.hexdigest(),base="zero")
    if component is not None:
        graph.update(parent_scalar_count=seed.count(0)-1, target_scalar_count=seed.count(0))
    joint.s.write(work/"graph.json",graph)
    return graph


def solve(work, engine, threads):
    graph=joint.s.read(work/"graph.json")
    frames=work/"frames"
    def build(item):
        g=bytes.fromhex(item)
        frame=joint.E.make(g,engine,frames)
        joint.E.check(frame,engine,frames)
        joint.s.write(frames/(item+".json.gz"),frame)
    with ThreadPoolExecutor(threads) as pool:
        list(pool.map(build,sorted(graph["dimensions"])))
    for t in graph["targets"]:
        joint.matrix_for_target(t,frames,engine,work)
    graph,matrix,rhs,pool,keys=joint.assemble(work,"all",frames)
    info,vals=physical.certify(matrix,rhs,lambda row: print(row,flush=True))
    assert info["kind"] in ("exact_solution","exact_separator")
    den=graph["denominator"]
    result=dict(targets=graph["targets"],target_coordinates=len(keys),
                input_sha256=graph["rhs"],rank=info)
    if "parent_scalar_count" in graph:
        result.update(parent_scalar_count=graph["parent_scalar_count"],selected_parent_coordinates=len(pool))
    else:
        result["all_incident_parent_coordinates"] = len(pool)
    if info["kind"]=="exact_solution":
        delta={}
        for i,c in vals.items():
            _,raw=pool[i]
            w=joint.s.expanded(joint.s.unpack_hex(raw))
            add_scaled(delta,parity_project({w:1},1),c/den)
        result["correction"]=[[joint.s.hexes(joint.s.packed(w)),str(c)] for w,c in sorted(delta.items())]
        result["polynomial_replay"]=verify(work,engine,result["correction"])
    else:
        result.update(functional=[[keys[i],str(c)] for i,c in sorted(vals.items())],
                      pairing=str(Fraction(info["integer_rhs_pairing"],den)))
    joint.s.write(work/"content_certificate.json.gz",result)
    return result


def verify(work,engine,correction):
    graph=joint.s.read(work/"graph.json")
    residual={t:{joint.s.unpack_hex(w):-Fraction(int(c),graph["denominator"])
                 for w,c in joint.s.read(work/(t+".rhs.json.gz"))} for t in graph["targets"]}
    for raw,c in correction:
        w=joint.s.unpack_hex(raw)
        for target,k in formal.q(w).items():
            t=bytes(sorted(b"".join(target))).hex()
            if t in residual:
                add_scaled(residual[t],{target:k},Fraction(c))
    checks=physical.polynomials([residual[t] for t in graph["targets"]],engine,work)
    assert all(c["nonzero_monomials"]==0 for c in checks)
    return checks


def verify_certificate(work, engine, frames_archive=None):
    graph = joint.s.read(work/"graph.json")
    saved = joint.s.read(work/"content_certificate.json.gz")
    assert saved["input_sha256"] == graph["rhs"] and saved["targets"] == graph["targets"]
    if saved["rank"]["kind"] == "exact_solution":
        checks = verify(work, engine, saved["correction"])
        assert checks == saved["polynomial_replay"]
        return checks
    if saved["rank"]["kind"] != "exact_separator":
        raise ValueError("unknown certificate type")
    frames = work/"frames"
    frames.mkdir(exist_ok=True)
    if frames_archive:
        entries = joint.s.read(frames_archive)
        assert set(entries) == {item+".json.gz" for item in graph["dimensions"]}
        for filename, entry in entries.items():
            joint.s.write(frames/filename, entry)
    for item in graph["dimensions"]:
        entry = joint.s.read(frames/(item+".json.gz"))
        assert entry["content"] == item
        joint.E.check(entry, engine, frames)
    for target in graph["targets"]:
        joint.matrix_for_target(target, frames, engine, work)
    _, matrix, rhs, _, keys = joint.assemble(work, "all", frames)
    coordinates = {tuple(k):Fraction(c) for k,c in saved["functional"]}
    assert len(coordinates) == len(saved["functional"]) and set(coordinates) <= set(keys)
    den = lcm(*(c.denominator for c in coordinates.values()))
    vector = [int(den*coordinates.get(key,0)) for key in keys]
    assert not any(joint.multiply(matrix,vector))
    pairing = Fraction(sum(x*y for x,y in zip(vector,rhs,strict=True)), den*graph["denominator"])
    assert pairing and pairing == Fraction(saved["pairing"])
    return dict(annihilated_parent_coordinates=matrix.shape[0], pairing=str(pairing))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("prepare","solve","verify"))
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--engine",type=Path)
    p.add_argument("--threads",type=int,default=2)
    p.add_argument("--frames-archive",type=Path,help="archived colour frames for a separating functional")
    selection=p.add_mutually_exclusive_group()
    selection.add_argument("--upper-band",action="store_true",
                   help="include all six inhomogeneous six-scalar and all thirteen seven-scalar contents")
    selection.add_argument("--layer-component",help="seed hex content for the parent-restricted scalar layer")
    a=p.parse_args()
    if not __debug__ or (a.mode!="prepare" and not a.engine):
        p.error("enable assertions and supply a colour_columns executable for the matrix operations")
    if a.mode=="prepare":
        print(prepare(a.data,a.work,a.upper_band,a.layer_component)["dimensions"])
    elif a.mode=="solve":
        print(solve(a.work,a.engine.resolve(),a.threads)["rank"])
    else:
        print(verify_certificate(a.work,a.engine.resolve(),a.frames_archive))
