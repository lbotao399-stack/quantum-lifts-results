"""Determine the first quantum product coefficient modulo formal classical boundaries."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import mixed_incoming as m
from dual633_pairing import integral
from gaugino_primary_order import charge
from mixed_formal import degree
from orthosymplectic_classical_page import minor
from primitive_replay import load_primitive
from so8_level24_chain import solve_image
from sp3_omega_replay import pair
from symplectic_primitive import HERE
from trace_descent import action, highest_weights, settings


LOWER = (1,1,5,2,2)


def grouped_basis(length):
    basis = highest_weights(LOWER, length)
    grouped = defaultdict(list)
    for fermions in sorted(basis):
        for row in basis[fermions]:
            label, = {degree(w) for w in row}
            grouped[label].append(row)
    return grouped


def coordinate_words(rows):
    """Select trace coefficients separating every formal highest-weight row over Q."""
    reduced = {}
    for row in rows:
        v = {w:c % m.P for w,c in row.items() if c % m.P}
        while v:
            w = min(v)
            c = v[w]
            if w not in reduced:
                inverse = pow(c,-1,m.P)
                reduced[w] = {x:n*inverse % m.P for x,n in v.items()}
                break
            for x,n in reduced[w].items():
                new = (v.get(x,0)-c*n) % m.P
                if new:
                    v[x] = new
                elif x in v:
                    del v[x]
        else:
            raise ValueError("dependent formal highest-weight basis")
    return sorted(reduced)


def classical_matrix(rows, words):
    ids = {w:i for i,w in enumerate(words)}
    matrix = np.zeros((len(rows),len(words)),dtype=np.int64)
    for j,row in enumerate(rows):
        for w,n in m.q(row).items():
            if w in ids:
                matrix[j,ids[w]] = n
    return matrix


def rank_minor(matrix):
    if not matrix.size or not np.any(matrix):
        return dict(rows=[],columns=[],rank=0,determinant=1)
    return minor(matrix % m.P, m.P)


def formal_correction(polynomial, archive):
    first = action(polynomial,1,6,"Sp")
    assert charge(first) == (LOWER,10)
    previous, centre, following = [grouped_basis(n) for n in (8,9,10)]
    rhs = defaultdict(dict)
    for w,c in first.items():
        rhs[degree(w)][w] = -c
    labels = sorted(set(previous) | set(centre) | set(following) | set(rhs))
    correction = Counter()
    blocks = []
    for label in labels:
        src = centre[label]
        if not src and not rhs[label]:
            continue
        centre_words = coordinate_words(src)
        target_words = coordinate_words(following[label])
        incoming = classical_matrix(previous[label],centre_words)
        outgoing = classical_matrix(src,target_words)
        before, after = rank_minor(incoming),rank_minor(outgoing)
        assert before["rank"] + after["rank"] == len(src), (label,len(src),before["rank"],after["rank"])
        coefficients = {}
        if rhs[label]:
            ids = after["rows"]
            values = [rhs[label].get(w,0) for w in target_words]
            sol,_ = solve_image(outgoing[ids],values)
            coefficients = {ids[i]:c for i,c in sol.items()}
            for i,c in coefficients.items():
                for w,k in src[i].items():
                    correction[w] += c*k
        blocks.append(dict(trace_degrees=label,source_eight=len(previous[label]),centre_nine=len(src),
                           following_ten=len(following[label]),incoming=before,outgoing=after,
                           first_coefficients=[[i,str(c)] for i,c in sorted(coefficients.items())]))
        print("formal block",label,len(src),before["rank"],after["rank"],len(coefficients),flush=True)
    correction = {w:c for w,c in correction.items() if c}
    assert m.q(correction) == {w:-c for w,c in first.items()}
    data = dict(weight=LOWER,prime=m.P,first_action_terms=len(first),blocks=blocks,
                correction_terms=len(correction),correction_denominator=integral(correction)[1],
                formal_residual_terms=0)
    m.write(archive,data)
    return correction,data


def run(primitive_path, archive, work, workers):
    settings()
    primitive,polynomial = load_primitive(primitive_path)
    correction,data = formal_correction(polynomial,archive)
    first_primitive = action(primitive,1,6,"Sp")
    closed = Counter(correction)
    closed.subtract(first_primitive)
    closed = {w:c for w,c in closed.items() if c}
    lower = HERE.parents[1] / "sp3_omega/lowered.json.gz"
    page = json.loads(gzip.decompress(lower.read_bytes()))
    lib = m.library(work,6,m.P)
    p1, q1k, diff = [pair(row,page,lib,workers) for row in (correction,first_primitive,closed)]
    assert (p1-q1k) % m.P == diff
    return dict(primitive_sha256=hashlib.sha256(primitive_path.read_bytes()).hexdigest(),
                formal_correction_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                lowered_page_sha256=hashlib.sha256(lower.read_bytes()).hexdigest(),
                formal_blocks=len(data["blocks"]),formal_source=sum(b["centre_nine"] for b in data["blocks"]),
                formal_incoming_rank=sum(b["incoming"]["rank"] for b in data["blocks"]),
                formal_outgoing_rank=sum(b["outgoing"]["rank"] for b in data["blocks"]),
                formal_first_action_terms=data["first_action_terms"],formal_correction_terms=len(correction),
                formal_residual_terms=0,reduced_difference_terms=len(closed),prime=m.P,
                product_pairing=p1,primitive_first_pairing=q1k,difference_pairing=diff)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("primitive","archive","work","out"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    args = p.parse_args()
    if not __debug__: p.error("enable assertions")
    args.work.mkdir(parents=True,exist_ok=True)
    result = run(args.primitive,args.archive,args.work,args.workers)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
