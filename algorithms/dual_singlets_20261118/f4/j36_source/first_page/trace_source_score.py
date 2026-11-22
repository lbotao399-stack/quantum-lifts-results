"""Evaluate the corrected F4 row on explicit spin/flavour singlets and their classical cycles."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
import gzip
from itertools import permutations, product
import json
from math import prod
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np

from flag_boundaries import HERE, adjoint
from ghost_trace_quotient import rotate, sign
from single_trace_flag import native26
from trace_source_cycles import GLUE, PAIRS, odd


def section(f):
    ghosts = tuple(k for k in f if k in (3,4))
    rest = f[len(ghosts):]
    # Match the cyclic scalar word by averaging its rotations relative to the
    # ghost block. A pair of ghosts is projected to the dotted determinant.
    circles = [rest[i:]+rest[:i] for i in range(len(rest))] if ghosts and rest else [rest]
    out = Counter()
    for g in permutations(ghosts):
        for w in circles:
            out[g+w] += sign(g)*Fraction(1,len(circles)*(2 if len(g)==2 else 1))
    return tuple(out.items())


@cache
def numeric(w):
    a,b = rotate(w),rotate(w[::-1])
    if len(w)%2 and a==b:
        return 0, ()
    return (-1 if b<a and len(w)%2 else 1),min(a,b)


def templates(source):
    k = source.count(GLUE)
    count = 4-k
    out = Counter()
    for terms in product(*(section(f) for f in source if f!=GLUE)):
        factors = [v for v,_ in terms]
        coeff = prod(c for _,c in terms)
        word = sum(factors,())
        site = next(i for i,v in enumerate(word) if v>=5)
        missing, = set(range(3))-set(PAIRS[word[site]-5])
        orientation = (-1 if missing!=1 else 1)*(-1)**sum(map(odd,word[:site]))
        for active in (i for i,v in enumerate(word) if v==missing):
            for plus,minus in product(permutations(range(5,5+count)),permutations(range(5+count,5+2*count))):
                ghosts = {3:iter(plus),4:iter(minus)}
                labels = tuple(3 if i==site else 4 if i==active else
                               next(ghosts[v]) if v in ghosts else v for i,v in enumerate(word))
                assert sorted(labels)==list(range(len(word)))
                phase = orientation * sign(tuple(v for v in labels if v>=5))
                traces = []
                first = 0
                for f in factors:
                    end = first+len(f)
                    s,canonical = numeric(labels[first:end])
                    phase *= s
                    traces.append(canonical)
                    first = end
                out[tuple(sorted(traces))] += phase*coeff
    return k,{p:v for p,v in out.items() if v}


def source_patterns(data):
    all_patterns = defaultdict(dict)
    for i,raw in enumerate(data["monomials"]):
        k,polys = templates(tuple(map(tuple,raw)))
        for p,v in polys.items():
            all_patterns[k].setdefault(p,{})[i] = v
    return all_patterns


def native_data():
    data = json.loads(gzip.decompress((HERE/"inputs.json.gz").read_bytes()))
    twice = np.zeros((52,52),dtype=np.int64)
    for i,j,n in data["inverse"]:
        v = 2*Fraction(n)
        assert v.denominator==1
        twice[i,j]=v.numerator
    dual = (twice@native26().reshape(52,-1)).reshape(52,26,26).transpose(0,2,1)
    assert np.array_equal(np.einsum("aij,bji->ab",dual,dual),12*twice)
    adj = adjoint()
    omega = defaultdict(Counter)
    for i,j in zip(*np.nonzero(twice)):
        x,y = adj[i],adj[j]
        for u,a in zip(*np.nonzero(x)):
            for v,b in zip(*np.nonzero(y)):
                value = int(twice[i,j]*x[u,a]*y[v,b])
                omega[int(u),int(v)][int(a),int(b)] += value
                omega[int(v),int(u)][int(a),int(b)] += value
    return twice,dual,{key:tuple((a,b,v) for (a,b),v in row.items() if v) for key,row in omega.items()}


def parent_rows(frame,twice,omega):
    base = {tuple(key):int(x) for key,x in frame["entries"]}
    for k in (2,3,4):
        n = 5-k
        parents = Counter()
        next_frame = Counter()
        for key,x in base.items():
            plus,minus,scalars = key[:n],key[n:2*n],key[2*n:]
            for i,u in enumerate(plus):
                for j,v in enumerate(minus):
                    rest = plus[:i]+plus[i+1:]+minus[:j]+minus[j+1:]
                    signed = (-1)**(i+n+j-1)*x
                    if k<4 and twice[u,v]:
                        next_frame[rest+scalars] += 3*int(twice[u,v])*signed
                    for a,b,c in omega.get((u,v),()):
                        parents[bytes((*scalars,a,b,*rest))] += signed*c
        yield k,{key:v for key,v in parents.items() if v}
        base = {key:v for key,v in next_frame.items() if v}


def write_problem(folder,dual,k,rows,patterns):
    # The C++ pass contracts integer native traces on the distinct colour
    # assignments. Python retains the unbounded rational dual coefficients.
    values = list(dict.fromkeys(rows.values()))
    assert len(values)<65536
    lookup = {x:i.to_bytes(2,"little") for i,x in enumerate(values)}
    with (folder/"query").open("wb") as out:
        out.write(f"{len(rows)} {5+2*(4-k)} {len(values)}\n".encode())
        for key,value in rows.items():
            out.write(key+lookup[value])
    with (folder/"patterns").open("w") as out:
        print(len(patterns),file=out)
        for p in patterns:
            print(len(p),*(v for t in p for v in (len(t),*t)),file=out)
    with (folder/"rep").open("w") as out:
        for m in dual:
            i,j = np.nonzero(m)
            print(len(i),file=out)
            for a,b in zip(i,j):
                print(a,b,m[a,b],file=out)
        weights = adjoint()[:4].diagonal(axis1=1,axis2=2).T
        assert 9*np.abs(weights).max()<127
        for w in weights:
            print(sum(int(x)*127**i for i,x in enumerate(w)),file=out)
    return values


def run():
    data = json.loads(gzip.decompress((HERE/"trace_source_cycles.json.gz").read_bytes()))
    frame = json.loads(gzip.decompress((HERE/"trace_score_frame.json.gz").read_bytes()))
    twice,dual,omega = native_data()
    patterns = source_patterns(data)
    results = Counter()
    audit = []
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        exe = folder/"evaluate"
        subprocess.run(["c++","-std=c++17","-O3","-ftrapv",str(HERE/"trace_source_contract.cpp"),"-o",str(exe)],check=True)
        for k,rows in parent_rows(frame,twice,omega):
            monomials = list(patterns[k])
            weights = write_problem(folder,dual,k,rows,monomials)
            subprocess.run([str(exe)],cwd=folder,check=True)
            totals = [0]*len(monomials)
            pairs = 0
            with (folder/"values").open() as inp:
                for value,record in zip(weights,inp,strict=True):
                    entries = list(map(int,record.split()))
                    assert len(entries)%2==0
                    pairs += len(entries)//2
                    for j,z in zip(entries[::2],entries[1::2]):
                        totals[j] += value*z
            den = frame["denominator"]*4*2**(5+2*(4-k))
            for p,v in zip(monomials,totals):
                for i,c in patterns[k][p].items():
                    results[i] += c*Fraction(v,den)
            audit.append(dict(factored_traces=k,colour_assignments=len(rows),distinct_row_weights=len(weights),
                              patterns=len(monomials),nonzero_weight_pattern_buckets=pairs))
    def scores(cols):
        return [str(sum(results[i]*x for i,x in c)) for c in cols]
    return dict(native_metric="Tr26/6",intermediate_covector="corrected folded eta",
                singlet_scores=scores(data["singlet_columns"]),classical_cycle_scores=scores(data["cycle_columns"]),
                monomial_scores=[str(results[i]) for i in range(len(data["monomials"]))],work=audit)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(run(),indent=2))
