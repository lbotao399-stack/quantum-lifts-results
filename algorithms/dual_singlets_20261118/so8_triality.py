"""Triality on the rational reflection-even SO(8) singlet complex at level 18."""
import argparse
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb, factorial
from pathlib import Path

import numpy as np
from flint import fmpq_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from fortuity.full.letters import is_ferm, letters_up_to
from terminal_interval_ranks import IntegralTerminalFrame

HERE = Path(__file__).resolve().parent / "so8_level18"


def representations():
    """Images of twice E_ij-E_ji acting on the two chiral eight-dimensional spaces."""
    eye = np.eye(8, dtype=np.int64)
    units = np.zeros((8, 8, 8), dtype=np.int64)
    units[0] = eye
    for i in range(1, 8):
        units[i, i, 0] = 1
        units[i, 0, i] = -1
    for i, j, k in ((1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)):
        for a,b,c in ((i,j,k),(j,k,i),(k,i,j)):
            units[a,c,b], units[a,b,c] = 1,-1
    for i in range(8):
        for j in range(8):
            anticommutator = units[i]@units[j].T + units[j]@units[i].T
            assert np.array_equal(anticommutator, 2*int(i==j)*eye)
    pairs = list(combinations(range(8), 2))
    vector = np.array([np.outer(eye[i],eye[j])-np.outer(eye[j],eye[i]) for i,j in pairs])
    spinors = [np.array([u[i]@u[j].T for i,j in pairs]) for u in (units, units.transpose(0,2,1))]
    upper = np.triu_indices(8,1)
    for spin in spinors:
        assert not (spin + spin.transpose(0,2,1)).any()
        for i,j in combinations(range(28),2):
            bracket = vector[i]@vector[j] - vector[j]@vector[i]
            assert np.array_equal(spin[i]@spin[j]-spin[j]@spin[i],
                                  2*np.tensordot(bracket[upper],spin,axes=1))
    return spinors


def sparse_matrix(rows, columns):
    indices, values, starts = [], [], [0]
    for row in rows:
        for j, n in row:
            indices.append(j)
            values.append(n)
        starts.append(len(indices))
    return csr_matrix((np.array(values, dtype=np.int64),indices,starts),shape=(len(rows),columns))


def matrix_values(rows, words, frame, degree, length):
    all_used = {letter for word in words for trace in word for letter in trace}
    norm = max(int(np.abs(frame.L[letter][1]).sum(axis=-1).max()) for letter in all_used)
    # Limits every intermediate trace/exterior product, including a prefix
    # whose degree is smaller than the final one.
    bound = factorial(degree) * max(1,norm)**length * 8**max(map(len,words))
    assert bound < 2**63
    evaluated = frame.eval_words(words, dtype=np.int64)
    total = max(sum(abs(n) for _,n in row) for row in rows)
    if int(np.abs(evaluated).max()) * total < 2**63:
        return sparse_matrix(rows, len(words)) @ evaluated
    return np.asarray([[sum(n*int(evaluated[j,k]) for j,n in row)
                        for k in range(evaluated.shape[1])] for row in rows], dtype=object)


def projector(key, data, spinors):
    n, r = int(key[1]), int(key[4])
    words = [m.decode(w) for w in data["words"]]
    rows, length = data["rows"], data["letters"]
    assert length == 9-n-r and all(sum(map(len, w))==length for w in words)
    alphabet = sorted(l for l in letters_up_to(2*r,2*r) if l[0]+l[1]<=2*r)
    upper = np.triu_indices(8,1)
    p = fmpq_mat(len(rows),len(rows))
    all_values = []
    for rec in data["controls"]:
        d = rec["degree"]
        points = rec["columns"]//comb(d+2,d)
        ids = [i for i,row in enumerate(rows) if m.nferm(words[row[0][0]])==d]
        assert len(ids)==rec["character_dimension"]
        chosen = sorted({j for i in ids for j,_ in rows[i]})
        assert all(m.nferm(words[j])==d for j in chosen)
        remap = {j:k for k,j in enumerate(chosen)}
        local = [[(remap[j],v) for j,v in rows[i]] for i in ids]
        frame = IntegralTerminalFrame(4,points,d+2,[d],[],seed=1,group="SO_odd")
        rng = np.random.default_rng(rec["seed"])
        vector = {}
        for l in alphabet:
            odd = int(is_ferm(l))
            a = np.triu(rng.integers(-1,2,size=(points,d+2 if odd else 1,8,8),dtype=np.int64),1)
            vector[l] = (odd,a-a.swapaxes(-2,-1))
        components = []
        for spin in (None,*spinors):
            frame.L = vector if spin is None else {
                l:(odd,np.tensordot(a[...,upper[0],upper[1]],spin,axes=1))
                for l,(odd,a) in vector.items()}
            components.append(matrix_values(local,[words[j] for j in chosen],frame,d,length))
        vector_values, left, right = (fmpq_mat(v.tolist()) for v in components)
        full = np.zeros((len(rows),rec["columns"]),dtype=np.int64)
        full[ids] = components[0]
        all_values.append(full)
        candidates = m.pivotcols(vector_values)
        assert len(candidates)==len(ids)
        average = (vector_values+(left+right)/(2**length))/3
        def sub(mat):
            return fmpq_mat([[mat[i,j] for j in candidates] for i in range(mat.nrows())])
        block = sub(average) * sub(vector_values).inv()
        assert block*vector_values == average and block*block == block
        for i,src in enumerate(ids):
            for j,tgt in enumerate(ids):
                p[src,tgt] = block[i,j]
    # Identical byte-level control to the archived rational basis construction.
    buf=io.BytesIO()
    np.save(buf,np.column_stack(all_values),allow_pickle=False)
    assert hashlib.sha256(buf.getvalue()).hexdigest()==data["vector_values_sha256"]
    return p


def diagonal(blocks):
    out = fmpq_mat(sum(b.nrows() for b in blocks),sum(b.ncols() for b in blocks))
    offset = 0
    for b in blocks:
        for i in range(b.nrows()):
            for j in range(b.ncols()):
                out[offset+i,offset+j] = b[i,j]
        offset += b.nrows()
    return out


def combine(a,b):
    return fmpq_mat(a.tolist()+b.tolist())


def read_map(record):
    out=fmpq_mat(*record["shape"])
    for i,row in enumerate(record["rows"]):
        for j,x in row:
            out[i,j]=x
    return out


def classical_ranks(data, n, full):
    source_offset = target_offset = 0
    ranks = []
    for r in sorted({int(key.split("_r")[1]) for key in data["bases"]}):
        sources = len(data["bases"].get(f"n{n+1}_r{r}", {}).get("rows", []))
        targets = len(data["bases"].get(f"n{n}_r{r}", {}).get("rows", []))
        classical = fmpq_mat(sources, targets)
        for i in range(sources):
            for j in range(targets):
                classical[i, j] = full[source_offset+i, target_offset+j]
        for i in range(source_offset, source_offset+sources):
            assert all(not full[i,j] for j in range(target_offset))
        if sources or targets:
            ranks.append(classical.rank())
        source_offset += sources
        target_offset += targets
    assert (source_offset, target_offset) == (full.nrows(), full.ncols())
    assert sum(ranks) == full.rank()
    return ranks


def bmn_companion(data, projectors, maps, spinors):
    # The r=0 classical maps are the upper-left blocks of the total maps.
    dims=[projectors[f"n{n}_r0"].nrows() for n in range(3)]
    q0=[fmpq_mat([[q[i,j] for j in range(dims[n])] for i in range(dims[n+1])])
        for n,q in enumerate(maps[:2])]
    pinv=projectors["n1_r0"]
    even=data["representatives"]["O18_plus_g2_O1"]["coefficients"][:dims[1]]
    co=fmpq_mat([even])
    assert co*q0[0]==fmpq_mat(1,dims[0])
    h=dims[1]-q0[0].rank()-q0[1].rank()
    trivial=pinv.rank()-(pinv*q0[0]).rank()-(q0[1]*pinv).rank()
    assert (h,trivial)==(1,0)
    assert combine(q0[1],co).rank()==q0[1].rank()+1
    assert combine(q0[1],co*pinv).rank()==q0[1].rank()

    record=data["bases"]["n1_r0"]
    words=[m.decode(w) for w in record["words"]]
    expr=Counter()
    for coefficient,row in zip(even,record["rows"],strict=True):
        for j,n in row:
            expr[words[j]]+=Fraction(coefficient)*n
    expr={w:n for w,n in expr.items() if n}
    formal=m.q(expr)
    assert formal and all(sum(map(len,w))==9 for w in formal)
    first=min(formal)

    # Evaluate the two spinor traces on one odd exterior component. Assigning
    # that component the same exterior generator for every matter fermion is
    # faithful to these degree-one products. Reflection exchanges chiralities.
    eff=sorted({l for w in expr for tr in w for l in tr})
    rng=np.random.default_rng(187231)
    up=np.triu_indices(8,1)
    colour={}
    for l in eff:
        a=np.triu(rng.integers(-1,2,(8,8),dtype=np.int64),1)
        colour[l]=a-a.T
    def value(mats):
        vals={}
        for tr in sorted({tr for w in expr for tr in w}):
            v=np.eye(8,dtype=np.int64)
            for l in tr:
                v=v@mats[l]
            vals[tr]=int(v.trace())
        return sum(c*math.prod(vals[tr] for tr in w) for w,c in expr.items())/(2**8)
    def three(col):
        reps=[{l:2*a for l,a in col.items()}]+[
            {l:np.tensordot(a[up],s,axes=1) for l,a in col.items()} for s in spinors]
        max_norm=max(int(np.abs(a).sum(axis=-1).max()) for rep in reps for a in rep.values())
        assert 8*max_norm**8<2**63
        return [value(rep) for rep in reps]
    original=three(colour)
    reflection=np.diag([-1]+[1]*7)
    reflected=three({l:reflection@a@reflection for l,a in colour.items()})
    assert reflected == [original[i] for i in (0,2,1)] and original[1]!=original[2]
    return dict(even_bmn_dimensions=dims,even_bmn_outgoing=q0[0].rank(),
                even_bmn_incoming=q0[1].rank(),even_bmn_cohomology=h,
                bmn_triality_invariant_cohomology=trivial,
                bm_even_terms=len(expr),universal_classical_image_terms=len(formal),
                universal_image_example=[m.encode(first),str(formal[first])],
                seed=187231,ordered_vector_and_spinor_evaluations=list(map(str,original)),
                odd_spinor_difference=str(original[1]-original[2]),
                reflected_difference=str(reflected[1]-reflected[2]))

def outgoing_target(data, block, spinors, incoming):
    """Project the second archived n=0, r=1 basis word-combination onto the standard representation."""
    from gaugino_primary_order import apply
    from rank_four_lowered_replay import derivative

    record = data["bases"]["n0_r1"]
    words = [m.decode(w) for w in record["words"]]
    rows = [{words[j]:n for j,n in row} for row in record["rows"]]
    index = 1
    coefficients = [int(i==index)-block[index,i] for i in range(len(rows))]
    target = Counter()
    for coefficient, row in zip(coefficients, rows, strict=True):
        for word,n in row.items():
            target[word] += Fraction(str(coefficient))*n
    target = {w:n for w,n in target.items() if n}
    assert not m.q(target)
    offset = len(data["bases"]["n0_r0"]["rows"])
    coordinate = fmpq_mat([[0]*offset+coefficients+[0]*(incoming.ncols()-offset-len(rows))])
    assert combine(incoming,coordinate).rank()==incoming.rank()+1

    # The two minus signs in the physical transformation (169) cancel.
    ranks = [len(target)]
    scalar = target
    for spin in (1,0):
        scalar = apply(scalar,derivative,(spin,0))
        ranks.append(len(scalar))
    planes = np.zeros((3,8,8), dtype=np.int64)
    for i in range(3):
        planes[i,2*i,2*i+1],planes[i,2*i+1,2*i] = 1,-1

    def scalar_value(word, matrices, replace=False):
        result = 1
        for tr in word:
            mat = np.eye(8,dtype=np.int64)
            for a,b,theta in tr:
                if theta:
                    assert len(theta)==1 and not a+b
                else:
                    assert replace and a+b==1
                mat = mat @ matrices[theta[0] if theta else 0]
            result *= int(mat.trace())
        return result

    value = sum(n*scalar_value(w,planes) for w,n in scalar.items())
    # Independent computation directly in each representation, before triality projection
    # or word canonicalization. All original eight-letter polynomials have one of each gaugino.
    u,v = np.triu_indices(8,1)
    matrices = [2*planes]+[np.tensordot(planes[:,u,v],spin,axes=1) for spin in spinors]
    probes = []
    for rep in matrices:
        total = 0
        for w,n in rows[index].items():
            order = [l for tr in w for l in tr if not l[2]]
            assert sorted(order)==[(0,1,()),(1,0,())] and sum(map(len,w))==8
            total += n*(-1)**int(order[0]==(1,0,()))*scalar_value(w,rep,replace=True)
        probes.append(Fraction(total,2**8))
    assert value==probes[0]-sum(probes)/3 and value!=0

    rref,rank = incoming.rref()
    piv = m.pivotcols(incoming)
    free = sorted(set(range(incoming.ncols()))-set(piv))
    reduced = [coordinate[0,j]-sum(coordinate[0,col]*rref[i,j] for i,col in enumerate(piv))
               for j in free]
    return dict(derivative_bidegree=[1,1],input_basis_row_zero_based=index,
                projected_cochain_coordinates=[str(x) for x in coefficients],
                quotient_free_columns=free,quotient_coordinates=list(map(str,reduced)),
                classical_terms_after_each_action=ranks,
                vacuum_planes=[[1,0,0,0],[0,1,0,0],[0,0,1,0]],
                direct_vector_and_two_spinor_traces=list(map(str,probes)),
                twice_lowered_vacuum_value=str(value))

def run(inputs):
    spinors = representations()
    data = json.loads(gzip.decompress(inputs.read_bytes()))
    by_degree = {}
    for key, rows in data["bases"].items():
        by_degree[key] = projector(key,rows,spinors)
        print(key,by_degree[key].rank(),flush=True)
    p = [diagonal([v for key,v in sorted(by_degree.items()) if key.startswith(f"n{n}_")]) for n in range(4)]
    maps = [read_map(rec) for rec in data["maps"]]
    assert [(q.nrows(),q.ncols()) for q in maps] == [(173,81),(188,173),(169,188)]
    for upper,lower in zip(maps[1:],maps):
        assert upper*lower == fmpq_mat(upper.nrows(),lower.ncols())
    for n,q in enumerate(maps):
        assert p[n+1]*q == q*p[n]
    degree = []
    for n in range(3):
        outgoing = maps[n-1] if n else fmpq_mat(81,0)
        incoming = maps[n]
        dim = p[n].nrows()-outgoing.rank()-incoming.rank()
        invariant = p[n].rank()-(p[n]*outgoing).rank()-(incoming*p[n]).rank()
        degree.append(dict(n=n,cochains=p[n].nrows(),invariant_cochains=p[n].rank(),
                           cohomology=dim,invariant_cohomology=invariant,
                           reflection_even_standard=dim-invariant))
    classes = {}
    for name,rec in data["representatives"].items():
        n=rec["n"]
        v=fmpq_mat([rec["coefficients"]])
        if n:
            assert v*maps[n-1] == fmpq_mat(1,maps[n-1].ncols())
        q=maps[n]
        def quotient_rank(x):
            return combine(q,x).rank()-q.rank()
        classes[name]=dict(whole=quotient_rank(v),invariant=quotient_rank(v*p[n]),
                           complementary=quotient_rank(v-v*p[n]))
    return dict(input_sha256=hashlib.sha256(inputs.read_bytes()).hexdigest(),
                group="SO(8)", level=18, spinor_brackets_checked=2*comb(28,2),
                ranks_of_maps=[q.rank() for q in maps],
                classical_ranks_by_grade=[classical_ranks(data,n,q) for n,q in enumerate(maps)],
                intertwining_coordinates=sum(q.nrows()*q.ncols() for q in maps),
                derivative_block_invariant_dimensions={k:v.rank() for k,v in by_degree.items()},
                reflection_even_cohomology=degree,distinguished_classes=classes,
                odd_bmn_companion=bmn_companion(data,by_degree,maps,spinors),
                outgoing_standard_target=outgoing_target(data,by_degree["n0_r1"],spinors,maps[0]))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",type=Path,default=HERE/"perturbative_complex.json.gz")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    result=run(args.inputs)
    m.write(args.out,result)
    print(result["reflection_even_cohomology"],result["distinguished_classes"])
