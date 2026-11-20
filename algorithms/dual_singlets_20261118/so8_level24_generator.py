"""Separate the odd triality image of the level-24 BMN class from all classical boundaries."""
import argparse
import gzip
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat

import mixed_incoming as m
from f4_bmn_primary import LETTERS, SIGNS, WEIGHT, bmn_convention
from fortuity.full.letters import is_ferm
from so8_triality import representations
from terminal_interval_ranks import IntegralTerminalFrame
from trace_descent import dimensions, highest_weights, multiply_rows


def even_class(path):
    result = Counter()
    for coefficient, traces in json.loads(gzip.decompress(path.read_bytes())):
        conversion = math.prod(SIGNS[i] for tr in traces for i in tr)
        sign, word = m.so_canon(tuple(tuple(LETTERS[i] for i in tr) for tr in traces))
        result[word] += sign * conversion * coefficient
    return {w:n for w,n in result.items() if n}


def evaluate(rows, degree, seed, points, spinors=(), operator=None):
    words = sorted({w for row in rows for w in row})
    length, = {sum(map(len,w)) for w in words}
    assert all(m.nferm(w)==degree for w in words)
    # The exact vector matrices have at most seven nonzero unit entries per row.
    assert math.factorial(degree)*7**length*8**max(map(len,words)) < 2**63
    slots = degree+2
    upper = np.triu_indices(8,1)
    pieces, chiral = [], [[],[]]
    distinguished = sorted(operator) if operator else []
    for start in range(0,points,6):
        count = min(6,points-start)
        rng = np.random.default_rng(seed+7919*start)
        frame = IntegralTerminalFrame(8,count,slots,[degree],[])
        for l in LETTERS:
            odd = int(is_ferm(l))
            raw = rng.integers(-1,2,(count,slots if odd else 1,8,8),dtype=np.int64)
            a = np.triu(raw,1)
            frame.L[l] = odd,a-a.swapaxes(-1,-2)
        pieces.append(multiply_rows(rows,words,frame.eval_words(words,dtype=np.int64)))
        for i,spin in enumerate(spinors):
            # The spin matrices represent twice the generators; divide before evaluation.
            mod = m.Frame(8,count,slots,[degree],[],p=m.P)
            for l,(odd,a) in frame.L.items():
                mod.L[l] = odd,np.tensordot(a[...,upper[0],upper[1]],spin,axes=1)*pow(2,-1,m.P)%m.P
            vals = mod.eval_words(distinguished,dtype=np.int64)
            chiral[i].append(multiply_rows([operator],distinguished,vals).ravel()%m.P)
    values = np.column_stack(pieces)
    return values, [np.concatenate(v) for v in chiral if v]


def combine(rows, coefficients):
    result = Counter()
    for row,n in zip(rows,coefficients,strict=True):
        if n:
            for w,v in row.items():
                result[w] += int(n)*v
    return {w:n for w,n in result.items() if n}


def small_rational(x, prime):
    if not x:
        return Fraction(0)
    bound = math.isqrt(prime//2)
    a,b,u,v = prime,x,0,1
    while abs(b)>bound:
        k=a//b
        a,b,u,v=b,a-k*b,v,u-k*v
    if v and abs(v)<=bound and math.gcd(b,v)==1:
        return Fraction(b,v)
    return None


def relation_rows(values, square, pivots, columns, free):
    p=m.P
    target=values[np.ix_(free,columns)].T
    coefficients=square.transpose().solve(nmod_mat((target%p).tolist(),p)).transpose()
    out=[]
    exact=None
    for j,idx in enumerate(free):
        lifted=[small_rational(int(x),p) for x in coefficients.tolist()[j]]
        if all(x is not None for x in lifted):
            den=math.lcm(*(x.denominator for x in lifted))
            rel=[0]*len(values)
            rel[idx]=den
            for i,v in zip(pivots,lifted):
                rel[i]=-int(v*den)
            if not np.any(np.asarray(rel,dtype=object)@values.astype(object)):
                out.append(rel)
                continue
        if exact is None:
            exact=fmpz_mat(values[np.ix_(pivots,columns)].T.tolist()).solve(fmpz_mat(target.tolist())).transpose()
        rational=[Fraction(str(exact[j,i])) for i in range(len(pivots))]
        den=math.lcm(*(x.denominator for x in rational))
        rel=[0]*len(values)
        rel[idx]=den
        for i,v in zip(pivots,rational):
            rel[i]=-int(v*den)
        assert not np.any(np.asarray(rel,dtype=object)@values.astype(object))
        out.append(rel)
    return out


def run(path,seed):
    p=m.P
    bmn_convention()
    op=even_class(path)
    assert len(op)==13449 and all(m.nferm(w)==2 and sum(map(len,w))==10 for w in op)
    assert all(not combine([m.raise_word(w,k) for w in op],op.values()) for k in (1,2,3))
    by_degree=highest_weights(WEIGHT,10)
    finite=dimensions(WEIGHT,10,8,"SO")
    stable=dimensions(WEIGHT,10,22,"SO")
    assert stable=={d:len(v) for d,v in by_degree.items()}=={0:36,2:523}
    assert finite=={0:36,2:516}
    # No scalar trace identity in this singlet; a two-fermion polynomial lifts
    # universally, if at all, within the two-fermion formal basis.
    rows=by_degree[2]
    source,chiral=evaluate(rows+[op],2,seed,92,representations(),op)
    table,vector=source[:-1],source[-1]
    mod=nmod_mat((table%p).tolist(),p)
    columns=m.pivotcols(mod)
    pivots=m.pivotcols(mod.transpose())
    assert len(columns)==len(pivots)==finite[2]
    free=sorted(set(range(len(rows)))-set(pivots))
    square=nmod_mat((table[np.ix_(pivots,columns)]%p).tolist(),p)
    relations=relation_rows(table,square,pivots,columns,free)
    images=[m.q(combine(rows,rel)) for rel in relations]
    # A representative of 2 O - rho O - rho^{-1} O on this integral colour frame.
    projected=(2*(vector%p)-chiral[0]-chiral[1])%p
    particular=square.transpose().solve(nmod_mat([[int(projected[j])] for j in columns],p))
    co=[0]*len(rows)
    for k,i in enumerate(pivots):
        co[i]=int(particular[k,0])
    assert np.array_equal(np.asarray(co,dtype=object)@(table%p).astype(object)%p,projected)
    obstruction=m.q(combine(rows,co))
    words=sorted({w for row in images+[obstruction] for w in row})
    matrix=[[row.get(w,0) for w in words] for row in images]
    rank=fmpz_mat(matrix).rank()
    augmented=np.array([[row.get(w,0)%p for w in words] for row in images+[obstruction]],dtype=np.int64)
    minor=m.minor(augmented,p)
    assert rank==5 and minor["rank"]==6 and minor["rows"][-1]==len(images)
    q0=m.q(op)
    targets=highest_weights(WEIGHT,11)[1]
    target,_=evaluate(targets+[q0],1,seed+4099,115)
    target_rank=nmod_mat((target[:-1]%p).tolist(),p).rank()
    assert target_rank==dimensions(WEIGHT,11,8,"SO")[1]==325 and not np.any(target[-1])
    i=int(np.flatnonzero(chiral[0]!=chiral[1])[0])
    return dict(level=24,weight=WEIGHT,source_length=10,spinor_generators_checked=2*math.comb(28,2),
                seed=seed,prime=p,operator_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                terms=len(op),universal_variation_terms=len(q0),classical_target_rank=target_rank,
                nonzero_classical_image_coordinates=int(np.count_nonzero(target[-1])),
                stable_dimensions=stable,eight_colour_dimensions=finite,
                frame_determinant=int(square.det()),relation_free_rows=free,
                relations=[[ [i,n] for i,n in enumerate(row) if n] for row in relations],
                rational_relation_image_rank=rank,augmented_modular_rank=minor["rank"],
                minor_rows=minor["rows"],minor_words=[m.encode(words[j]) for j in minor["columns"]],
                minor_determinant=minor["determinant"],
                first_odd_probe_column=i,probe_vector_spinors=[int(vector[i]%p),int(chiral[0][i]),int(chiral[1][i])])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--operator",type=Path,default=Path(__file__).with_name("so8_level24")/"even_bmn.json.gz")
    parser.add_argument("--seed",type=int,default=240088)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    data=run(args.operator,args.seed)
    m.write(args.out,data)
    print(data["rational_relation_image_rank"],data["augmented_modular_rank"],data["minor_determinant"])
