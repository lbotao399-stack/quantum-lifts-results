"""Exact SO(9) pure-gaugino boundary rank at length twelve and the product source."""
import argparse
import math
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np
from flint import fmpz_mat

import mixed_incoming as m
from fortuity.hwcoh.evalg import table
from mixed33_exact import read
from rank_four_s5 import decode_rows, dimensions
from rank_four_source import exact_quantum
from q1d_exact import pair_terms


class IntegralFrame(m.Frame):
    def necklace_values(self, necks):
        result, stack, previous = {}, [], ()
        for tr in sorted(set(necks)):
            prefix = 0
            while prefix < min(len(previous),len(tr)) and previous[prefix] == tr[prefix]:
                prefix += 1
            del stack[prefix:]
            if not stack:
                degree,mat = self.L[tr[0]]
                stack.append((mat,degree))
            for letter in tr[len(stack):]:
                mat,degree = stack[-1]
                stack.append(self.mat_mul_letter(mat,degree,letter))
            mat,degree = stack[-1]
            result[tr] = mat.trace(axis1=2,axis2=3),degree
            previous = tr
        return result

    def mat_mul_letter(self, arr, degree, letter):
        parity, right = self.L[letter]
        assert parity == 1
        i,j,sign,st = table(self.k,degree,1)
        product = (arr[:,i] @ right[:,j]) * sign[None,:,None,None]
        return np.add.reduceat(product,st,axis=1), degree+1

    def scal_mul(self, left, dl, right, dr):
        i,j,sign,st = table(self.k,dl,dr)
        return np.add.reduceat(left[:,i]*right[:,j]*sign[None,:],st,axis=1), dl+dr


G = ((1,0,()),(0,1,()))


def pfaffian(matrix):
    @lru_cache(None)
    def value(ids):
        if not ids:
            return 1
        return sum((-1)**(j-1)*int(matrix[ids[0],ids[j]])*value(ids[1:j]+ids[j+1:])
                   for j in range(1,len(ids)))
    return value(tuple(range(len(matrix))))


def sample(task):
    seed, words, index, orientation = task
    frame = IntegralFrame(4,1,12,[12],G,seed=1,group="SO_odd")
    rng = np.random.default_rng(seed)
    for l in G:
        blocks = np.zeros((1,12,9,9),dtype=np.int64)
        for mat in blocks[0]:
            perm = rng.permutation(9)
            for i,j in zip(perm[:8:2],perm[1:8:2]):
                s = 2*int(rng.integers(2))-1
                mat[i,j],mat[j,i] = s,-s
        frame.L[l] = 1,blocks
    values = frame.eval_words(words, dtype=np.int64)[:,0]
    gram = np.einsum("iab,jba->ij",frame.L[G[0]][1][0],frame.L[G[1]][1][0])
    direct = math.factorial(6)*pfaffian(gram-gram.T)
    assert orientation*int(values[index]) == direct
    return values.tolist()


def unpack(path):
    data = read(path)
    letters = [(i,j,tuple(theta)) for i,j,theta in data["letters"]]
    result = {}
    for raw,n in data["words"]:
        encoded = bytes.fromhex(raw)
        assert encoded[-1]==0
        word = tuple(tuple(letters[i-1] for i in trace) for trace in encoded[:-1].split(b"\x00"))
        assert word not in result and m.so_canon(word) == (1,word)
        result[word] = Fraction(n,data["denominator"])
    return result


def times_scalar(rows):
    out = {}
    for w,n in rows.items():
        sign,word = m.so_canon((G,)+w)
        if sign:
            out[word] = out.get(word,0)+sign*n
    return {w:n for w,n in out.items() if n}


def run(source, primitive, out, workers):
    mu = (6,6,0,0,0)
    dim = dimensions("SO_odd",mu,12)[12]
    pure = decode_rows(m.prepare(mu,12),"sources")
    lower = m.prepare(mu,11)
    bd = decode_rows(lower,"images")
    sign,w = m.so_canon((G,)*6)
    words = sorted({term for row in pure+bd for term in row} | {w})
    # Each exterior coefficient matrix is a signed partial permutation.
    # A product has entries of magnitude at most one and trace at most nine.
    bound = math.factorial(12)*9**6
    assert bound < 2**63 and all(sum(map(len,v))==12 and len(v)<=6 for v in words)
    seeds = list(range(9321601,9321601+dim+8))
    with ProcessPoolExecutor(workers) as pool:
        values = list(pool.map(sample,[(s,words,words.index(w),sign) for s in seeds]))
    evaluations = fmpz_mat(values).transpose()
    def evaluated(rows):
        return fmpz_mat([[r.get(v,0) for v in words] for r in rows])*evaluations
    basis,images,joint = (evaluated(rs) for rs in (pure, bd, bd+[{w:sign}]))
    assert basis.rank()==dim
    rank,aug = images.rank(),joint.rank()
    assert (rank,aug)==(dim-1,dim)
    f,h = unpack(source),unpack(primitive)
    s5,w5 = m.so_canon((G,)*5)
    qf = exact_quantum(f)
    rhs = m.q(h)
    rhs[w5] = rhs.get(w5,0)+s5
    assert qf == rhs
    letters = {l for v in f for tr in v for l in tr}
    assert all(not pair_terms(l,g) and not pair_terms(g,l) for l in letters|set(G) for g in G)
    product = times_scalar(f)
    image = times_scalar(qf)
    subtract = m.q(times_scalar(h))
    subtract[w] = subtract.get(w,0)+sign
    assert image == subtract
    cert=dict(group="SO(9)",weight=mu,source_weight=[5,5,1,1,1],source_length=12,
              terminal_dimension=dim,lower_dimension=dimensions("SO_odd",mu,11)[11],
              formal_terminal=len(pure),formal_lower=len(bd),formal_words=len(words),
              seeds=seeds,integer_bound=bound,exact_probe_rank=basis.rank(),boundary_rank=rank,
              with_S6_rank=aug,independent_pfaffians=len(seeds),product_source_terms=len(product),
              first_action=[[[m.encode(v),str(n)] for v,n in sorted(d.items())] for d in (qf,image)],
              product_primitive_terms=len(times_scalar(h)),
              boundary_evaluations=[[int(x) for x in row] for row in images.tolist()],
              S6_evaluations=[int(x) for x in joint.tolist()[-1]])
    m.write(out,cert)
    print(dim,rank,aug,"product source",len(product),"quantum terms",len(image),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--primitive",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    assert __debug__
    run(**vars(p.parse_args()))
