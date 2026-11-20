"""Bound the part of the F4 first-order source left by universal classical cycles."""
import argparse
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from flint import nmod_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from f4_colour import basis, frame
from f4_singlet_characters import singlets
from mixed33_exact import read
from mixed_middle import build, pin


def coefficient_matrix(rows, width, prime):
    ri, ci, values = [], [], []
    for j,row in enumerate(rows):
        for i,n in row:
            ri.append(j)
            ci.append(i)
            values.append(n % prime)
    mat = csr_matrix((np.array(values, dtype=np.int64), (ri,ci)),shape=(len(rows),width))
    assert int(max(np.diff(mat.indptr)))*(prime-1)**2 < 2**63
    return mat


def formal_columns(rows, prime):
    pivot = {}
    for row in rows:
        x = {i:v % prime for i,v in row if v % prime}
        while x:
            j = min(x)
            c = x[j]
            if j not in pivot:
                inv = pow(c,-1,prime)
                pivot[j] = [(i,n*inv % prime) for i,n in x.items()]
                break
            for i,n in pivot[j]:
                v = (x.get(i,0)-c*n) % prime
                if v:
                    x[i] = v
                else:
                    x.pop(i,None)
    return sorted(pivot)


def probe(rows, words, dims, prime, lib, workers, salt, colour_frame=frame):
    co = coefficient_matrix(rows,len(words),prime)
    degrees = [m.nferm(w) for w in words]
    assert set(degrees) == set(dims)
    blocks, frames = [], []
    for degree,dimension in sorted(dims.items()):
        idx = [i for i,n in enumerate(degrees) if n==degree]
        ws = [words[i] for i in idx]
        letters = {l for w in ws for tr in w for l in tr}
        k = degree+2
        pts = math.ceil((dimension+18)/math.comb(k,degree))
        batches = []
        for start in range(0,pts,4):
            count = min(4,pts-start)
            seed = salt+137*degree+start*187
            f = colour_frame(letters,count,k,degree,prime,seed)
            chunks = np.array_split(np.arange(len(ws)),min(workers*2,len(ws)))
            with ThreadPoolExecutor(workers) as pool:
                values = np.concatenate(list(pool.map(
                    lambda jj:m.values([ws[i] for i in jj],f,lib),chunks)))
            tests = np.linspace(0,len(ws)-1,3,dtype=int)
            np.testing.assert_array_equal(values[tests],f.eval_words([ws[i] for i in tests],dtype=np.int64))
            batches.append(np.asarray(co[:,idx] @ values) % prime)
            frames.append(dict(degree=degree,k=k,points=count,seed=seed,words=len(ws)))
        blocks.append(np.concatenate(batches,axis=1))
        print("probes",degree,len(ws),pts,"dimension",dimension,flush=True)
    return np.concatenate(blocks,axis=1), frames


def run(plan_path, work, workers, prime):
    work.mkdir(exist_ok=True,parents=True)
    z = read(plan_path)
    assert z["mu"] == [4,4,1,1,1] and z["length"] == 10
    # All 11 terminal letters have weight one. Independently enumerate the
    # formal terminal highest weights to bound the rational formal Q0 rank.
    terminal = work/"terminal"
    terminal.mkdir(exist_ok=True)
    fplan = read(build(4,terminal,workers,length=11,power=5))
    assert all(not v for v in fplan["images"])
    formal_dimension = fplan["formal_singlets"]
    words = list(map(m.decode,z["target_dictionary"]))
    formal_ids = formal_columns(z["images"],prime)
    assert len(formal_ids) == formal_dimension == 358
    formal = coefficient_matrix(z["images"],len(words),prime)[:,formal_ids].toarray()
    group_weight = tuple(z["mu"])
    dims_out,dims_in = (singlets("F4",group_weight,n) for n in (11,10))
    lib = m.library(work,27,prime)
    images, recipes_a = probe(z["images"],words,dims_out,prime,lib,workers,507100)
    classic = m.minor(images,prime)
    print("classical image",classic["rank"],"of",sum(dims_out.values()),flush=True)
    src, recipes_b = probe(z["sources"],list(map(m.decode,z["source_dictionary"])),
                           dims_in,prime,lib,workers,507260)
    srcminor = m.minor(src,prime)
    joint = m.minor(np.column_stack((src,formal)),prime)
    kernel_upper = sum(dims_in.values())-classic["rank"]
    stable_lower = joint["rank"]-formal_dimension
    assert 0 <= stable_lower <= kernel_upper
    mats,idx = basis()
    m.write(work/"ranks.json",dict(group="F4",weight=group_weight,length=z["length"],prime=prime,
            plan_sha256=pin(plan_path),formal_terminal_singlets=formal_dimension,
            formal_image_word_columns=formal_ids,source_character=dims_in,target_character=dims_out,
            classical=classic,source=srcminor,with_formal_image=joint,
            kernel_upper=kernel_upper,restricted_formal_kernel_lower=stable_lower,
            unaccounted_kernel_upper=kernel_upper-stable_lower,frames=recipes_a+recipes_b,
            derivation_commutator_indices=idx,basis_sha256=m.digest(mats.astype("<i2").tobytes())))
    print("source",srcminor["rank"],"joint",joint["rank"],"classical kernel upper",kernel_upper,
          "universal lower",stable_lower,"remaining upper",kernel_upper-stable_lower,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=12)
    p.add_argument("--prime",type=int,default=1000003)
    assert __debug__
    run(**vars(p.parse_args()))
