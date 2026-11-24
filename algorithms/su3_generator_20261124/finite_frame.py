"""Bounded integer evaluations and complete finite-colour relation matrices."""
import argparse
from collections import defaultdict
import ctypes
from fractions import Fraction
import json
from math import ceil, factorial
from pathlib import Path
import pickle
import subprocess
import time

import numpy as np
from flint import fmpq_mat, nmod_mat
from scipy.sparse import csr_matrix
from trace_algebra import skeleton, is_ferm

HERE=Path(__file__).resolve().parent
PRIME=1000003


def pivots(matrix):
    reduced,rank=matrix.rref(); answer=[]; column=0
    for i in range(rank):
        while not reduced[i,column]: column+=1
        answer.append(column); column+=1
    return answer


def library(path):
    path=Path(path)
    if not path.exists(): subprocess.run(["c++","-O3","-std=c++17","-shared","-fPIC",str(HERE/"integer_frame.cpp"),"-o",str(path)],check=True)
    lib=ctypes.CDLL(str(path)); lib.trace_frame_error.restype=ctypes.c_char_p
    lib.trace_frame_evaluate.argtypes=[ctypes.c_int]*5+[ctypes.c_void_p]*3+[ctypes.c_int]+[ctypes.c_void_p]*3
    lib.trace_frame_evaluate.restype=ctypes.c_int
    return lib


def matrices(letters,degree,points,seed):
    random=np.random.default_rng(seed); colours=3; values={}
    for letter in letters:
        odd=int(is_ferm(letter)); number=degree+1 if odd else 1
        data=np.zeros((points,number,colours,colours),np.int64)
        for matrix in data.reshape(-1,colours,colours):
            for _ in range(2):
                order=random.permutation(colours)
                if random.integers(2):
                    matrix[order,np.roll(order,1)]+=random.integers(-1,2,size=colours)
                else:
                    value=int(random.integers(-1,2)); matrix[order[0],order[0]]+=value; matrix[order[1],order[1]]-=value
        assert not data.trace(axis1=-2,axis2=-1).any()
        assert abs(data).sum(axis=-1).max(initial=0)<=2
        values[letter]=(odd,data)
    return values


def word_values(words,degree,points,seed,lib,letters=None):
    alphabet=sorted(letters if letters is not None else {letter for word in words for trace in word for letter in trace})
    assert all(sum(is_ferm(l) for t in w for l in t)==degree for w in words)
    # Every matrix factor has infinity norm at most two. The coefficient of
    # any exterior monomial is bounded by degree! times the norm product;
    # each trace contributes at most three. The bound also covers prefixes.
    bound=max((factorial(degree)*2**sum(map(len,w))*3**len(w) for w in words),default=0)
    assert bound < 2**63
    values=matrices(alphabet,degree,points,seed)
    ids={letter:i+1 for i,letter in enumerate(alphabet)}; parity=[0]; offsets=[0]; data=[]; offset=0
    for letter in alphabet:
        odd,array=values[letter]; parity.append(odd); offsets.append(offset)
        data.append(array.reshape(-1)); offset+=array.size
    parity=np.asarray(parity,np.int32); offsets=np.asarray(offsets,np.uint64); data=np.concatenate(data)
    raw=bytearray(); starts=[0]
    for word in words:
        for trace in word: raw.extend(ids[letter] for letter in trace); raw.append(0)
        starts.append(len(raw))
    starts=np.asarray(starts,np.uint64); raw=np.frombuffer(bytes(raw),np.uint8)
    output=np.empty((len(words),points*(degree+1)),np.int64)
    result=lib.trace_frame_evaluate(degree+1,points,degree,0,len(ids),parity.ctypes.data,offsets.ctypes.data,data.ctypes.data,len(words),starts.ctypes.data,raw.ctypes.data,output.ctypes.data)
    assert result==0,lib.trace_frame_error()
    assert abs(output).max(initial=0)<=bound
    return output,bound


def evaluate(rows,degree,points,seed,lib,letters=None):
    words=sorted({w for row in rows for w in row}); ids={word:i for i,word in enumerate(words)}
    entries=[]; columns=[]; starts=[0]
    for row in rows:
        for word,c in row.items(): entries.append(c); columns.append(ids[word])
        starts.append(len(entries))
    coefficients=csr_matrix((np.asarray(entries,np.int64),columns,starts),shape=(len(rows),len(words)))
    values,bound=word_values(words,degree,points,seed,lib,letters)
    row_bound=max((sum(map(abs,row.values())) for row in rows),default=0)*bound
    assert row_bound<2**63
    result=coefficients@values
    assert abs(result).max(initial=0)<=row_bound
    return result,dict(word_bound=bound,combination_bound=row_bound)


def rational_array(matrix):
    return [[str(matrix[i,j]) for j in range(matrix.ncols())] for i in range(matrix.nrows())]


def build(plan,work):
    work.mkdir(parents=True,exist_ok=True); lib=library(work.parent/"integer_frame.so")
    groups=defaultdict(list)
    for i,row in enumerate(plan["sources"]): groups[skeleton(row)].append(i)
    receipt=[]
    for index,(key,indices) in enumerate(sorted(groups.items())):
        target=work/(f"frame_{index:02d}.pkl")
        if target.exists():
            block=pickle.loads(target.read_bytes()); assert block["skeleton"]==key and block["indices"]==indices
            receipt.append((len(indices),len(block["basis_rows"]))); continue
        started=time.monotonic(); rows=[plan["sources"][i] for i in indices]
        degree=plan["length"]-3; points=ceil((len(rows)+10)/(degree+1)); seed=240511+997*plan["length"]+17191*index
        alphabet=sorted({letter for row in rows for word in row for trace in word for letter in trace})
        values,bounds=evaluate(rows,degree,points,seed,lib,alphabet)
        modular=nmod_mat(values.tolist(),PRIME); columns=pivots(modular)
        basis_rows=pivots(nmod_mat(values[:,columns].transpose().tolist(),PRIME))
        assert len(basis_rows)==len(columns)
        minor=values[np.ix_(basis_rows,columns)]
        determinant=int(nmod_mat(minor.tolist(),PRIME).det()); assert determinant
        rank=len(columns)
        print(json.dumps(dict(stage="frame",length=plan["length"],block=index,formal=len(rows),rank=rank,seconds=time.monotonic()-started)),flush=True)
        # Coordinates of every formal row, determined in an injective frame.
        # Completeness is accepted only after the total ranks attain the
        # independently known complete character dimension.
        coordinate=fmpq_mat(minor.transpose().tolist()).solve(fmpq_mat(values[:,columns].transpose().tolist())).transpose()
        assert coordinate*fmpq_mat(minor.tolist())==fmpq_mat(values[:,columns].tolist())
        rational=rational_array(coordinate)
        block=dict(skeleton=key,indices=indices,alphabet=alphabet,points=points,seed=seed,degree=degree,
                   columns=columns,basis_rows=basis_rows,determinant_mod_prime=determinant,prime=PRIME,
                   bounds=bounds,minor=minor,values=values[:,columns],coordinates=rational)
        target.write_bytes(pickle.dumps(block,protocol=5)); receipt.append((len(indices),rank))
        print(json.dumps(dict(stage="relations",length=plan["length"],block=index,relations=len(rows)-rank,
                              seconds=time.monotonic()-started,max_coordinate_length=max(len(x) for r in rational for x in r))),flush=True)
    total=sum(r for _,r in receipt); expected={6:1169,7:1896,8:1588}[plan["length"]]
    assert total==expected,(total,expected)
    summary=dict(length=plan["length"],formal=sum(n for n,_ in receipt),finite=total,relations=sum(n-r for n,r in receipt),blocks=receipt)
    (work/"receipt.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("plan",type=Path); parser.add_argument("work",type=Path)
    args=parser.parse_args(); build(pickle.loads(args.plan.read_bytes()),args.work)
