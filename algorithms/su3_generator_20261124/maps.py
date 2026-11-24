"""Exact stable and finite-colour differentials in complete native spaces."""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import pickle
import time
from trace_algebra import conjugate_word, combine


def half(row):
    return {word:Fraction(c) for word,c in row.items() if word<=conjugate_word(word)[1]}


def stable_map(source,target):
    basis=[half(row) for row in target["sources"]]; lead={max(row):i for i,row in enumerate(basis)}
    assert len(lead)==len(basis)
    answer=[]
    for index,row in enumerate(source["images"]):
        remaining=half(row); coordinates={}
        while remaining:
            word=max(remaining); i=lead[word]; coefficient=remaining[word]/basis[i][word]
            coordinates[i]=coefficient
            for key,value in basis[i].items():
                changed=remaining.get(key,0)-coefficient*value
                if changed: remaining[key]=changed
                else: remaining.pop(key,None)
        # Full oriented trace polynomial comparison, independent of dropping
        # redundant conjugate words for the triangular solve.
        assert combine(coordinates,target["sources"])==row
        answer.append(coordinates)
        if index%500==0: print(json.dumps(dict(stage="stable_map",length=source["length"],done=index,total=len(source["sources"]))),flush=True)
    return answer


def finite_coordinates(work,number):
    coordinates=[None]*number; basis=[]; offset=0
    for file in sorted(work.glob("frame_*.pkl")):
        block=pickle.loads(file.read_bytes()); size=len(block["basis_rows"])
        basis.extend(block["indices"][j] for j in block["basis_rows"])
        for original,row in zip(block["indices"],block["coordinates"]):
            coordinates[original]={offset+j:Fraction(value) for j,value in enumerate(row) if value!="0"}
        offset+=size
    assert all(row is not None for row in coordinates)
    assert all(coordinates[j]=={i:1} for i,j in enumerate(basis))
    return basis,coordinates


def sparse_combine(recipe,rows):
    answer={}
    for j,c in recipe.items():
        for k,v in rows[j].items():
            value=answer.get(k,0)+c*v
            if value: answer[k]=value
            else: answer.pop(k,None)
    return answer


def build(work,length):
    start=time.monotonic(); source=pickle.loads((work/f"stable{length}.pkl").read_bytes())
    target=pickle.loads((work/f"stable{length+1}.pkl").read_bytes())
    stable=stable_map(source,target)
    source_basis,source_coordinate=finite_coordinates(work/f"frame{length}",len(source["sources"]))
    target_basis,target_coordinate=finite_coordinates(work/f"frame{length+1}",len(target["sources"]))
    finite=[sparse_combine(stable[j],target_coordinate) for j in source_basis]
    # Every one of the complete source relations must have zero finite image.
    for j,row in enumerate(stable):
        assert sparse_combine(row,target_coordinate)==sparse_combine(source_coordinate[j],finite)
    result=dict(length=length,stable=stable,finite=finite,source_basis=source_basis,target_basis=target_basis)
    (work/f"map{length}.pkl").write_bytes(pickle.dumps(result,protocol=5))
    print(json.dumps(dict(stage="maps_complete",length=length,stable_shape=[len(stable),len(target["sources"])],finite_shape=[len(finite),len(target_basis)],stable_nonzero=sum(map(len,stable)),finite_nonzero=sum(map(len,finite)),seconds=time.monotonic()-start)),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path);parser.add_argument("length",type=int)
    args=parser.parse_args();build(args.work,args.length)
