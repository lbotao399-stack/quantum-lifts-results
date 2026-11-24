"""Extract a closed nonproduct representative and its full rational dual."""
import argparse
from fractions import Fraction
from math import gcd,lcm
from pathlib import Path
import json,pickle,time
from flint import fmpq,fmpq_mat,nmod_mat
from maps import stable_map,finite_coordinates,sparse_combine
from finite_frame import pivots,PRIME


def modular(rows,width):
    result=nmod_mat(len(rows),width,PRIME)
    for i,row in enumerate(rows):
        for j,c in row.items():result[i,j]=c.numerator*pow(c.denominator,-1,PRIME)%PRIME
    return result


def rational(rows,width):
    result=fmpq_mat(len(rows),width)
    for i,row in enumerate(rows):
        for j,c in row.items():result[i,j]=fmpq(c.numerator,c.denominator)
    return result


def rational_reconstruct(value,prime=PRIME):
    if value==0:return Fraction(0)
    bound=int((prime//2)**0.5)
    r0,r1=prime,int(value);t0,t1=0,1
    while abs(r1)>bound:r0,r1,t0,t1=r1,r0-(r0//r1)*r1,t1,t0-(r0//r1)*t1
    if not t1 or abs(t1)>bound or gcd(r1,t1)!=1:raise ValueError("reconstruction failed")
    result=Fraction(r1,t1)
    assert result.numerator*pow(result.denominator,-1,prime)%prime==value
    return result


def build(work):
    started=time.monotonic();incoming=pickle.loads((work/"map6.pkl").read_bytes());outgoing=pickle.loads((work/"map7.pkl").read_bytes())
    source=pickle.loads((work/"stable7.pkl").read_bytes());products=pickle.loads((work/"products.pkl").read_bytes())
    _,coordinate=finite_coordinates(work/"frame7",4991)
    stable_products=stable_map(dict(length="products",sources=products["complete"],images=products["complete"]),source)
    stable_proper=stable_map(dict(length="proper",sources=products["proper"],images=products["proper"]),source)
    finite_products=[sparse_combine(row,coordinate) for row in stable_products]
    finite_proper=[sparse_combine(row,coordinate) for row in stable_proper]
    # Full stable identities plus the already verified commutative finite
    # coordinate maps certify the finite identities without a dense product.
    for row in stable_products:assert not sparse_combine(row,outgoing["stable"])
    for row in incoming["stable"]:assert not sparse_combine(row,outgoing["stable"])
    print(json.dumps(dict(stage="cycles_and_nilpotence",seconds=time.monotonic()-started)),flush=True)
    matrix=modular(outgoing["finite"],1588).transpose();reduced,rank=matrix.rref();assert rank==1037
    pivot=[];column=0
    for i in range(rank):
        while not reduced[i,column]:column+=1
        pivot.append(column);column+=1
    free=[j for j in range(1896) if j not in pivot]
    kernels=[]
    for j in free:
        row={j:Fraction(1)}
        row.update({p:Fraction(-int(reduced[i,j])%PRIME) for i,p in enumerate(pivot) if reduced[i,j]})
        kernels.append(row)
    base=incoming["finite"]+finite_proper
    all_cycles=modular(base+kernels,1896).transpose();selected=pivots(all_cycles)
    assert len(selected)==859
    outside=[j-len(base) for j in selected if j>=len(base)];assert len(outside)==1
    chosen=outside[0]
    try:candidate={j:rational_reconstruct(int(c)%PRIME) for j,c in kernels[chosen].items()}
    except ValueError:candidate=None
    if candidate is None or sparse_combine(candidate,outgoing["finite"]):
        independent_equations=pivots(matrix.transpose())
        square=fmpq_mat([[fmpq(outgoing["finite"][p].get(i,Fraction(0)).numerator,outgoing["finite"][p].get(i,Fraction(0)).denominator) for p in pivot] for i in independent_equations])
        right=fmpq_mat([[fmpq(-outgoing["finite"][free[chosen]].get(i,Fraction(0)).numerator,outgoing["finite"][free[chosen]].get(i,Fraction(0)).denominator)] for i in independent_equations])
        solved=square.solve(right);candidate={free[chosen]:Fraction(1)}
        candidate.update({p:Fraction(str(solved[i,0])) for i,p in enumerate(pivot) if solved[i,0]})
    assert not sparse_combine(candidate,outgoing["finite"])
    print(json.dumps(dict(stage="representative",terms=len(candidate),free_column=free[chosen],seconds=time.monotonic()-started)),flush=True)
    independent=[base[j] for j in selected if j<len(base)]
    assert len(independent)==858
    spanning=independent+[candidate];columns=pivots(modular(spanning,1896))
    assert len(columns)==859
    square=fmpq_mat([[fmpq(row.get(j,Fraction(0)).numerator,row.get(j,Fraction(0)).denominator) for j in columns] for row in spanning])
    right=fmpq_mat([[0]]*858+[[1]]);solved=square.solve(right)
    dual={j:Fraction(str(solved[i,0])) for i,j in enumerate(columns) if solved[i,0]}
    pairing=lambda row:sum(c*dual.get(j,0) for j,c in row.items())
    assert all(pairing(row)==0 for row in base)
    assert pairing(candidate)==1
    assert all(pairing(row)==0 for row in finite_products)
    result=dict(candidate=candidate,dual=dual,stable_products=stable_products,stable_proper=stable_proper,
                finite_products=finite_products,finite_proper=finite_proper,selected=selected,dual_columns=columns,
                outgoing_rank=rank,free_column=free[chosen])
    (work/"generator.pkl").write_bytes(pickle.dumps(result,protocol=5))
    print(json.dumps(dict(stage="generator_complete",representative_terms=len(candidate),dual_terms=len(dual),seconds=time.monotonic()-started)),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("work",type=Path);args=parser.parse_args();build(args.work)
