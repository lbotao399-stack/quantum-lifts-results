"""Reconstruct signed BMN primary indices over the integers with an explicit error bound."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from math import comb, isqrt
from pathlib import Path
import subprocess


# Lie algebra dimension, maximum root coordinate, sum of absolute positive root coordinates.
ROOTS = {"A2": (8,2,4), "C2": (10,2,4), "G2": (14,2,6),
         "B3": (21,1,5), "C3": (21,2,6), "B4": (36,1,7),
         "C4": (36,2,8), "D4": (28,1,6), "F4": (52,2,22),
         "D4_2": (28,2,6), "D4_3": (28,3,10), "E6_2": (78,2,22)}


def prime_moduli(grid, bound):
    modulus = 1
    candidate = ((998244353-1)//grid)*grid+1
    while modulus <= 2*bound:
        if all(candidate % k for k in range(3,isqrt(candidate)+1,2)):
            yield candidate
            modulus *= candidate
        candidate -= grid


def compute(binary, group, degree, maximum, workers=1, ordinary=False):
    dimension, coordinate, shift = ROOTS[group]
    grid = 1 << (coordinate*degree + shift).bit_length()
    bound = 6*(2**(3*dimension))*comb(degree+4*dimension+3,4*dimension+3)
    primes = list(prime_moduli(grid,bound))
    minimum = 0 if ordinary else 2
    keys = [(a,b,c) for a in range(minimum,maximum+1) for b in range(minimum,a+1)
            for c in range(minimum,b+1) if a+b+c<=degree]

    def evaluate(prime):
        command = [str(binary),group,str(degree),str(maximum),str(prime),str(grid)]
        raw = subprocess.check_output(command+(["ordinary"] if ordinary else []))
        parsed = [tuple(map(int,line.split())) for line in raw.splitlines()]
        if [row[:3] for row in parsed] != keys or any(len(row)!=4 or not 0<=row[3]<prime for row in parsed):
            raise ValueError("Incomplete modular index table")
        return prime, [row[3] for row in parsed], hashlib.sha256(raw).hexdigest()

    indices = [0]*len(keys)
    modulus = 1
    hashes = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for prime,residues,digest in pool.map(evaluate,primes):
            inverse = pow(modulus,-1,prime)
            for i,r in enumerate(residues):
                indices[i] += modulus*((r-indices[i])*inverse % prime)
            modulus *= prime
            hashes[str(prime)] = digest
    indices = [v if v<=modulus//2 else v-modulus for v in indices]
    if modulus<=2*bound or any(abs(v)>bound for v in indices):
        raise ArithmeticError("Integer reconstruction failed its bound")
    return dict(group=group,maximum_level=2*degree,maximum_first_part=maximum,
                index_type="ordinary" if ordinary else "primary",
                grid=grid,coefficient_bound=str(bound),modulus=str(modulus),
                modular_sha256=hashes,rows=[[*k,v] for k,v in zip(keys,indices)])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary",type=Path)
    p.add_argument("group",choices=ROOTS)
    p.add_argument("degree",type=int)
    p.add_argument("maximum_first_part",type=int)
    p.add_argument("--workers",type=int,default=1)
    p.add_argument("--ordinary",action="store_true")
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if not 2<=a.maximum_first_part<=a.degree or a.workers<1:
        p.error("Require a positive number of workers and 2 <= maximum first part <= degree")
    table=compute(a.binary.resolve(),a.group,a.degree,a.maximum_first_part,a.workers,a.ordinary)
    tmp=a.out.with_suffix(".part")
    tmp.write_text(json.dumps(table,separators=(",",":"))+"\n")
    tmp.replace(a.out)
    negative=sorted((sum(row[:3]),*row) for row in table["rows"] if row[3]<0)
    print(a.group,len(table["rows"]),"coefficients;",len(negative),"negative; first",negative[:1],flush=True)
