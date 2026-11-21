"""Solve a normal-degree equation at several primes and verify the reconstructed rational primitive."""
import argparse
from fractions import Fraction
import gzip
from math import gcd, isqrt, lcm
import json
from pathlib import Path
import subprocess

from su3_chart_denominators import recover as recover_common_denominator
from su3_subregular_correction import accumulate
from su3_subregular_differential import ODD, qmonomial, variation
from su3_subregular_image import write


def equation(source, degree, previous):
    data=json.loads(gzip.decompress(source.read_bytes()))
    assert degree <= data["degree"]
    rhs={bytes.fromhex(w):Fraction(-n*data["sign"],data["denominator"])
         for w,n in data["coefficients"] if len(w)<=2*degree}
    for archive in previous:
        for word, c in json.loads(gzip.decompress(archive.read_bytes())):
            image = {w:n for w,n in qmonomial(bytes.fromhex(word)).items() if len(w)<=degree}
            accumulate(rhs, image, -Fraction(c))
    assert rhs and all(len(w)==degree for w in rhs)
    return rhs


def export(path, rhs, mixed_rules=None):
    denominator=1
    for n in rhs.values():
        denominator=lcm(denominator,n.denominator)
    common=0
    for n in rhs.values():
        common=gcd(common,int(n*denominator))
    unit=Fraction(common,denominator)
    if mixed_rules is None:
        rules=[(k,*pair,n) for k in range(128) for pair,n in variation(k).items() if len(pair)==2]
    else:
        rules=[(k,len(term),*term,n) for k,term,n in mixed_rules]
    with path.open("w") as out:
        print(*map(int,ODD),file=out)
        print(len(rules),file=out)
        for rule in rules:
            print(*rule,file=out)
        print(len(rhs),file=out)
        for w,n in sorted(rhs.items()):
            value=n/unit
            assert value.denominator==1 and abs(value)<2**63
            print(w.hex(),int(value),file=out)
    return unit


def rational(value, modulus):
    bound=isqrt((modulus-1)//2)
    x,y,previous,current=modulus,value,0,1
    while y>bound:
        quotient=x//y
        x,y=y,x-quotient*y
        previous,current=current,previous-quotient*current
    if not current or abs(current)>bound or gcd(y,current)!=1 or gcd(current,modulus)!=1:
        return None
    result=Fraction(y,current)
    assert (result.numerator-value*result.denominator)%modulus==0
    return result


def common_denominator_candidate(residues, modulus):
    # Rank candidate denominators by the height of *all* centered numerators.
    # The sparse characteristic-zero multiplication below decides whether it is a solution.
    half=modulus//2
    best=modulus
    selected=1
    for denominator in range(1,65):
        height=0
        for n in residues.values():
            height=max(height,abs((denominator*n+half)%modulus-half))
            if height>=best:
                break
        if height<best:
            best,selected=height,denominator
    if best>=modulus//256:
        return None
    return {w:Fraction((selected*n+half)%modulus-half,selected) for w,n in residues.items()}


def check_candidate(values, unit, rhs, degree, image_of=None):
    candidate={w:v*unit for w,v in values.items() if v}
    check=dict(rhs)
    for w,c in candidate.items():
        image = image_of(w) if image_of is not None else {v:n for v,n in qmonomial(w).items() if len(v)==degree}
        accumulate(check,image,-c)
    return candidate if not check else None


def solve(engine, work, rhs, output, degree, depth, *, mixed_rules=None, image_of=None,
          minimum_roots=0, maximum_roots=None):
    if (mixed_rules is None) != (image_of is None):
        raise ValueError("Mixed Taylor rules require their rational replay map")
    if not 0 <= minimum_roots <= 3 or (minimum_roots and mixed_rules is None):
        raise ValueError("Restricting off-diagonal ghost degree requires the mixed chart")
    if maximum_roots is not None and (mixed_rules is None or not minimum_roots <= maximum_roots <= 4):
        raise ValueError("Invalid off-diagonal ghost quotient")
    work.mkdir(parents=True,exist_ok=True)
    matrix=work/"equation.txt"
    unit=export(matrix,rhs,mixed_rules)
    residues={}
    modulus=1
    records=[]
    for prime in (2147483647,1000000007,1000000009,998244353,1000033,1000003):
        path=work/f"solution-{prime}.txt"
        command=[str(engine.resolve()),str(matrix),str(path),str(prime),str(depth)]
        if mixed_rules is not None:
            command.append("mixed")
        if minimum_roots or maximum_roots is not None:
            command.append(str(minimum_roots))
        if maximum_roots is not None:
            command.append(str(maximum_roots))
        subprocess.run(command,check=True)
        with path.open() as stream:
            dimensions=list(map(int,next(stream).split()))
            new={bytes.fromhex(w):int(n) for w,n in (line.split() for line in stream)}
        old=modulus
        inv=pow(old,-1,prime)
        for w in set(new)|residues.keys():
            x=residues.get(w,0)
            residues[w]=x+old*((new.get(w,0)-x)*inv%prime)
        modulus*=prime
        records.append(dict(prime=prime,equations=dimensions[0],parents=dimensions[1],
                            peeled=dimensions[2],other_pivots=dimensions[3]))
        values=common_denominator_candidate(residues,modulus)
        candidate=check_candidate(values,unit,rhs,degree,image_of) if values is not None else None
        if candidate is None:
            values={w:rational(n,modulus) for w,n in residues.items()}
            if all(v is not None for v in values.values()):
                candidate=check_candidate(values,unit,rhs,degree,image_of)
        if candidate is None and len(records) >= 2:
            recovered = recover_common_denominator(residues, modulus, unit, rhs, degree, image_of)
            if recovered is not None:
                candidate, denominator, trials = recovered
                records[-1].update(verified_common_denominator=denominator, exact_denominator_trials=trials)
        if candidate is not None:
            write(output,[[w.hex(),str(c)] for w,c in sorted(candidate.items())])
            return dict(input_monomials=len(rhs),rational_unit=str(unit),terms=len(candidate),residual=0,systems=records)
    raise ValueError("the modular solutions did not reconstruct a verified rational solution")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image",type=Path,required=True)
    p.add_argument("--degree",type=int,required=True)
    p.add_argument("--initial",type=Path,action="append",default=[])
    p.add_argument("--engine",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--recipe",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--depth",type=int,default=1)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires assertions")
    rhs=equation(args.image,args.degree,args.initial)
    result=solve(args.engine,args.work,rhs,args.recipe,args.degree,args.depth)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
