"""Project the two spacetime Taylor directions on the SU(3) scalar chart and replay their signs."""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import gcd
from pathlib import Path

from su3_subregular_chart_replay import verify
from su3_subregular_differential import ODD, variation
from su3_subregular_image import write
from su3_subregular_root_coupled import read_coefficients
from su3_subregular_root_linear import decode
from su3_subregular_sparse import export, rational


REFERENCE = Path(__file__).resolve().parent / "taylor_hooks/su3_scalar_levi/subregular"
EXCHANGED = bytes((k & ~12) | ((k & 4) << 1) | ((k & 8) >> 1) for k in range(256))


@lru_cache(maxsize=50000)
def exchange(word):
    mapped = word.translate(EXCHANGED)
    odd = [k for k in mapped if ODD[k]]
    sign = (-1)**sum(k > other for i,k in enumerate(odd) for other in odd[i+1:])
    return sign, bytes(sorted(mapped))


def add(poly, word, coefficient):
    value = poly.get(word,0) + coefficient
    if value:
        poly[word] = value
    else:
        poly.pop(word,None)


def exchange_poly(poly):
    return {new:sign*c for w,c in poly.items() for sign,new in [exchange(w)]}


def quotient(poly):
    result = {}
    for w,c in poly.items():
        sign,other = exchange(w)
        if w == other and sign == 1:
            continue
        add(result,min(w,other),c if w <= other else -sign*c)
    return result


def expand(poly):
    result = {}
    for w,c in poly.items():
        sign,other = exchange(w)
        half = Fraction(c,2)
        add(result,w,half)
        add(result,other,-sign*half)
    return result


def check_rules():
    tested = 0
    for k in range(128):
        if k >= 64 and k % 16 == 0:
            continue
        if EXCHANGED[EXCHANGED[k]] != k or ODD[k] != ODD[EXCHANGED[k]]:
            raise ValueError("Invalid exchange of jet labels")
        if exchange_poly(variation(k)) != variation(EXCHANGED[k]):
            raise ValueError("Exchange does not commute with the classical supercharge")
        tested += 1
    return tested


def prepare(level, output):
    output.mkdir(parents=True,exist_ok=True)
    directory = REFERENCE/"complete"
    image = directory/f"j{level}_image.json.gz"
    data = decode(image)
    coefficients = {bytes.fromhex(w):n for w,n in data["coefficients"]}
    if any(coefficients.get(other,0) != -sign*n for w,n in coefficients.items() for sign,other in [exchange(w)]):
        raise ArithmeticError("The complete first image has a spin-symmetric component")
    qimage = quotient(coefficients)
    if expand(qimage) != coefficients:
        raise ArithmeticError("The orbit basis failed to reconstruct the complete first image")
    image_size, orbit_size = len(coefficients), len(qimage)
    del coefficients,qimage,data
    high = decode(directory/f"j{level}_high_d_primitive.json.gz")
    residual_path = output/"high_d_residual.json.gz"
    receipt = verify(image,high,4,1,REFERENCE,residual_path=residual_path)
    if json.loads(json.dumps(receipt)) != json.loads((directory/f"j{level}_high_d.json").read_bytes()):
        raise ArithmeticError("The high-d primitive did not reproduce its reference")
    data = decode(residual_path)
    rhs = {bytes.fromhex(w):Fraction(-n,data["denominator"]) for w,n in data["coefficients"]}
    if any(any(k//16==5 for k in w) for w in rhs):
        raise ArithmeticError("The remaining equation contains d coefficients")
    rules = [(k,term,n) for k in range(128) if k//16!=5 for term,n in variation(k).items()]
    unit = export(output/"equation.txt",rhs,rules)
    info = dict(level=level,intertwined_generators=check_rules(),first_image_terms=image_size,
                first_image_orbits=orbit_size,remaining_terms=len(rhs),remaining_orbits=len(quotient(rhs)),
                equation_unit=str(unit))
    (output/"spin.json").write_text(json.dumps(info,indent=2)+"\n")
    return info


def replay(level, recipe, out):
    directory = REFERENCE/"complete"
    high = quotient(read_coefficients(decode(directory/f"j{level}_high_d_primitive.json.gz")))
    low = read_coefficients(decode(recipe) if isinstance(recipe,Path) else recipe)
    if quotient(low) != low or any(any(k//16==5 for k in w) for w in low):
        raise ValueError("Expected the alternating orbit coordinates without d spectators")
    for w,c in low.items():
        add(high,w,c)
    polynomial = [[w.hex(),str(c)] for w,c in sorted(expand(high).items())]
    info = verify(directory/f"j{level}_image.json.gz",polynomial,4,0,REFERENCE)
    if info["residual_terms_below_d_count"]:
        raise ArithmeticError("A full chart coefficient survived")
    write(out,polynomial)
    return info


def reconstruct(level, inputs, preparation, out):
    info = json.loads(preparation.read_bytes())
    if info["level"] != level or info["intertwined_generators"] != check_rules():
        raise ValueError("The prepared chart differs from the requested one")
    unit = Fraction(info["equation_unit"])
    modulus = 1
    residues = {}
    systems = []
    for prime_text, filename in inputs:
        prime = int(prime_text)
        path = Path(filename)
        if prime < 3 or gcd(modulus*unit.denominator,prime) != 1:
            raise ValueError("Incompatible reconstruction moduli")
        with path.open() as stream:
            header = list(map(int,next(stream).split()))
            if len(header) != 4 or not (1 <= header[2]+header[3] <= min(header[:2])):
                raise ValueError("Invalid sparse solve dimensions")
            values = {}
            for line in stream:
                token, value = line.split()
                word, coefficient = bytes.fromhex(token), int(value)
                if word in values or not 0 < coefficient < prime:
                    raise ValueError("Invalid modular coordinate")
                values[word] = coefficient
        inv = pow(modulus,-1,prime)
        for w in residues.keys() | values.keys():
            old = residues.get(w,0)
            residues[w] = old + modulus*((values.get(w,0)-old)*inv % prime)
        modulus *= prime
        systems.append(dict(prime=prime, header=header, nonzero=len(values),
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))

    factor = unit.numerator * pow(unit.denominator,-1,modulus)
    coordinates = []
    for w, value in sorted(residues.items()):
        coefficient = rational(value*factor % modulus, modulus)
        if coefficient is None:
            raise ArithmeticError("More modular coordinates are required for rational reconstruction")
        if coefficient:
            coordinates.append([w.hex(),str(coefficient)])
    result = replay(level,coordinates,out)
    return dict(level=level, equation_unit=str(unit), modulus=str(modulus), systems=systems,
                reconstructed_orbit_coefficients=len(coordinates), replay=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level",type=int,choices=(32,34),required=True)
    parser.add_argument("--prepare",type=Path)
    parser.add_argument("--solution",type=Path)
    parser.add_argument("--native",nargs=2,action="append",metavar=("PRIME","FILE"))
    parser.add_argument("--preparation",type=Path)
    parser.add_argument("--receipt",type=Path)
    parser.add_argument("--out",type=Path)
    args = parser.parse_args()
    if (not __debug__ or sum(map(bool,(args.prepare,args.solution,args.native))) != 1
            or bool(args.out) != bool(args.solution or args.native) or bool(args.preparation) != bool(args.native)):
        parser.error("choose --prepare DIR, --solution FILE --out FILE or --native PRIME FILE --preparation JSON --out FILE")
    if args.prepare:
        result = prepare(args.level,args.prepare)
    elif args.native:
        result = reconstruct(args.level,args.native,args.preparation,args.out)
    else:
        result = replay(args.level,args.solution,args.out)
    if args.receipt:
        args.receipt.write_text(json.dumps(result,indent=2)+"\n")
    print(result)
