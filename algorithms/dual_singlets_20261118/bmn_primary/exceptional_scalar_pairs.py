"""Cartan evaluation certificates for two-scalar exceptional BMN classes."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import prod
from pathlib import Path
from random import Random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from exceptional_trace import folded_roots, root_data
from exceptional_weyl_pairs import DEGREES, molien
from f4_singlet_characters import CARTAN, reflect, roots

TARGETS = {"F4": (4,4), "E6": (4,4), "E7": (6,4), "E8": (6,6)}


def inverse_killing_form(adjoint):
    """Invert the sum of root outer products in simple coroot coordinates."""
    rank = len(adjoint[0])
    form = [[Fraction(sum(w[i]*w[j] for w in adjoint)) for j in range(rank)]
            + [Fraction(i==j) for j in range(rank)] for i in range(rank)]
    for i in range(rank):
        scale = form[i][i]
        form[i] = [v/scale for v in form[i]]
        for j,row in enumerate(form):
            if j != i:
                form[j] = [x-row[i]*y for x,y in zip(row,form[i])]
    return [row[rank:] for row in form]


@lru_cache(None)
def colour_data(group):
    if group == "F4":
        simple, adjoint = CARTAN[group], roots(group)
        short = set(folded_roots()[0])
    else:
        simple, adjoint = root_data(group)
        short = set()
    rank = len(simple)
    form = inverse_killing_form(adjoint)
    norm = {w:sum(w[i]*form[i][j]*w[j] for i,j in product(range(rank),repeat=2))
            for w in adjoint}
    assert len(set(norm.values())) == 1+bool(short)
    if short:
        assert all(norm[w] == min(norm.values()) for w in short)
    pairs = []
    total = Counter()
    for i,alpha in enumerate(adjoint):
        for j,beta in enumerate(adjoint):
            gamma = tuple(x+y for x,y in zip(alpha,beta))
            if gamma not in norm:
                continue
            step, p = beta, 0
            while True:
                step = tuple(x-y for x,y in zip(step,alpha))
                if step not in norm:
                    break
                p += 1
            factor = (p+1)**2*norm[alpha]*norm[beta]/(2*norm[gamma])
            pairs.append((i,j,factor))
            total[gamma] += factor
    # Including the two Cartan/root inputs gives the adjoint Killing Casimir.
    assert all(total[w]+2*norm[w] == 1 for w in adjoint)
    weights = {}
    if group == "E6":
        # Either minuscule endpoint of the length-two arms gives a 27.
        dominant = tuple(int(i==3) for i in range(rank))
        orbit, pending = {dominant}, [dominant]
        while pending:
            v = pending.pop()
            for s in enumerate(simple):
                w = reflect(v,s)
                if w not in orbit:
                    orbit.add(w)
                    pending.append(w)
        assert len(orbit) == 27
        weights[5] = tuple(sorted(orbit))
    return adjoint, pairs, weights


def monomials(degrees, a, b):
    letters = [(d-y,y) for d in degrees for y in range(d+1) if d-y<=a and y<=b]
    def visit(i, u, v):
        if not (u or v):
            yield ()
            return
        for j in range(i,len(letters)):
            x,y = letters[j]
            if x<=u and y<=v:
                for tail in visit(j,u-x,v-y):
                    yield ((x,y),)+tail
    return list(visit(0,a,b))


def evaluation(group, point, prime, maximum, target):
    x,y = point
    adjoint, pairs, other = colour_data(group)
    def restrict(weights):
        return [(sum(i*j for i,j in zip(w,x)),sum(i*j for i,j in zip(w,y))) for w in weights]
    root_values = restrict(adjoint)
    values = {}
    for d in DEGREES[group]:
        if d>maximum:
            break
        ev = restrict(other[d]) if d in other else root_values
        for k in range(d+1):
            trace = sum(u**(d-k)*v**k for u,v in ev)
            values[d-k,k] = trace if prime is None else trace%prime
    a,b = target
    coupled = 0
    for i,j,factor in pairs:
        u,v = root_values[i]
        s,t = root_values[j]
        coefficient = factor if prime is None else factor.numerator*pow(factor.denominator,-1,prime)%prime
        coupled += coefficient*(u*t-v*s)**b*u**(a-b)
    return values, coupled if prime is None else coupled%prime


def exact_determinant(rows):
    matrix = [list(map(Fraction,r)) for r in rows]
    det = 1
    for i,row in enumerate(matrix):
        # These rows follow the modular pivot order, so every pivot is nonzero.
        pivot = row[i]
        if not pivot:
            raise ValueError("Rational pivot vanished")
        det *= pivot
        for j in range(i+1,len(matrix)):
            ratio = matrix[j][i]/pivot
            matrix[j] = [u-ratio*v for u,v in zip(matrix[j],row)]
    return det


def minor(columns, prime):
    if not columns:
        return {"rank":0, "rows":[], "columns":[], "determinant":1}
    matrix = [list(r) for r in zip(*columns)]
    rows = list(range(len(matrix)))
    pivots = []
    determinant = 1
    for col in range(len(columns)):
        n = len(pivots)
        pivot = next((i for i in range(n,len(rows)) if matrix[i][col]), None)
        if pivot is None:
            continue
        matrix[pivot],matrix[n] = matrix[n],matrix[pivot]
        rows[pivot],rows[n] = rows[n],rows[pivot]
        top = matrix[n]
        determinant = determinant*top[col]%prime
        inverse = pow(top[col],-1,prime)
        for i in range(n+1,len(rows)):
            row = matrix[i]
            a = row[col]*inverse%prime
            if a:
                for j in range(col,len(columns)):
                    row[j] = (row[j]-a*top[j])%prime
        pivots.append(col)
    return {"rank":len(pivots), "rows":rows[:len(pivots)], "columns":pivots, "determinant":determinant}


def calculate(group, prime, seed, table):
    a,b = TARGETS[group]
    degrees = DEGREES[group]
    dims = molien(group, a+b, table)
    # All possible factors of positive bidegree in a polynomial of weight (a,b).
    weights = sorted({tuple(sorted((u,v),reverse=True)) for u in range(a+1) for v in range(b+1)})
    words = {w:monomials(degrees,*w) for w in weights}
    rng = Random(seed)
    npoints = max(map(len,words.values()))+6
    points = [[[rng.randrange(-12,13) for _ in range(int(group[1:]))] for _ in range(2)]
              for _ in range(npoints)]
    evaluated = [evaluation(group,p,prime,a+b,(a,b)) for p in points]
    records = []
    for w in weights:
        columns = [[prod(values[key] for key in term)%prime for values,_ in evaluated] for term in words[w]]
        basis = minor(columns,prime)
        row = {"bidegree":w,"weyl_dimension":dims[w],"number_of_products":len(columns),"minor":basis}
        if w == (a,b):
            row["augmented_minor"] = minor(columns+[[v for _,v in evaluated]], prime)
            row["products"] = words[w]
            row["weyl_highest_weight_multiplicity"] = dims[w]-dims[a+1,b-1]
            cert = row["augmented_minor"]
            assert basis["rank"] == len(columns) and cert["rank"] == len(columns)+1 == dims[w]
            exact = []
            for i in cert["rows"]:
                ev, coupled = evaluation(group,points[i],None,a+b,(a,b))
                exact.append([prod(ev[key] for key in term) for term in words[w]]+[coupled])
            det = exact_determinant(exact)
            assert det.numerator*pow(det.denominator,-1,prime)%prime == cert["determinant"]
            row["rational_augmented_determinant"] = str(det)
        else:
            assert basis["rank"] == dims[w], (group,w,basis["rank"],dims[w])
        records.append(row)
    return {"group":group,"flavour_partition":[a,b,0],"level":2*(a+b),"prime":prime,
            "seed":seed,"cartan_points":points,"records":records}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prime",type=int,choices=(1000003,1000033),default=1000003)
    p.add_argument("--seed",type=int,default=70319)
    p.add_argument("--out",type=Path,required=True)
    args = p.parse_args()
    table = json.loads((HERE/"exceptional_weyl_classes.json").read_bytes())["groups"]
    runs = [calculate(group,args.prime,args.seed,table) for group in TARGETS]
    for s in runs:
        last = next(r for r in s["records"] if "augmented_minor" in r)
        print(s["group"],last["weyl_dimension"],last["number_of_products"],last["minor"]["rank"],
              last["augmented_minor"]["rank"],last["augmented_minor"]["determinant"],
              "lower exceptions", [r["bidegree"] for r in s["records"] if r is not last
                                    and r["weyl_dimension"]!=r["minor"]["rank"]])
    args.out.write_text(json.dumps(runs,separators=(",",":"))+"\n")
