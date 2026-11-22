"""Test the six-scalar weight-zero remainder on the three invariant source cochains."""
from itertools import combinations, permutations
import json
from random import Random

from colour import Colour, add, compact
from source_symmetry import multiplicity, singlet_trace


def sign(t):
    return (-1)**sum(x > y for i,x in enumerate(t) for y in t[i+1:])


def probes(number):
    rng = Random(112336)
    for _ in range(number):
        yield [{a:rng.choice((-2,-1,1,2)) for a in rng.sample(range(52),9)} for _ in range(3)]


def row(c, scalars, recipes):
    def bracket(x, y):
        out = {}
        for i,a in x.items():
            for j,b in y.items():
                for k,z in c.br[i][j].items():
                    add(out,k,a*b*z)
        return compact(out)

    def metric(x, y):
        return sum(a*b*c.metric.get((i,j),0) for i,a in x.items() for j,b in y.items())

    out = [0]*len(recipes)
    for i in range(3):
        split = bracket(scalars[(i+1)%3],scalars[(i+2)%3])
        for p in permutations(range(3)):
            fields = [split,scalars[i],*(scalars[j] for j in p)]
            pairs = {(a,b):bracket(fields[a],fields[b]) for a,b in combinations(range(5),2)}
            for k,(which,a,b,d,e,f) in enumerate(recipes):
                if which == "metric":
                    value = metric(fields[a],fields[b])*metric(fields[d],pairs[e,f])
                else:
                    value = metric(pairs[a,b],bracket(fields[d],pairs[e,f]))
                out[k] += sign(p)*value
    return out


def run():
    c = Colour()
    recipes = []
    for a,b in combinations(range(5),2):
        remaining = sorted(set(range(5))-{a,b})
        recipes.append(("metric",a,b,*remaining))
        for d in remaining:
            e,f = (v for v in remaining if v != d)
            recipes.append(("tree",a,b,d,e,f))
    inputs = list(probes(3))
    data = [row(c,s,recipes) for s in inputs]
    for cols in combinations(range(len(recipes)),3):
        minor = [[r[k] for k in cols] for r in data]
        determinant = sum(sign(p)*minor[0][p[0]]*minor[1][p[1]]*minor[2][p[2]] for p in permutations(range(3)))
        if determinant:
            break
    else:
        raise ArithmeticError("The three projected images are dependent on these probes")
    dimensions = [multiplicity((3,3),k,singlet_trace) for k in (4,5,6)]
    assert dimensions == [0,3,4]
    return dict(invariant_dimensions=dimensions, scalar_flavour=[2,2,2],
                recipes=[recipes[k] for k in cols], probes=inputs, minor=minor, determinant=determinant)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(run(),indent=2))
