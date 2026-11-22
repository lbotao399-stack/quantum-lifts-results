"""Bare fifth-glueball colour factor on the doublets of a long-root F4 instanton."""
from collections import Counter
from itertools import combinations
import json
from math import factorial, prod

from colour import Colour, Q, add, compact, wedge


def calculate():
    c = Colour()
    e, = [i for i in c.positive if all(not c.br[j][i] for j in c.data["raising"])]
    f = c.opposite[e]
    h = c.br[e][f]
    factor = Q(2, sum(v*c.weights[e][i] for i,v in h.items()))
    weights = [sum(factor*v*c.weights[j][i] for i,v in h.items()) for j in range(52)]
    counts = sorted((int(w), k) for w,k in Counter(weights).items())
    assert counts == [(-2,1),(-1,14),(0,22),(1,14),(2,1)]
    high = [i for i,w in enumerate(weights) if w == 1]
    low = [{a:factor*v for a,v in c.br[f][i].items()} for i in high]
    form = {}
    for i,x in enumerate(high):
        raised = {}
        for a,v in low[i].items():
            for b,y in c.br[e][a].items():
                add(raised,b,v*y)
        assert compact(raised) == {x:1}
        for j,y in enumerate(high):
            n = sum(c.metric.get((x,a),0)*v for a,v in low[j].items())
            opposite = sum(c.metric.get((y,a),0)*v for a,v in low[i].items())
            assert n == -opposite and int(n) == n
            if i < j and n:
                form[x,y] = int(n)
    assert len(form) == 7 and sorted(k for pair in form for k in pair) == high
    row = {():1}
    for _ in range(5):
        out = {}
        for monomial,v in row.items():
            for pair,n in form.items():
                terms = (*monomial,*pair)
                _,sg = wedge(terms)
                add(out,tuple(sorted(terms)),sg*2*v*n)
        row = compact(out)
    check = {}
    for selection in combinations(form,5):
        terms = tuple(x for pair in selection for x in pair)
        _,sg = wedge(terms)
        check[tuple(sorted(terms))] = sg*factorial(5)*2**5*prod(form[k] for k in selection)
    assert row == check and len(row) == 21
    return dict(raising_colour=e, opposite_colour=f, lowering_coefficient=int(factor),
                weights=counts, doublets=high, alternating_form=[[*p,v] for p,v in form.items()],
                fifth_power=[[list(m),v] for m,v in sorted(row.items())],
                native_metric="Tr26/6", physical_instanton_supercharge_pairing="open")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
