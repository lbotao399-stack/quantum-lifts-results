"""Replay the first classical bracket joining the two pairs of Taylor axes over Q."""
from collections import Counter
from fractions import Fraction
from functools import cache
import gzip
import json
from math import factorial, gcd, lcm
from pathlib import Path

from currents import PairComplex

HERE = Path(__file__).resolve().parent
COLOUR = json.loads((HERE / "colour.json").read_bytes())
WEIGHTS = list(map(tuple, COLOUR["weights"]))
BRACKETS = COLOUR["brackets"]
SIMPLE = ((2, -1), (-3, 2))
RAISING = [WEIGHTS.index(a) for a in SIMPLE]
LOWERING = [WEIGHTS.index(tuple(-n for n in a)) for a in SIMPLE]


def orient(word):
    if len(set(word)) != len(word):
        return (), 0
    sign = (-1) ** sum(a > b for i, a in enumerate(word) for b in word[i+1:])
    return tuple(sorted(word)), sign


def combine(pieces):
    result = Counter()
    for coefficient, vector in pieces:
        for w, n in vector.items():
            result[w] += coefficient*n
    return {w:n for w,n in result.items() if n}


def colour_action(vector, generator, dual=False):
    table = BRACKETS[generator]
    if dual:
        reverse = [[] for _ in WEIGHTS]
        for c, column in enumerate(table):
            for dest, n in column:
                reverse[dest].append((c, -n))
        table = reverse
    result = Counter()
    for word,n in vector.items():
        for i,(d,c,a) in enumerate(word):
            for dest,m in table[c]:
                changed, sign = orient(word[:i]+((d,dest,a),)+word[i+1:])
                if sign:
                    result[changed] += sign*m*n
    return {w:n for w,n in result.items() if n}


def taylor_action(vector, raise_weight):
    result = Counter()
    for word,n in vector.items():
        for i,(d,c,a) in enumerate(word):
            m = d-a if raise_weight else a
            if m:
                dest = (d,c,a+(1 if raise_weight else -1))
                updated, sign = orient(word[:i]+(dest,)+word[i+1:])
                if sign:
                    result[updated] += sign*m*n
    return {w:n for w,n in result.items() if n}


def determinant_cycle(source):
    initial = combine((n,{key:sign}) for w,n in source.items() for key,sign in [orient(w)])
    pieces = [(1,initial)]
    raised = initial
    for k in range(1,6):
        raised = taylor_action(raised,True)
        lowered = raised
        for _ in range(k):
            lowered = taylor_action(lowered,False)
        pieces.append((Fraction((-1)**k, factorial(k)*factorial(k+1)),lowered))
    result = combine(pieces)
    common = lcm(*(Fraction(n).denominator for n in result.values()))
    integral = {w:int(n*common) for w,n in result.items()}
    divisor = gcd(*integral.values())
    return {w:n//divisor for w,n in integral.items()}


def boundary(vector):
    result = Counter()
    for word,n in vector.items():
        for i,(d,c,a) in enumerate(word):
            for j in range(i+1,len(word)):
                e,b,t = word[j]
                rest = word[:i]+word[i+1:j]+word[j+1:]
                for g,m in BRACKETS[c][b]:
                    changed,sign = orient(((d+e,g,a+t),)+rest)
                    if sign:
                        result[changed] += n*m*sign*(-1)**(i+j-1)
    return {w:n for w,n in result.items() if n}


def cap(vector, colour):
    return combine((n*(-1)**i, {w[:i]+w[i+1:]:1})
                   for w,n in vector.items() for i,letter in enumerate(w) if letter == (1,colour,0))


def value(dual, vector):
    return sum(n*dual.get(w,0) for w,n in vector.items())


def solve(rows, rhs):
    """Use the first independent columns and set all free coefficients to zero."""
    matrix = [list(map(Fraction,r+[k])) for r,k in zip(rows,rhs)]
    pivots = []
    for j in range(len(rows[0])):
        top = len(pivots)
        pivot = next((i for i in range(top,len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            continue
        matrix[top],matrix[pivot] = matrix[pivot],matrix[top]
        divisor = matrix[top][j]
        matrix[top] = [n/divisor for n in matrix[top]]
        for i,r in enumerate(matrix):
            if i != top and r[j]:
                factor = r[j]
                matrix[i] = [n-factor*m for n,m in zip(r,matrix[top])]
        pivots.append(j)
    assert all(not any(r) for r in matrix[len(pivots):])
    result = [Fraction(0)]*len(rows[0])
    for i,j in enumerate(pivots):
        result[j] = matrix[i][-1]
    assert all(sum(n*m for n,m in zip(row,result)) == b for row,b in zip(rows,rhs))
    return result


@cache
def paths(a,b):
    if min(a,b) < 0:
        return ()
    if a == b == 0:
        return ((),)
    return tuple(s+(i,) for i in (0,1) for s in paths(a-(i==0),b-(i==1)))


def mixed_component(highest, lowest):
    removed = WEIGHTS.index((-1,0))
    mixed = WEIGHTS.index((0,-1))
    remainder = tuple(letter for letter in lowest if letter != (1,removed,0))

    @cache
    def low_vector(path,dual):
        if not path:
            return {lowest:1}
        generator = (LOWERING if dual else RAISING)[path[-1]]
        return colour_action(low_vector(path[:-1],dual),generator,dual)

    @cache
    def high_vector(path):
        if not path:
            return highest
        return colour_action(high_vector(path[:-1]),LOWERING[path[-1]])

    terms, pieces = [], []
    for u in range(14):
        for v in range(14):
            bracket = dict(BRACKETS[u][v]).get(mixed,0)
            if not bracket:
                continue
            word,sign = orient(((1,v,0),)+remainder)
            if not sign:
                continue
            dx,dy = (n-m for n,m in zip(WEIGHTS[v],WEIGHTS[removed]))
            strings = paths(2*dx+3*dy,dx+2*dy)
            if not strings:
                continue
            frame = [low_vector(s,False) for s in strings]
            pairing = [[value(low_vector(s,True),y) for s in strings] for y in frame]
            coefficients = solve(pairing,[sign*y.get(word,0) for y in frame])
            for s,c in zip(strings,coefficients):
                if c:
                    terms.append([u,v,list(s),str(bracket*c)])
                    pieces.append((bracket*c,cap(high_vector(s),u)))
    return combine(pieces),terms


def vector(rows):
    return {tuple(map(tuple,w)):n for w,n in rows}


def graded(word):
    return (len(word),sum(a for d,c,a in word),sum(d-a for d,c,a in word),
            sum(WEIGHTS[c][0] for d,c,a in word),sum(WEIGHTS[c][1] for d,c,a in word))


def verify():
    raw = json.loads(gzip.decompress((HERE/"cross_witness.json.gz").read_bytes()))
    source = vector(raw["determinant_cycle"])
    highest = determinant_cycle(source)
    assert len(highest)==6573 and all(graded(w)==(8,5,5,8,0) for w in highest)
    assert not boundary(highest)
    assert all(not colour_action(highest,g) for g in RAISING)
    assert all(not taylor_action(highest,r) for r in (True,False))

    colours = [WEIGHTS.index(w) for w in ((-3,1),(-1,0),(0,-1))]
    lowest = tuple(sorted((1,c,a) for c in colours for a in (0,1)))
    assert graded(lowest)==(6,3,3,-8,0) and not boundary({lowest:1})
    assert all(not colour_action({lowest:1},g) for g in LOWERING)
    assert all(not taylor_action({lowest:1},r) for r in (True,False))

    residual,terms = mixed_component(highest,lowest)
    assert all(graded(w)==(7,5,4,7,1) for w in residual) and not boundary(residual)
    dual = vector(raw["residual_dual"])
    assert all(graded(w)==(7,5,4,7,1) for w in dual)
    pair = PairComplex(5,7,COLOUR,4)
    letters = [(a+b,c,a) for a,b in pair.modes for c in range(14)]
    dims = [sum(1 for _ in pair.words(k,(7,1))) for k in (6,7,8)]
    for indices in pair.words(8,(7,1)):
        predecessor = tuple(letters[i] for i in indices)
        assert value(dual,boundary({predecessor:1})) == 0
    pairing = value(dual,residual)
    assert dims==[83,314,230] and pairing==-405141750
    return dict(highest_cycle_terms=len(highest),mixed_colour=[0,-1],
                low_removed_colour=[-1,0],terms=terms,residual_terms=len(residual),
                residual_dimensions=dims,dual_terms=len(dual),predecessors_tested=dims[-1],
                exact_pairing=str(pairing),physical_classical_target_upper=0)


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(),indent=2))
