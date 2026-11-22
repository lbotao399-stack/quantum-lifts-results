"""Rebuild the rational first images of scalar-dressed SO(11) and SO(13) glueballs."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import sys

from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from mixed_incoming import so_canon
from orthogonal_glueball_first_pair import GHOSTS, PHI, PfaffianCoefficient, orientations
from orthogonal_glueball_relation_span import determinant
from orthogonal_glueball_third_pair import encode
from orthogonal_trace_dual import pullback
from rank_four_lowered_replay import terminal
from terminal_commutators import words

A = PHI[0]
ALPHABET = PHI + GHOSTS
RECIPES = [
    (dict(paths=['mp','mp','ps'], single=[[0,1],[1,2,0]],
          matrices=[[1,-2,-2,-1],[2,2,2,-1],[0,2,0,2]],
          scalar_directions=[[],[],[[-1,1,1]]]), (1,2)),
    (dict(paths=['mp','mp'], single=[[0,1,2],[1,2,0]],
          matrices=[[-1,-1,-1,1],[-2,0,-2,-2]], scalar_directions=[[],[]]), (1,2)),
    (dict(paths=['mp','mps'], single=[[2,1,0],[1,2,0]],
          matrices=[[2,1,2,1],[-1,-2,0,0]], scalar_directions=[[],[[-2,2,-1]]]), (0,2)),
]
# Source row, letter, and insertion geometry for the six quartic columns.
INSERTIONS = [(0,0,'middle'),(0,1,'ends'),(0,3,'ends'),(0,4,'ends'),
              (1,0,'middle'),(2,0,'middle')]


def unit(size, index):
    return [int(j == index) for j in range(size)]


def scalar_projector(order, axes):
    """Differentiate the scalar-marker determinant and then along the extra flavour-one axes."""
    project = PfaffianCoefficient(order, 0)
    result = {}
    divisors = set()
    for old, coefficient in project.weights.items():
        counts = [(old >> (3*i)) & 7 for i in range(9)]
        for axis in axes:
            counts[axis] += 1
            assert counts[axis] < 8
            coefficient *= counts[axis]
        weight = sum(n << (3*i) for i,n in enumerate(counts))
        result[weight] = coefficient
        partial = {0}
        for i,n in enumerate(counts):
            for _ in range(n):
                partial |= {v + (1 << (3*i)) for v in partial}
        divisors.update(partial)
    project.weights, project.divisors = result, divisors
    return project


def pfaffian(spec, axes, order, length, colours, closed=False):
    # determinant supplies the ordered input catalogue. Set the actual degree
    # after replacing its default six-scalar projector by the required one.
    rel = determinant(6, spec)
    rel.project = scalar_projector(order, axes)
    rel.closed = closed
    rel.letters = length
    rel.pf_degree = length - sum(len(v)-1 for v in spec['paths']) - 2*int(closed)
    assert rel.pf_degree > colours
    return rel


def target_coordinates(weight):
    basis = words(weight)
    columns = [dict() for _ in basis]
    for i,monomial in enumerate(basis):
        choices = []
        for trace in monomial:
            letters = tuple(x for x,n in zip(GHOSTS+PHI,trace) for _ in range(n))
            choices.append(orientations(letters))
        for ordered in product(*choices):
            sign, word = so_canon(ordered)
            if sign:
                coordinate, check = terminal(word)
                assert check == monomial and coordinate
                columns[i][word] = coordinate
    return basis, columns


def first_duals(colours, weight, indices):
    basis, coords = target_coordinates(weight)
    rows, checks = [], []
    for j in indices:
        row, info = pullback(colours, coords[j], 1)
        rows.append(row)
        checks.append(info['checked_columns'])
    return rows, checks


def evaluate(row, coefficient):
    return sum(c*coefficient(w) for w,c in row.items())


def first_matrix(colours, recipes, power):
    forms, checks = first_duals(colours, (power,power,2,0,0), range(3))
    matrix = []
    for spec,axes in recipes:
        rel = pfaffian(spec, axes, 1, 2*power+3, colours)
        matrix.append([str(evaluate(row, lambda w: rel.coefficient(encode(w)))) for row in forms])
        rel.clear()
    return matrix,checks


def replacements(letter, geometry):
    if geometry == 'middle':
        return ((A,letter,A),)
    assert geometry == 'ends'
    return ((A,A,letter),(letter,A,A))


@lru_cache(8192)
def insert(word, letter, geometry):
    out = Counter()
    for i,trace in enumerate(word):
        for j,x in enumerate(trace):
            if x == letter:
                for seq in replacements(x,geometry):
                    sign, w = so_canon(word[:i]+(trace[:j]+seq+trace[j+1:],)+word[i+1:])
                    if sign:
                        out[w] += sign
    return out


@lru_cache(8192)
def insertion_parents(word, letter, geometry):
    candidates = set()
    for i,trace in enumerate(word):
        if len(trace) < 4:
            continue
        for j in range(len(trace)):
            rotation = trace[j:]+trace[:j]
            if rotation[:3] in replacements(letter,geometry):
                sign, w = so_canon(word[:i]+((letter,)+rotation[3:],)+word[i+1:])
                if sign:
                    candidates.add(w)
    return [(w,insert(w,letter,geometry)[word]) for w in candidates if insert(w,letter,geometry)[word]]


def quartic_first():
    weight = (6,6,4,0,0)
    basis = words(weight)
    assert len(basis) == 9
    # Multiplication of the quadratic three rows by Tr(A^2) fills these slots.
    square = (0,0,2,0,0)
    lower = words((6,6,2,0,0))
    assert list(basis[:3]) == [(square,)+v for v in lower]
    forms, checks = first_duals(11,weight,range(3,9))
    matrix = []
    for index,letter,geometry in INSERTIONS:
        spec,axes = RECIPES[index]
        rel = pfaffian(spec,axes,1,15,11)
        @lru_cache(None)
        def coefficient(word):
            return sum(n*rel.coefficient(encode(v)) for v,n in insertion_parents(word,ALPHABET[letter],geometry))
        matrix.append([str(evaluate(row,coefficient)) for row in forms])
        coefficient.cache_clear()
        rel.clear()
    return matrix,checks


def path_types(plus, minus, scalars, budget):
    """Enumerate all internally sorted nonsingle open paths and even closed cycles."""
    types = []
    for p,m,s in product(range(2),range(2),range(scalars+1)):
        word = 'm'*m+'p'*p+'s'*s
        if len(word) < 2:
            continue
        for closed in (False,True):
            cost = len(word)-1+int(closed)
            if cost <= budget and (not closed or len(word)%2 == 0):
                types.append((closed,word,cost,p,m,s))
    types.sort()
    def visit(start,remaining,p,m,s,paths):
        if min(p,m,s) < 0 or max(0,p-3)+max(0,m-3) > 2*remaining:
            return
        if max(p,m) <= 3:
            yield paths,p,m,s
        for i in range(start,len(types)):
            closed,word,cost,np,nm,ns = types[i]
            if cost <= remaining:
                yield from visit(i,remaining-cost,p-np,m-nm,s-ns,paths+((closed,word),))
    return tuple(visit(0,budget,plus,minus,scalars,()))


def quartic_tree():
    basis,coords = target_coordinates((6,6,4,0,0))
    selected,checked = [],0
    shapes = path_types(6,6,4,4)
    assert len(shapes) == 8
    for paths,plus,minus,_ in shapes:
        closed = [v for c,v in paths if c]
        assert closed in ([],['mp'])
        opens = [v for c,v in paths if not c]
        marked = sum(v.count('s') for v in opens)
        for matrices in product(*(range(3) if len(v)%2 else [3] for v in opens)):
            for ps,ms in product(combinations(range(3),plus),combinations(range(3),minus)):
                spec = dict(paths=opens, single=[ps,ms], matrices=[unit(4,i) for i in matrices],
                            scalar_directions=[[unit(3,0)]*v.count('s') for v in opens])
                for axes in combinations_with_replacement(range(3),4-marked):
                    if sum((0,1,-1,0)[i] for i in ps+ms+axes+matrices):
                        continue
                    rel = pfaffian(spec,(0,)*marked+axes,0,16,11,closed=bool(closed))
                    row = [str(rel.commuting_coefficient(w)) for w in basis]
                    ordered = [str(evaluate(col,lambda w: rel.coefficient(encode(w)))) for col in coords]
                    assert row == ordered
                    rel.clear()
                    checked += 1
                    if fmpq_mat(selected+[row]).rank() > len(selected):
                        selected.append(row)
    assert checked == 49 and len(selected) == 3
    return checked,selected


def run():
    assert not path_types(6,6,2,2) and not path_types(7,6,1,2)
    assert not path_types(7,7,2,2) and not path_types(7,7,0,0)
    quadratic, lower_checks = first_matrix(11,RECIPES,6)
    assert str(fmpq_mat(quadratic).det()) == '-144000'
    quartic, higher_checks = quartic_first()
    assert str(fmpq_mat(quartic).det()) == '843448320000000'
    count, identities = quartic_tree()
    assert path_types(6,6,5,3) == ((((False,'mp'),)*3,3,3,5),)
    spec = dict(paths=['mp']*3,single=[[0,1,2],[0,1,2]],matrices=[unit(4,3)]*3,scalar_directions=[[]]*3)
    thirteen, thirteen_checks = first_matrix(13,[(spec,(0,0)),(spec,(1,2))],7)
    assert thirteen == [['-60','120','1080'],['-30','60','540']]
    return dict(so11_quadratic=quadratic, quadratic_checked_columns=lower_checks,
                so11_quartic_six_columns=quartic, quartic_checked_columns=higher_checks,
                quartic_tree_polarizations=count, quartic_tree_basis=identities,
                so13_quadratic=thirteen, so13_checked_columns=thirteen_checks)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with assertions enabled')
    print(json.dumps(run(),indent=2))
