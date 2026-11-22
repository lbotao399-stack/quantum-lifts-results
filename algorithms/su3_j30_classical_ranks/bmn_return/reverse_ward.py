"""Solve and replay the first-order reverse Taylor Ward identity in the relative complex."""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if (HERE / "fortuity").is_dir():
    sys.path.insert(0, str(HERE))
else:
    runtime = HERE.parents[1] / "su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work"
    sys.path[:0] = [str(runtime/"so8-bmn/vendor"), str(runtime/"so8-protected")]
from fortuity.full.letters import is_ferm, Q_letter
from q1d_exact import pair_terms


def add(row, key, x):
    if x:
        row[key] = row.get(key, 0)+x
        if not row[key]:
            del row[key]


@cache
def canonical(word):
    odd = [t for t in word if is_ferm(t[0])]
    sign = (-1)**sum(x > y for i, x in enumerate(odd) for y in odd[i+1:])
    return (sign if len(odd) == len(set(odd)) else 0), tuple(sorted(word))


def term(row, word, x, key=None):
    sign, target = canonical(word)
    add(row, (key, target), sign*x)


def letter_degree(l):
    return l[0]+l[1]+len(l[2])


def letters(degree):
    return [(a, b, t) for n in range(1, degree+1) for a in range(n+1)
            for b in range(n-a+1) for t in combinations(range(3), n-a-b)]


class Brackets:
    """Coordinate tensors in an invariant metric and its inverse."""
    def __init__(self, bracket, metric, inverse):
        self.n = len(bracket)
        self.indices = range(self.n)
        self.bracket = bracket
        self.metric = metric
        self.inverse = inverse
        for a, b in product(self.indices, repeat=2):
            assert sum(metric[a][j]*inverse[j][b] for j in self.indices) == (a == b)

    @cache
    def split(self, a):
        return [(b, c, row[a]) for b, rows in enumerate(self.bracket)
                for c, row in enumerate(rows) if a in row]

    @cache
    def rotate(self, a, u):
        return [(b, row[a]) for b, row in enumerate(self.bracket[u]) if a in row]

    @cache
    def join(self, a, b):
        row = {}
        for u, x in enumerate(self.inverse[a]):
            for v, y in enumerate(self.inverse[b]):
                if x and y:
                    for k, z in self.bracket[u][v].items():
                        for h, g in enumerate(self.metric[k]):
                            add(row, h, x*y*z*g)
        return list(row.items())

    @cache
    def loop(self, a, b):
        row = {}
        for u, v in product(self.indices, repeat=2):
            scale = self.inverse[u][v]
            if scale:
                for i, x in self.rotate(a, u):
                    for j, y in self.rotate(b, v):
                        add(row, (i, j), scale*x*y)
        return list(row.items())


def cross_product():
    bracket = [[{} for _ in range(3)] for _ in range(3)]
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        bracket[i][j][k] = 1
        bracket[j][i][k] = -1
    eye = [[int(i == j) for j in range(3)] for i in range(3)]
    return Brackets(bracket, eye, eye)


def traceless_matrices():
    units = [(i, j) for i in range(3) for j in range(3) if i != j]
    basis = [{ij: 1} for ij in units]+[{(0,0):1, (1,1):-1}, {(1,1):1, (2,2):-1}]
    metric = [[sum(x*b.get((j,i), 0) for (i,j), x in a.items())
               for b in basis] for a in basis]
    inverse = [row[:] for row in metric]
    inverse[6][6] = inverse[7][7] = Fraction(2, 3)
    inverse[6][7] = inverse[7][6] = Fraction(1, 3)
    def commutator(a, b):
        out = {}
        for (i,j), x in a.items():
            for (k,l), y in b.items():
                if j == k:
                    add(out, (i,l), x*y)
                if l == i:
                    add(out, (k,j), -x*y)
        coords = [out.get(m, 0) for m in units]+[out.get((0,0), 0), -out.get((2,2), 0)]
        assert sum(out.get((i,i), 0) for i in range(3)) == 0
        return {i: x for i, x in enumerate(coords) if x}
    return Brackets([[commutator(a, b) for b in basis] for a in basis], metric, inverse)


class Ward:
    def __init__(self, colour, spin, flavour):
        self.colour, self.spin, self.flavour = colour, spin, flavour

    @cache
    def single(self, token, tree):
        l, a = token
        if tree:
            row = {}
            for e, x, y in Q_letter(l):
                for b, c, f in self.colour.split(a):
                    term(row, ((x,b), (y,c)), Fraction(e*f, 2))
            return [(w, v) for (_, w), v in row.items()]
        dots = list(l[:2])
        if not dots[self.spin] or self.flavour in l[2]:
            return []
        dots[self.spin] -= 1
        z = (*dots, tuple(sorted((*l[2], self.flavour))))
        return [(((z, a),), (-1)**sum(i > self.flavour for i in l[2]))]

    def one(self, word, tree):
        out = {}
        phase = 1
        for i, t in enumerate(word):
            for new, v in self.single(t, tree):
                term(out, word[:i]+new+word[i+1:], phase*v)
            if is_ferm(t[0]):
                phase = -phase
        return [(w, v) for (_, w), v in out.items()]

    @cache
    def shifted(self, x, y):
        counts = [(i in x[2])+(i in y[2])-(i != self.flavour) for i in range(3)]
        if min(counts)<0 or max(counts)>1:
            return None
        return (x[0]+y[0]+(self.spin == 1), x[1]+y[1]+(self.spin == 0),
                tuple(i for i in range(3) if counts[i]))

    def two(self, word, correction):
        out = {}
        for i, j in combinations(range(len(word)), 2):
            (x,a), (y,b) = word[i], word[j]
            parity = [is_ferm(t[0]) for t in word]
            sign = (-1)**(parity[i]*sum(parity[:i])+parity[j]*(sum(parity[:j])-parity[i]))
            tail = word[:i]+word[i+1:j]+word[j+1:]
            if correction:
                z = self.shifted(x, y)
                swap = -(-1)**(parity[i]*parity[j])
                if z is None or (x == y and swap == -1):
                    continue
                key = ("pair", min(x,y), max(x,y))
                for c, v in self.colour.join(a, b):
                    term(out, ((z,c),)+tail, sign*v*(swap if x > y else 1), key)
            else:
                for v, s, t in pair_terms(x, y):
                    for (c,d), f in self.colour.loop(a, b):
                        term(out, ((s,c), (t,d))+tail, sign*v*f)
        return out

    def gauge(self, word):
        # Rightmost rotation is the infinitesimal simultaneous gauge action.
        out = {}
        for u in self.colour.indices:
            for i, (l,a) in enumerate(word):
                for b, x in self.colour.rotate(a, u):
                    first = word[:i]+((l,b),)+word[i+1:]
                    for j, (m,c) in enumerate(first):
                        z = self.shifted(m, (0,0,()))
                        if z is None:
                            continue
                        for v, scale in enumerate(self.colour.inverse[u]):
                            if scale:
                                for d, y in self.colour.rotate(c, v):
                                    term(out, first[:j]+((z,d),)+first[j+1:], scale*x*y, ("gauge",m))
        return out

    def equations(self, word):
        rows = defaultdict(dict)
        for (key, t), x in self.two(word, True).items():
            for z, y in self.one(t, True):
                add(rows[z], key, x*y)
        for t, x in self.one(word, True):
            for (key, z), y in self.two(t, True).items():
                add(rows[z], key, x*y)
        for (_, t), x in self.two(word, False).items():
            for z, y in self.one(t, False):
                add(rows[z], "rhs", -x*y)
        for t, x in self.one(word, False):
            for (_, z), y in self.two(t, False).items():
                add(rows[z], "rhs", -x*y)
        for (key, z), x in self.gauge(word).items():
            add(rows[z], key, x)
        return rows.values()


def examples(degree, singles, pairs):
    basis = letters(degree)
    for x in basis:
        for a in singles:
            yield ((x,a),)
    for x, y in combinations_with_replacement(basis, 2):
        if letter_degree(x)+letter_degree(y) <= degree:
            for a, b in pairs:
                w = ((x,a), (y,b))
                if canonical(w)[0]:
                    yield w


def system(ward, words, coefficients=None):
    pivots = {}
    rows = 0
    for w in words:
        for eq in ward.equations(w):
            if not eq:
                continue
            rows += 1
            rhs = eq.pop("rhs", 0)
            if coefficients is not None:
                assert sum(v*coefficients.get(k, 0) for k, v in eq.items()) == rhs, (w, eq, rhs)
                continue
            while eq:
                k = min(eq)
                c = eq[k]
                if k not in pivots:
                    pivots[k] = ({t: Fraction(v,c) for t,v in eq.items()}, Fraction(rhs,c))
                    break
                other, b = pivots[k]
                rhs -= c*b
                for t,v in other.items():
                    add(eq, t, -c*v)
            else:
                assert rhs == 0, (w,rhs)
    solution = {}
    for k in sorted(pivots, reverse=True):
        eq, rhs = pivots[k]
        value = rhs-sum(v*solution.get(t, 0) for t,v in eq.items() if t != k)
        if value:
            solution[k] = value
    return rows, len(pivots), solution


def encode(solution):
    return [[*k, str(v)] for k, v in sorted(solution.items())]


def decode(record):
    def letter(row):
        return row[0], row[1], tuple(row[2])
    return {(tag, *map(letter, items[:-1])): Fraction(items[-1]) for tag, *items in record}


STEPS = ((0,1), (0,2), (1,1), (1,2))


def solve(spin, flavour):
    ward = Ward(cross_product(), spin, flavour)
    rows, rank, values = system(ward, examples(7, (0,), ((0,0), (0,1))))
    return dict(spin=spin, flavour=flavour, degree=7, equations=rows, rank=rank, coefficients=encode(values))


def replay(record, verification, first_colour=None):
    spin, flavour = record["spin"], record["flavour"]
    values = decode(record["coefficients"])
    if verification == "cross":
        colour, singles, pairs = cross_product(), (0,), ((0,0), (0,1))
    else:
        colour = traceless_matrices()
        singles = (0,6)
        pairs = ((0,0), (0,2), (0,3), (0,5), (6,0), (6,6), (6,7))
        if verification == "full-sl3":
            singles = colour.indices if first_colour is None else (first_colour,)
            pairs = tuple(product(singles, colour.indices))
    ward = Ward(colour, spin, flavour)
    count, _, _ = system(ward, examples(7, singles, pairs), values)
    return dict(spin=spin, flavour=flavour, tested=verification, equations=count)


def colour_slice(task):
    record, colour = task
    return replay(record, "full-sl3", colour)["equations"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("solve", "cross", "sl3", "full-sl3"))
    p.add_argument("packet", type=Path)
    p.add_argument("--step", type=int, choices=range(4))
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    if not __debug__ or args.workers < 1:
        p.error("enable assertions and use a positive worker count")
    steps = STEPS if args.step is None else (STEPS[args.step],)
    if args.mode == "solve":
        result = [solve(*s) for s in steps]
        args.packet.write_text(json.dumps(result, separators=(",",":"))+"\n")
        print([(r["equations"],r["rank"],len(r["coefficients"])) for r in result])
    else:
        data = json.loads(args.packet.read_bytes())
        selected = [r for r in data if (r["spin"],r["flavour"]) in steps]
        assert [(r["spin"],r["flavour"]) for r in selected] == list(steps)
        if args.mode == "full-sl3" and args.workers > 1:
            tasks = [(r, a) for r in selected for a in range(8)]
            with ProcessPoolExecutor(args.workers) as pool:
                counts = list(pool.map(colour_slice, tasks))
            records = [dict(spin=r["spin"], flavour=r["flavour"], tested=args.mode,
                            equations=sum(counts[8*i:8*i+8]), by_input_colour=counts[8*i:8*i+8])
                       for i, r in enumerate(selected)]
        else:
            records = [replay(r, args.mode) for r in selected]
        print(json.dumps(records, indent=2))
