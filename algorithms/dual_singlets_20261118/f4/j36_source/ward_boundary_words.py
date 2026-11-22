"""Differentiate the F4 test trace as jet words and evaluate its colour components."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache, lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
from zipfile import ZipFile

from verify import ARCHIVE, contents, differentiate, insert


def tree_words():
    ghosts = ((1, 0, 0), (0, 1, 0))
    scalar = tuple((0, 0, 1 << i) for i in range(3))
    source = Counter()
    for i, j, k in permutations(range(3)):
        orientation = (-1)**((i > j) + (i > k) + (j > k))
        for a, c in product(range(2), repeat=2):
            jet = (*ghosts[a][:2], 1 << i)
            word = (jet, scalar[j], ghosts[1-a], scalar[k], ghosts[c], ghosts[1-c])
            insert(source, (word,), orientation * (-1)**(a+c))
    boundary = differentiate((coefficient, word) for word, coefficient in source.items() if coefficient)
    # Expand the four exterior assignments and collect the resulting signed trace permutations.
    layouts = Counter()
    for (trace,), value in boundary.items():
        for us in ((0, 1), (1, 0)):
            for vs in ((2, 3), (3, 2)):
                u, v = iter(us), iter(vs)
                order = tuple(next(u) if slot == ghosts[0] else next(v) if slot == ghosts[1]
                              else 4+scalar.index(slot) for slot in trace)
                odd = [a for a in order if a < 4]
                inversions = sum(x > y for i, x in enumerate(odd) for y in odd[i+1:])
                layouts[order] += (-1)**inversions * value
    return {k: v for k, v in layouts.items() if v}, len(boundary)


def replay(z):
    colour = json.loads(contents(z, 'collected-O0/color-basis-native.json'))
    metric = {(i,j): Fraction(x) for i,j,x in colour['inverseMetric']}
    matrices = []
    for row in colour['nativeDualTranspose26']:
        matrix = defaultdict(list)
        for i,j,x in row:
            value = Fraction(x) * 2
            assert value.denominator == 1
            matrix[i].append((j, int(value)))
        matrices.append(matrix)

    @lru_cache(maxsize=32000)
    def matrix_word(word):
        factor = matrices[word[-1]]
        if len(word) == 1:
            return {(i,j): v for i, row in factor.items() for j,v in row}
        values = defaultdict(int)
        for (i,k), a in matrix_word(word[:-1]).items():
            for j,b in factor.get(k, ()):
                values[i,j] += a*b
        return {k: x for k,x in values.items() if x}

    layouts, word_count = tree_words()
    @cache
    def coefficient_core(colours):
        total = 0
        for layout, c in layouts.items():
            word = tuple(colours[i] for i in layout)
            first = matrix_word(word[:3])
            if first:
                last = matrix_word(word[3:])
                total += c * sum(a*last.get((j,i),0) for (i,j),a in first.items())
        return Fraction(total, 2**7)

    # Wedge-multiply three bilinears directly. The first ghost occupies the chosen
    # plus slot and the second one occupies the permuted minus slot of that factor.
    @cache
    def spectators(u, v):
        result = Fraction()
        for order in permutations(range(3)):
            value = 1
            for i,j in enumerate(order):
                value *= metric.get((u[i], v[j]), 0)
            parity = sum(x>y for i,x in enumerate(order) for y in order[i+1:])
            result -= 6*(-1)**parity*value
        return result

    def coefficient(k):
        answer = Fraction()
        for plus in combinations(range(5),3):
            for minus in combinations(range(5),3):
                value = spectators(tuple(k[i] for i in plus), tuple(k[5+j] for j in minus))
                if value:
                    u = [i for i in range(5) if i not in plus]
                    v = [i for i in range(5) if i not in minus]
                    inversions = sum(i>j for i in plus for j in u) + sum(i>j for i in minus for j in v)
                    core = tuple(k[i] for i in u) + tuple(k[5+j] for j in v) + k[10:]
                    answer += value*(-1)**inversions*coefficient_core(core)
        return answer

    keys = {tuple(k) for k,_ in json.loads(contents(z, 'cochains/inputs/exact-raw-anchors.json'))['anchors']}
    keys.update(tuple(k) for k,_ in json.loads(contents(z, 'exact-score/eta0-complete-exact-anchors.json')))
    rows = [[list(k), str(coefficient(k))] for k in sorted(keys)]
    digest = sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()
    assert len(rows) == 2981 and digest == 'bd35f8393b86358cc4f329b3ee9a9b3f51b3267b620be85a31000a0c3cc96fbd'
    return dict(canonical_tree_words=word_count, ordered_colour_layouts=len(layouts), components=len(rows), sha256=digest)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    if not __debug__ or sha256(args.archive.read_bytes()).hexdigest() != ARCHIVE:
        parser.error('use the pinned archive with Python assertions enabled')
    with ZipFile(args.archive) as archive:
        print(json.dumps(replay(archive), indent=2))
