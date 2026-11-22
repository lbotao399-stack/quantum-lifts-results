"""Weight bases of exterior chains of positive two-variable G2 currents."""
from collections import defaultdict
from functools import cache
from itertools import combinations


class PairComplex:
    def __init__(self, high, length, data, second):
        self.high, self.second, self.length = high, second, length
        self.weights = data["weights"]
        self.brackets = data["brackets"]
        maximum = high+second - length + 2
        self.modes = [(a, n-a) for n in range(1, maximum+1) for a in range(n+1)]

    @cache
    def exterior(self, size):
        terms = defaultdict(list)
        for inds in combinations(range(14), size):
            w = tuple(map(sum, zip(*(self.weights[i] for i in inds))))
            terms[w].append(inds)
        return terms

    @cache
    def patterns(self, length):
        if length > self.high+self.second:
            return []
        patterns = []
        def visit(start, a, b, slots, chosen):
            if a < 0 or b < 0 or slots < 0 or a+b < slots:
                return
            if a+b == slots:
                base = ([(0, b)] if b else []) + ([(1, a)] if a else [])
                patterns.append(base + chosen)
                return
            excess = a+b-slots
            for mode in range(start, len(self.modes)):
                u, v = self.modes[mode]
                if u+v-1 > excess:
                    break
                cap = min(14, slots, excess//(u+v-1), a//u if u else 14, b//v if v else 14)
                for n in range(1, cap+1):
                    visit(mode+1, a-u*n, b-v*n, slots-n, chosen+[(mode, n)])
        visit(2, self.high, self.second, length, [])
        return patterns

    def words(self, length, weight):
        for parts in self.patterns(length):
            groups = sorted(parts, key=lambda p: len(self.exterior(p[1])))
            @cache
            def reaches(k, x, y):
                if k == len(groups):
                    return x == y == 0
                if k == len(groups)-1:
                    return (x, y) in self.exterior(groups[k][1])
                return any(reaches(k+1, x-a, y-b) for a, b in self.exterior(groups[k][1]))
            def emit(k, x, y, word):
                if k == len(groups):
                    yield tuple(sorted(word))
                    return
                mode, n = groups[k]
                table = self.exterior(n)
                weights = ((x, y),) if k == len(groups)-1 else table
                for a, b in weights:
                    if reaches(k+1, x-a, y-b):
                        for cols in table.get((a, b), ()):
                            yield from emit(k+1, x-a, y-b, word+tuple(14*mode+c for c in cols))
            yield from emit(0, *weight, ())
            reaches.cache_clear()
