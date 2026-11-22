"""Rank over Q using sparse integer row elimination with primitive row normalization."""
from collections import defaultdict
from heapq import heapify, heappop, heappush
from math import gcd

def rational_rank(columns):
    rows = defaultdict(dict)
    incidence = defaultdict(set)
    for col, entries in enumerate(columns):
        for i, v in entries.items():
            if v:
                rows[i][col] = v
                incidence[col].add(i)
    queue = [(len(row), i) for i, row in rows.items()]
    heapify(queue)
    pivots = 0
    while queue:
        size, i = heappop(queue)
        if i not in rows or size != len(rows[i]):
            continue
        pivot = rows.pop(i)
        j = min(pivot, key=lambda k: (len(incidence[k]), abs(pivot[k]), k))
        for k in pivot:
            incidence[k].remove(i)
        a = pivot[j]
        for other in tuple(incidence[j]):
            row = rows[other]
            b = row[j]
            divisor = gcd(a, b)
            m, n = a // divisor, b // divisor
            if m != 1:
                for k in row:
                    row[k] *= m
            for k, c in pivot.items():
                v = row.get(k, 0) - n*c
                if not v:
                    del row[k]
                    incidence[k].remove(other)
                else:
                    if k not in row:
                        incidence[k].add(other)
                    row[k] = v
            if row:
                if m != 1 or abs(n) != 1:
                    divisor = gcd(*row.values())
                    if divisor > 1:
                        for k in row:
                            row[k] //= divisor
                heappush(queue, (len(row), other))
            else:
                del rows[other]
        pivots += 1
    return pivots
