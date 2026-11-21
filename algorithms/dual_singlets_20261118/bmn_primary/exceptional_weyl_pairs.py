"""Two-vector Molien multiplicities from pinned CharliePy conjugacy classes."""
import argparse
import ast
from collections import Counter
import hashlib
import json
from math import prod
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEGREES = {
    "F4": (2, 6, 8, 12),
    "E6": (2, 5, 6, 8, 9, 12),
    "E7": (2, 6, 8, 10, 12, 14, 18),
    "E8": (2, 8, 12, 14, 18, 20, 24, 30),
}
RAW_SHA256 = "7793acd5b7deb1aca4636161d9c086801408da313f9a0d85e73a6421c1fa34c6"
RAW_URL = ("https://raw.githubusercontent.com/jay-taylor/charliepy/"
           "0a70d17a7f101ad5c9e8d07bfdf907d7b06083b3/charliepy/data/rawdata.py")


def denominator(group, word):
    rank = int(group[1:])
    edges = [(0,1), (1,2), (2,3)] if group == "F4" else (
        [(0,2), (1,3)] + [(i,i+1) for i in range(2,rank-1)])
    cartan = [[2*int(i == j) for j in range(rank)] for i in range(rank)]
    for i,j in edges:
        cartan[i][j] = cartan[j][i] = -1
    if group == "F4":
        cartan[1][2] = -2
    matrix = [[int(i == j) for j in range(rank)] for i in range(rank)]
    for i in word:
        # Dynkin coordinates of the image of each basis vector.
        for v in matrix:
            entry = v[i]
            for j in range(rank):
                v[j] -= entry*cartan[i][j]
    power = matrix
    traces = []
    coeff = [1]
    for k in range(1, rank+1):
        traces.append(sum(power[i][i] for i in range(rank)))
        value, rem = divmod(-sum(traces[j]*coeff[k-1-j] for j in range(k)), k)
        assert not rem
        coeff.append(value)
        power = [[sum(x*y for x,y in zip(row,col)) for col in zip(*matrix)] for row in power]
    return tuple(coeff)


def extract(path):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != RAW_SHA256:
        raise ValueError("Unexpected CharliePy data revision")
    fields = {n.targets[0].id: n.value for n in ast.parse(raw).body
              if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)}
    groups = {}
    for group, degrees in DEGREES.items():
        order = prod(degrees)
        mass = Counter()
        for _,centralizer,representative in ast.literal_eval(fields[group+"conjclasses"]):
            count, rem = divmod(order, centralizer)
            assert not rem
            mass[denominator(group, representative)] += count
        assert sum(mass.values()) == order
        groups[group] = [[list(poly), size] for poly,size in sorted(mass.items())]
    return {"source": RAW_URL, "source_sha256": RAW_SHA256, "groups": groups}


def molien(group, maximum, table):
    numerators = Counter()
    order = prod(DEGREES[group])
    assert sum(size for _,size in table[group]) == order
    for poly, size in table[group]:
        h = [1]
        for k in range(1, maximum+1):
            h.append(-sum(poly[i]*h[k-i] for i in range(1,min(k,len(poly)-1)+1)))
        for a in range(maximum+1):
            for b in range(a+1):
                numerators[a,b] += size*h[a]*h[b]
    assert all(n % order == 0 for n in numerators.values())
    return {w:n//order for w,n in numerators.items()}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("charliepy_rawdata", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.write_text(json.dumps(extract(args.charliepy_rawdata), separators=(",",":"))+"\n")
