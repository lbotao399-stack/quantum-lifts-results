"""Reconstruct the j=36 native tree preimage by exact trace-block elimination."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import gzip
import json
import sys
from zipfile import ZipFile

from verify import ARCHIVE, HERE, contents, source

sys.path.insert(0,str(HERE.parents[1]))
from mixed_incoming import differential, q, so_canon


@cache
def trace_weight(trace):
    return (sum(j[0] for j in trace),sum(j[1] for j in trace),
            *(sum(i in j[2] for j in trace) for i in range(3)))


@cache
def sequences(weight, length):
    if not length:
        return ((),) if not any(weight) else ()
    result = []
    slots = [range(x+1) for x in weight[:2]]+[range(min(x,1)+1) for x in weight[2:]]
    for first in product(*slots):
        if not 1 <= sum(first) <= sum(weight)-length+1:
            continue
        jet = (first[0],first[1],tuple(i for i in range(3) if first[i+2]))
        remaining = tuple(x-y for x,y in zip(weight,first))
        result.extend((jet,)+tail for tail in sequences(remaining,length-1))
    return tuple(result)


@cache
def one_trace(weight, excess):
    length = sum(weight)-excess
    if length < 2:
        return ()
    words = set()
    for word in sequences(weight,length):
        sign,canonical = so_canon((word,))
        if sign:
            words.add(canonical[0])
    return tuple(sorted(words))


def incoming(profile):
    basis = set()
    for i in range(len(profile)):
        for j in range(i,len(profile)):
            factors = [one_trace(weight,(k==i)+(k==j)) for k,weight in enumerate(profile)]
            for traces in product(*factors):
                sign,word = so_canon(traces)
                if sign:
                    basis.add(word)
    return sorted(basis)


def add(row, other, scale):
    for key,value in other.items():
        v = row.get(key,0) + scale*value
        if v:
            row[key] = v
        else:
            row.pop(key,None)


def solve(profile, target):
    pivots = {}
    basis = incoming(profile)
    for i,word in enumerate(basis):
        row = dict(differential(word))
        pre = {i:1}
        while row:
            key = min(row)
            if key not in pivots:
                scale = Fraction(1,row[key])
                pivots[key] = ({k:v*scale for k,v in row.items()},{k:v*scale for k,v in pre.items()})
                break
            left,right = pivots[key]
            scale = -row[key]
            add(row,left,scale)
            add(pre,right,scale)
    residual = dict(target)
    solution = {}
    while residual:
        key = min(residual)
        left,right = pivots[key]
        value = residual[key]
        add(residual,left,-value)
        add(solution,right,value)
    return {basis[i]:c for i,c in solution.items()}


def construct(archive):
    assert sha256(archive.read_bytes()).hexdigest() == ARCHIVE
    with ZipFile(archive) as z:
        raw = map(json.loads,contents(z,'collected-O0/O0-canonical.jsonl').splitlines())
        original = source(raw)
    blocks = defaultdict(dict)
    for traces,c in original.items():
        decoded = tuple(tuple((p,q,tuple(i for i in range(3) if mask & (1<<i))) for p,q,mask in t)
                        for t in traces)
        sign,word = so_canon(decoded)
        assert sign
        profile = tuple(sorted(trace_weight(t) for t in word))
        assert word not in blocks[profile]
        blocks[profile][word] = -sign*c
    pre = {}
    for profile,target in sorted(blocks.items()):
        local = solve(profile,target)
        assert q(local) == target
        pre.update(local)
    encoded = [[str(c),[[[p,q,sum(1<<i for i in s)] for p,q,s in t] for t in w]]
               for w,c in sorted(pre.items())]
    expected = json.loads(gzip.decompress((HERE/'tree_preimage.json.gz').read_bytes()))
    assert encoded == expected
    print(len(blocks),'blocks,',len(pre),'native preimage terms, all stored coefficients reproduced')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('assertions must be enabled')
    construct(args.archive)
