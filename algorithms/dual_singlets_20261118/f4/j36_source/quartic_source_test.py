"""Test the physical tree map on the 20 F4 source recipes and four dressings."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
from zipfile import ZipFile

from dressed_tree import add_scaled, dressings, nonzero, put, verify
from stable_trace_ranks import check
from verify import ARCHIVE, HERE, contents, native_splits, parity


def raised(polynomial, old, new):
    result = Counter()
    for traces, value in polynomial.items():
        for i, word in enumerate(traces):
            for j, letter in enumerate(word):
                if letter == old:
                    changed = word[:j] + (new,) + word[j+1:]
                    put(result, traces[:i] + (changed,) + traces[i+1:], value)
    return nonzero(result)


def determinant_weight_two(polynomial):
    assert polynomial
    for word in polynomial:
        letters = Counter(x for trace in word for x in trace)
        assert [letters[i] for i in range(3)] == [2,2,2]
    for i, j in permutations(range(3), 2):
        assert not raised(polynomial, i, j)


def input_check(archive):
    bank = json.loads(contents(archive, 'cochains/project-bank.json'))['quotient_seeds']
    pairs = {5:(0,1), 6:(0,2), 7:(1,2)}
    longest = 0
    for index, seed in enumerate(bank):
        split, expected, exchanged = Counter(), Counter(), Counter()
        for traces, coeff in seed['out']:
            word = tuple(map(tuple, traces))
            put(expected, word, int(coeff))
            swapped = tuple(tuple(7-x if x in (3,4) else x for x in t) for t in word)
            put(exchanged, swapped, int(coeff))
            assert sum(x == 3 for t in word for x in t) == sum(x == 4 for t in word for x in t) == 1
            longest = max(longest, *map(len, word))
        assert nonzero(exchanged) == nonzero(expected)
        determinant_weight_two(nonzero(expected))
        for traces, coeff in seed['source']:
            word = tuple(map(tuple, traces))
            (i, j, marker), = [(i,j,x) for i,t in enumerate(word) for j,x in enumerate(t)
                               if x in pairs]
            # The last seed is a dotted derivative of scalar marker-5, with
            # the opposite dotted gaugino. Its output includes both assignments.
            pair = pairs[marker] if index < 19 else (3, marker-5)
            for (x,y), sign in ((pair,1),(pair[::-1],-1)):
                changed = word[i][:j] + (x,y) + word[i][j+1:]
                out = word[:i] + (changed,) + word[i+1:]
                put(split, out, sign*int(coeff))
                if index == 19:
                    swapped = tuple(tuple(7-z if z in (3,4) else z for z in t) for t in out)
                    put(split, swapped, sign*int(coeff))
        assert nonzero(split) == nonzero(expected), index
    factors = dressings()
    for factor in factors:
        determinant_weight_two(factor)
        longest = max(longest, *(len(t) for w in factor for t in w))
    assert (len(bank), len(factors), longest) == (20,4,8)
    return dict(seeds=len(bank), dressings=len(factors), largest_trace_in_slots=longest)


def exterior_product(left, right):
    out = Counter()
    for a, x in left.items():
        for b, y in right.items():
            if a & b:
                continue
            smaller = b
            sign = 1
            while smaller:
                bit = smaller & -smaller
                if (a // (2*bit)).bit_count() % 2:
                    sign = -sign
                smaller -= bit
            out[a | b] += sign*x*y
    return nonzero(out)


def determinant(matrix):
    output = Counter()
    for perm in permutations(range(len(matrix))):
        swaps = sum(a>b for i,a in enumerate(perm) for b in perm[i+1:])
        term = {0: (-1)**swaps}
        for i,j in enumerate(perm):
            term = exterior_product(term, matrix[i][j])
        add_scaled(output, term, 1)
    return nonzero(output)


def polarization_control():
    # Odd-variable order: six scalar-to-dot parameters, two marked parameters,
    # four positive dotted colour components, four negative dotted components.
    # The even parts of the three scalar columns are e1,e2,e3.
    matrix = []
    for row in range(4):
        entries = []
        for column in range(4):
            term = {0:1} if row == column and column < 3 else {}
            term[(1 << (2*column)) | (1 << (8+row))] = 1
            term[(1 << (2*column+1)) | (1 << (12+row))] = 1
            entries.append(term)
        matrix.append(entries)
    three = determinant([row[:3] for row in matrix[:3]])
    four = determinant(matrix)
    product = exterior_product(exterior_product(three,three), exterior_product(four,four))
    coefficient = product[(1 << 16)-1]
    assert coefficient
    return coefficient


def jet_control():
    def transfer(jet, flavour, dot):
        p, q, mask = jet
        bit = 1 << flavour
        if not mask & bit:
            return 0, jet
        sign = (-1)**(mask.bit_count()+(mask & (bit-1)).bit_count())
        return sign*(1+(p,q)[dot]), (p+int(dot==0),q+int(dot==1),mask^bit)

    jets = [(p,q,mask) for p,q,mask in product(range(3),range(3),range(8))
            if 1 <= p+q+mask.bit_count() <= 2]
    for jet, flavour, dot in product(jets, range(3), range(2)):
        factor, image = transfer(jet, flavour, dot)
        residual = Counter()
        for c, x, y in native_splits(image):
            residual[x,y] += factor*c
        for c, x, y in native_splits(jet):
            a, left = transfer(x, flavour, dot)
            b, right = transfer(y, flavour, dot)
            residual[left,y] += c*a
            residual[x,right] += (-1)**parity(x)*c*b
        assert not nonzero(residual)
    assert len(jets) == 17
    return len(jets)*6


def run(path):
    assert sha256(path.read_bytes()).hexdigest() == ARCHIVE
    with ZipFile(path) as archive:
        source = input_check(archive)
        data = json.loads((HERE/'dressed-tree.json').read_bytes())
        rank, kernel = verify(archive, data)
    source['projected_rank'] = rank
    source['parameter_kernel_dimension'] = kernel
    source['zero_before_quartic_reduction'] = sum(row['formal_terms']==0 for row in data['kernel'])
    source['nonempty_quartic_reductions'] = kernel - source['zero_before_quartic_reduction']
    return dict(archive_sha256=ARCHIVE, source=source,
                top_exterior_coefficient=polarization_control(),
                tree_transfer_anticommutators=jet_control(),
                stable_trace_ranks=[check(d) for d in range(1,9)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    output = json.dumps(run(args.archive), indent=2) + '\n'
    if args.out:
        args.out.write_text(output)
    print(output, end='')
