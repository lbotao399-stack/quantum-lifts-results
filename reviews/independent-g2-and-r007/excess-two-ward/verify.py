"""Rational reverse-shear test at Taylor excess two using separating matrix units."""
import argparse
from functools import cache
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path
import sys


def matrix_colour(rw, n):
    # e_ij has trace dual e_ji. Colour indices in reverse_ward label coordinates.
    units = list(product(range(n), repeat=2))
    metric = [[int(i == l and j == k) for k, l in units] for i, j in units]
    brackets = []
    for i, j in units:
        row = []
        for k, l in units:
            value = {}
            if j == k:
                rw.add(value, i*n+l, 1)
            if i == l:
                rw.add(value, k*n+j, -1)
            row.append(value)
        brackets.append(row)
    return rw.Brackets(brackets, metric, metric)


def forward(rw, word, spin, flavour):
    out = {}
    sign = 1
    for i, (letter, colour) in enumerate(word):
        p, q, fs = letter
        if flavour in fs:
            z = list((p,q))
            z[spin] += 1
            jacobian = (-1)**sum(k > flavour for k in fs) * z[spin]
            new = ((*z, tuple(k for k in fs if k != flavour)), colour)
            rw.term(out, word[:i]+(new,)+word[i+1:], sign*jacobian)
        if rw.is_ferm(letter):
            sign = -sign
    return {w:v for (_,w),v in out.items()}


def action(rw, w, coeff):
    @cache
    def pair(word, correction):
        out = {}
        for (key,target), v in w.two(word,correction).items():
            rw.add(out,target,v*(coeff.get(key,0) if correction else 1))
        return out
    return pair


def anticommutator(rw, first, second, word):
    out = {}
    individual = 0
    for left,right in ((first,second),(second,first)):
        chain = {}
        for mid,c in right(word).items():
            for end,d in left(mid).items():
                rw.add(chain,end,c*d)
        individual += len(chain)
        for end,c in chain.items():
            rw.add(out,end,c)
    return out, individual


def main(repo):
    root = repo / "algorithms/su3_j30_classical_ranks/bmn_return"
    sys.path.insert(0,str(root))
    import reverse_ward as rw
    native = rw.cross_product()
    colours = matrix_colour(rw,5)
    cycle = [5*i+(i+1)%5 for i in range(5)]
    first_pair_colours = [(0,0),(0,1)]
    all_records = []
    for spin,flavour in product(range(2),range(3)):
        small = rw.Ward(native,spin,flavour)
        count, rank, coeff = rw.system(small,rw.examples(5,(0,),first_pair_colours))
        checked,_,_ = rw.system(small,rw.examples(5,range(3),list(product(range(3),repeat=2))),coeff)
        unit_ward = rw.Ward(colours,spin,flavour)
        pair = action(rw,unit_ward,coeff)
        q1 = lambda word: pair(word,False)
        correction = lambda word: pair(word,True)
        fwd = lambda word: forward(rw,word,spin,flavour)
        count_fwd = nontrivial_fwd = 0
        # A connected colour pair is bilinear, so its two inputs may be fixed
        # to two different edges of a four-cycle. All output components are kept.
        four = [1,7,13,15]
        for word in rw.examples(5,(),[(four[0],j) for j in four[1:]]):
            residual,individual = anticommutator(rw,q1,fwd,word)
            assert not residual, (spin,flavour,word,residual)
            count_fwd+=1
            nontrivial_fwd+=individual
        raw_tokens = [l for l in rw.letters(3) if l[:2]==(0,0)]
        cases = []
        for types in combinations_with_replacement(raw_tokens,3):
            if any(sum(i in l[2] for l in types) != 2-int(i==flavour) for i in range(3)):
                continue
            tested = nonzero = 0
            for j,k in permutations(cycle[1:],2):
                sg,word = rw.canonical(tuple(zip(types,(cycle[0],j,k))))
                assert sg
                residual,individual = anticommutator(rw,q1,correction,word)
                assert not residual, (spin,flavour,word,residual)
                tested+=1
                nonzero+=individual
            cases.append(dict(letters=[list(l[2]) for l in types], tests=tested,
                              separate_nonzero_components=nonzero))
        assert len(cases)==4 and sum(row["separate_nonzero_components"] for row in cases)>0
        all_records.append(dict(spin=spin,flavour=flavour,first_equations=count,
                                first_pivots=rank,first_expanded_equations=checked,
                                forward_pairs=count_fwd,forward_separate_components=nontrivial_fwd,
                                triple_patterns=cases))
    def par(n,maximum=None):
        if n==0:
            yield ()
        else:
            for k in range(min(n, maximum if maximum is not None else n),0,-1):
                for tail in par(n-k,k):
                    yield (k,)+tail
    possible = {}
    for n in range(1,15):
        possible[n]=[]
        for shape in par(n):
            if len(shape)>=3 and shape[2]>3:
                continue
            odd = tuple(max(0,sum(v>=i for v in shape)-2) for i in (1,2,3))
            if min(odd)>=2:
                possible[n].append(shape)
    assert not any(possible[d] for d in range(1,12))
    assert possible[12]==[(3,3,3,3)]
    return dict(first_and_second_shear_tests=all_records,
                admissible_source_shapes={str(n):v for n,v in possible.items() if v})


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(main(args.repo.resolve()),indent=2))
