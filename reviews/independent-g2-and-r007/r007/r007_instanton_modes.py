"""R-charge support of isolated charge-one instanton profiles at the SU(3) R-007 Ward endpoints."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb, prod
from pathlib import Path
import gzip

import r007_primary_replay as primary


def load(name, expected):
    raw = (Path(__file__).resolve().parent/name).read_bytes()
    assert sha256(raw).hexdigest() == expected, name
    return raw


def occupancy(masks):
    # Slots 0,1,2,3 represent flavours 4,1,2,3 in that order.
    return tuple(sum((m >> i) & 1 for m in masks) for i in range(4))


def modes(weight, per_flavour):
    # The determinant character is trivial on SU(4); test ALL common shifts.
    base = [x-min(weight) for x in weight]
    return sum(prod(comb(per_flavour,x+j) for x in base)
               for j in range(per_flavour-max(base)+1))


def raw_source():
    raw = load('dmk_o330_hw.json','8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34')
    # SD: anti-scalar in the two complementary matter slots, chiral fermion in its own slot.
    masks = (12,10,6,2,4,8)
    contents = Counter()
    for c,multi in json.loads(raw)['terms']:
        letters = tuple(l for tr in multi for l in tr)
        assert Fraction(c) and len(letters)==10 and sum(l>=3 for l in letters)==5
        site = [masks[l] for l in letters]
        assert occupancy(site) == (0,3,6,6)
        assert occupancy([15 ^ m for m in site]) == (10,7,4,4)
        contents[tuple(letters.count(i) for i in range(6))] += 1
    assert sum(contents.values())==6139
    return masks, contents


def raw_target():
    archive = json.loads(gzip.decompress(load('su30_0.json.gz',
        '81f936300745c652c4279bb769365b59df01fdcb321860d657b7b6d3e1188ebc')))['selected_recipes']
    recipe = json.loads(load('R007_recipe.json',
        'c82ec0745328b8482c12dec843bb25d25345f9f7dffa6776720a37971acba4df'))['source_linear_combination']
    a = Counter()
    for i,c in recipe:
        b = {():Fraction(c)}
        for g in archive[i]:
            b = primary.multiply(b,primary.factor(tuple(g)))
        a.update(b)
    output = primary.clean(a)
    sizes = [len(output)]
    for x,y in ((0,1),(1,1),(0,2),(1,2)):
        output = primary.act(output,('T',x,y))
        sizes.append(len(output))
    assert sizes==[348,1059,1864,4561,5080]
    contents = Counter()
    bmn_sets = ((0,), (1,), (2,), (1,2), (0,2), (0,1))
    for multi in output:
        letters = [l for tr in multi for l in tr]
        # A self-dual gaugino, scalar, matter fermion, curvature gives respectively
        # 3,2,1,0 minimum modes. The first letter here is the antichiral gaugino
        # of the defining BPS alphabet in a self-dual background.
        masks = [14 ^ sum(2 << i for i in l[2]) for l in letters]
        assert len(letters)==9 and occupancy(masks)==(0,4,7,7)
        assert occupancy([15^m for m in masks])==(9,5,2,2)
        assert [sum(l[i] for l in letters) for i in (0,1)] == [2,2]
        content = tuple(sum(l[2]==key for l in letters) for key in bmn_sets)
        contents[content]+=1
    return sizes, contents


def wick_screens(masks, src, target):
    # Normal-ordered intersite pairs between equal field species, with the right
    # insertion radially conjugated. Each bosonic-background matter propagator
    # removes one quartet from the sum of the endpoint minima. This is a
    # combinatorial availability bound only; there are no propagator evaluations.
    tests = {}
    for name, required, sign in [('self_dual',3,0),('anti_self_dual',5,1)]:
        available = 0
        total = sum(src.values())*sum(target.values())
        for left,n in src.items():
            for right,h in target.items():
                caps = tuple(min(a,b) for a,b in zip(left,right))
                if sum(caps)<required:
                    continue
                for taken in product(*(range(1+v) for v in caps)):
                    if sum(taken)!=required:
                        continue
                    rem = [m ^ (15 if sign else 0) for m,a,b in zip(masks,left,taken)
                           for _ in range(a-b)]
                    # This check uses the global bound only and makes no
                    # assumptions on the detailed spin dependence at either site.
                    if max(occupancy(rem))<=6:
                        available+=n*h
                        break
        tests[name] = dict(mandatory_quartets=required,raw_word_pairs=total,
                           word_pairs_passing_species_and_source_capacity=available)
    return tests


def run():
    masks, src = raw_source()
    sizes, target = raw_target()
    # Independent enumeration of all 7^4 occupation profiles of the full
    # 24-dimensional Grassmann space, with binomial weight multiplicities.
    exterior = Counter()
    for occ in product(range(7),repeat=4):
        key = tuple(x-min(occ) for x in occ)
        exterior[key] += prod(comb(6,x) for x in occ)
    assert sum(exterior.values())==2**24
    answers = {}
    for name,w in [('source',(0,3,6,6)),('target',(0,4,7,7))]:
        dual = tuple(max(w)-x for x in w)
        answers[name] = [modes(a,6) for a in (w,dual)]
        assert answers[name]==[exterior[a] for a in (w,dual)]
    assert answers==dict(source=[20,20],target=[0,0])
    half = [(11-r,7-r,4-r,4-r) for r in range(5)]
    assert all(modes(w,6)==0 for w in half)
    return dict(source_terms=sum(src.values()),target_terms=sizes[-1],target_rebuild=sizes,
                flavour_order=[4,1,2,3],modes_per_flavour=6,full_exterior_dimension=2**24,
                exterior_weight_multiplicities=answers,
                source_bare_minima=[[0,3,6,6],[10,7,4,4]],
                target_grade_two_minima=[[0,4,7,7],[9,5,2,2]],
                current_minima=[[0,1,1,1],[1,0,0,0]],
                uncontracted_Ward_minima=[[9]*4,[11]*4],
                matter_pair_content_tests=wick_screens(masks,src,target))


if __name__=='__main__':
    print(json.dumps(run(), indent=2))
