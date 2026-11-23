"""Deterministic two-internal-vertex Pfaffian graph identities.

This extends the existing path recipes, without changing field or Q0
conventions. All original source caches are read-only.
"""
from pathlib import Path
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import json
import math
import pickle
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
JETS = HERE.parent / 'fortuity_discovery_20261123/jets'
sys.path.insert(0, str(JETS))
import path_relation_completion as pc
import pfaffian_relations as pf
from flint import fmpz_mat


@lru_cache(None)
def ordered_subwords(letters, length):
    result = set()
    for indices in combinations(range(len(letters)), length):
        selected = tuple(letters[i] for i in indices)
        rest = tuple(l for i, l in enumerate(letters) if i not in indices)
        for order in set(permutations(selected)):
            result.add((order, rest))
    return tuple(sorted(result))


def graph_basis(contents, topology):
    result = set()
    for content in contents:
        letters = tuple(sorted(content))
        if topology == 'one4':
            for path, rest in ordered_subwords(letters, 4):
                for pairs in pc.paired_content(rest):
                    sign, paths = pc.canonical_paths((path,)+pairs)
                    if sign:
                        result.add(paths)
        elif topology == 'two3':
            for left, rest in ordered_subwords(letters, 3):
                for right, leftover in ordered_subwords(rest, 3):
                    if right < left:
                        continue
                    for pairs in pc.paired_content(leftover):
                        sign, paths = pc.canonical_paths((left,right)+pairs)
                        if sign:
                            result.add(paths)
        elif topology == 'cycle2':
            for cycle,rest in ordered_subwords(letters,2):
                sign,word=pf.m.so_canon((cycle,))
                if not sign:
                    continue
                for pairs in pc.paired_content(rest):
                    sign,paths=pc.canonical_paths(pairs)
                    if sign:
                        result.add((paths,word[0]))
        else:
            raise ValueError(topology)
    return sorted(result)


def raise_cycle(graph):
    paths,cycle=graph
    result=Counter()
    for new_paths,c in pc.raising_paths(paths).items():
        result[new_paths,cycle]+=c
    for word,c in pf.spin_action({(cycle,):1},True).items():
        result[paths,word[0]]+=c
    return {w:c for w,c in result.items() if c}


def expand_cycle(graph,rank):
    paths,cycle=graph
    result=Counter()
    for word,c in pc.expand_paths(paths,rank).items():
        sign,new_word=pf.m.so_canon(word+(cycle,))
        if sign:
            result[new_word]+=sign*c
    return {w:c for w,c in result.items() if c}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rank', type=int, default=4)
    parser.add_argument('--seconds', type=int, default=300)
    parser.add_argument('--topology', choices=['one4','two3','cycle2'], default='one4')
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('target completion cap')))
    signal.alarm(args.seconds)
    rank = args.rank
    length = 2 * rank + 4
    start = time.monotonic()
    output = HERE / f'complete_target_sp{rank}.pkl'
    source = JETS / ('target_snapshot_181.pkl' if rank == 4 else 'target_sp5_snapshot_181.pkl')
    rows, recipes = pickle.loads((output if output.exists() else source).read_bytes())
    echelon = pf.Echelon()
    for row in rows:
        assert echelon.add(row)
    dims, contents = pf.profile_inventory(rank, (3,3,2*rank-1,0,0), length)
    by_profile = {profile:pf.Echelon() for profile in dims}
    for row in rows:
        for profile, part in pf.profile_components(row).items():
            by_profile[profile].add(part)
    missing = [profile for profile,n in dims.items() if len(by_profile[profile].rows)<n]
    receipt = dict(rank=rank, native_weight=[3,3,2*rank-1,0,0], length=length,
                   exact_character_dimension=sum(dims.values()), starting_dimension=len(rows), additions=[])
    recipe_type='paths_cycle' if args.topology=='cycle2' else 'paths'
    memoized_expand = lru_cache(maxsize=512)(expand_cycle if args.topology=='cycle2' else pc.expand_paths)
    for profile in missing:
        basis = graph_basis(contents[profile], args.topology)
        images = [(raise_cycle(paths) if args.topology=='cycle2' else pc.raising_paths(paths)) for paths in basis]
        target = sorted({w for image in images for w in image})
        positions = {w:i for i,w in enumerate(target)}
        raising = fmpz_mat(len(target),len(basis))
        for j,image in enumerate(images):
            for w,c in image.items():
                raising[positions[w],j]=c
        kernel, dimension = raising.nullspace()
        print(json.dumps(dict(topology=args.topology, graph_basis=len(basis), raised=len(target),
                              highest_weights=dimension, elapsed_seconds=time.monotonic()-start)),flush=True)
        old_recipes = {tuple(r[1]) for r in recipes if r[0] == recipe_type}
        for j in range(dimension):
            coefficients = [int(kernel[i,j]) for i in range(len(basis))]
            divisor = math.gcd(*coefficients)
            recipe = tuple((basis[i],c//divisor) for i,c in enumerate(coefficients) if c)
            if recipe in old_recipes:
                continue
            polynomial = Counter()
            for paths,c in recipe:
                for w,n in memoized_expand(paths,rank).items():
                    polynomial[w] += c*n
            row = pf.primitive({w:c for w,c in polynomial.items() if c})
            assert not pf.spin_action(row,True)
            if row and echelon.add(row):
                assert by_profile[profile].add(row)
                rows.append(row)
                recipes.append((recipe_type,recipe))
                output.write_bytes(pickle.dumps((rows,recipes),protocol=4))
                addition = dict(kernel_column=j, identity_index=len(rows)-1, trace_terms=len(row),
                                path_terms=len(recipe),elapsed_seconds=time.monotonic()-start)
                receipt['additions'].append(addition)
                print(json.dumps(addition),flush=True)
            if j % 10 == 0:
                print(json.dumps(dict(column=j,dimension=len(rows),elapsed_seconds=time.monotonic()-start)),flush=True)
            if len(by_profile[profile].rows)==dims[profile]:
                break
    receipt.update(complete=len(rows)==sum(dims.values()), relation_dimension=len(rows),
                   elapsed_seconds=time.monotonic()-start)
    (HERE/f'complete_target_sp{rank}_{args.topology}.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    main()
