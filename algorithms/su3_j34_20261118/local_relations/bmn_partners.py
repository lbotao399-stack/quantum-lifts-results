"""Resolve the complete grade-zero classical partner spaces of the level-34 bosonic seed."""
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'dual_singlets_20261118'))
from rank_two_bmn_product_replay import highest, witness
from sp3_omega_products import weight as factor_weight
from su3_bmn_parity_products import labelled
from trace_descent import dimensions


def verify():
    weight = (8,5,4)
    raw = gzip.decompress((HERE/'public-bmn-column.json.gz').read_bytes())
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest() == 'ea35527b4bac2c2b4fdee4909d665eced216a71c'
    original = json.loads(raw)
    assert (original['N'],original['j'],tuple(original['lam'])) == (3,34,weight)
    records = json.loads(gzip.decompress((HERE/'bmn-partner-products.json.gz').read_bytes()))
    lower = {}
    for record in records:
        y = record['length']
        recipes = [tuple(map(tuple,r)) for r in record['inputs']]
        for recipe in recipes:
            assert sum(sum(g[:5]) for g in recipe) == y
            assert tuple(map(sum,zip(*map(factor_weight,recipe)))) == (0,0,*weight)
        cyclic,ordinary = [witness([highest(labelled(r,ordinary_diagonal=b),weight) for r in recipes])
                           for b in (False,True)]
        assert json.loads(json.dumps(cyclic)) == record['minor']
        assert ordinary['rank'] == cyclic['rank'] == len(recipes) == record['rational_product_rank']
        lower[y] = len(recipes)
    assert lower == {9:1,10:10,11:16,12:6,13:1}

    previous = 0
    upper = {}
    euler = [0,0]
    for n,row in enumerate(original['sectors']):
        y = sum(weight)-n
        cochains = [sum(dimensions((0,0,*weight),y,3,c).values()) for c in ('SU+','SU-')]
        assert n == row['n'] and row['check'] and sum(cochains) == row['m']
        upper[y] = sum(cochains)-previous-row['r_up']
        previous = row['r_up']
        assert upper[y] == row['h']
        for i in range(2):
            euler[i] += (-1)**y*cochains[i]
    assert set(upper) == set(range(8,18)) and euler == [14,-17]
    lengths = sorted(upper)
    options = []
    for y in lengths:
        pairs = [(a,b) for a in range(upper[y]+1) for b in range(upper[y]+1-a)
                 if (a,b)[y%2] >= lower.get(y,0)]
        options.append(pairs)
    solutions = [h for h in product(*options) if all(
        sum((-1)**y*pair[i] for y,pair in zip(lengths,h)) == euler[i] for i in (0,1))]
    h, = solutions
    result = {y:v for y,v in zip(lengths,h) if any(v)}
    assert result == {9:(0,1),10:(10,1),11:(2,16),12:(6,0),13:(0,1)}
    return dict(colour_euler=euler,cohomology_even_odd=result,protected_Cartan_products=lower)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('enable assertions')
    print(json.dumps(verify(),indent=2))
