"""Check the rational rank bounds, the Cartan character and the fourfold return."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import product
from rank_two_bmn_index import colour_profiles, multiplicity
from rank_two_bmn_product_replay import highest, witness
from sp3_omega_products import weight as factor_weight
from su3_bmn_parity_products import labelled
from trace_descent import dimensions

HERE = Path(__file__).resolve().parent
WEIGHT = (5,4,4)
PRIMES = (1000003,1000033)
SEEDS = (891033,891059)


def read_minor(record, name, height, width):
    minor = record[name]
    r = minor['rank']
    assert 0 < minor['det'] < record['prime']
    for indices,bound in ((minor['rows'],height),(minor['columns'],width)):
        assert len(set(indices)) == len(indices) == r
        assert all(0 <= i < bound for i in indices)
    return r


def check(blocks):
    chain = {y:[sum(dimensions((0,0,*WEIGHT),y,3,c).values()) for c in ('SU+','SU-')]
             for y in range(5,14)}
    ranks = {4:[0,0],13:[0,0]}
    for y in range(5,13):
        row = blocks['colour'][str(y)]
        assert (row['mode'],row['gauge'],row['sector'],row['length'],row['cochains']) == (
            'colour_image','SU3',[5,4,4,0,0],y,sum(chain[y]))
        assert row['prime'] == 1000033 and row['seed'] == 901261+13*y
        ranks[y] = [read_minor(row,c,row['target_evaluations'],row['cochains']) for c in ('even','odd')]
    upper = [[chain[y][i] - ranks[y-1][i] - ranks[y][i] for y in range(5,14)] for i in range(2)]
    euler = [sum((-1)**y * chain[y][i] for y in range(5,14)) for i in range(2)]
    assert euler == [2,-1] and upper == [[0,0,0,2,0,0,0,0,0],[0,0,0,0,1,0,0,0,0]]
    cartan = sum(multiplicity(WEIGHT,p) for p in colour_profiles(WEIGHT,9,colour_group='SU3'))
    assert cartan == 0
    determinants = []
    for p,seed in zip(PRIMES,SEEDS,strict=True):
        row = blocks['products'][str(p)]
        assert (row['mode'],row['gauge'],row['prime'],row['seed'],row['length'],row['sector']) == (
            'explicit_classical_image','SU3',p,seed,9,[5,4,4,0,0])
        a,b,c = row['cochains']
        assert [a,b,c] == [sum(chain[y]) for y in (8,9,10)] == [971,1360,1117]
        assert read_minor(row,'incoming',b,a) == 658 and read_minor(row,'outgoing',c,b) == 701
        assert read_minor(row,'augmented',b,a+1) == 659 and a in row['augmented']['columns']
        determinants.append(row['augmented']['det'])
    assert determinants == [563838,242547]
    return dict(euler=euler, bmn_cohomology=b-658-701, bmn_Cartan_dimension=cartan,
                augmented_determinants=determinants)


def check_source(records):
    flavour = (7,4,4)
    dims = [sum(dimensions((0,0,*flavour),10,3,c).values()) for c in ('SU+','SU-')]
    assert dims == [2768,2896]
    ranks = []
    for length in (9,10):
        r = records[str(length)]
        assert (r['mode'],r['gauge'],r['length'],r['sector'],r['prime'],r['seed']) == (
            'colour_image','SU3',length,[*flavour,0,0],1000003,913330+length)
        ranks.append([read_minor(r,c,r['target_evaluations'],r['cochains']) for c in ('even','odd')])
    assert ranks == [[1113,1073],[1651,1822]]
    upper = [d-i-o for d,i,o in zip(dims,*ranks,strict=True)]
    assert upper == [4,1]
    products = json.loads((HERE/'source-products.json').read_bytes())
    recipes = [tuple(map(tuple,r)) for r in products['inputs']]
    assert products['rational_product_rank'] == len(recipes) == 4 and products['length'] == 10
    for recipe in recipes:
        assert tuple(map(sum,zip(*map(factor_weight,recipe)))) == (0,0,*flavour)
        assert sum(sum(g[:5]) for g in recipe) == 10
    cyclic,ordinary = [witness([highest(labelled(r,ordinary_diagonal=b),flavour) for r in recipes])
                       for b in (False,True)]
    assert json.loads(json.dumps(cyclic)) == products['minor']
    assert cyclic['signed_integer_determinant'] == '2812385772'
    assert ordinary['signed_integer_determinant'] == '3528'
    return dict(source_even_odd_upper=upper,ordinary_Cartan_determinant=3528)


def rebuild(folder, native_product, blocks, sources, threads):
    taylor = product.ALG/'dual_singlets_20261118/taylor_hooks'
    env = dict(os.environ, PAR_THREADS=str(threads))
    def compile_and_run(program, jobs):
        exe = folder/program
        subprocess.run(['c++','-std=c++17','-O3','-pthread','-I'+str(taylor/'native'),
                        str(taylor/(program+'.cpp')),'-o',str(exe)],check=True)
        for args,expected in jobs:
            result = subprocess.check_output([str(exe),*map(str,args)],env=env)
            assert json.loads(result) == expected, args
    compile_and_run('unitary_parity', [
        (['SU3',*WEIGHT,0,0,y,901261+13*y,1000033,'image'], blocks['colour'][str(y)])
        for y in range(5,13)] + [
        (['SU3',7,4,4,0,0,y,913330+y,1000003,'image'], sources[str(y)]) for y in (9,10)])
    compile_and_run('su2_levi_image', [
        ([native_product,*WEIGHT,0,0,9,p,seed,'SU3'],blocks['products'][str(p)])
        for p,seed in zip(PRIMES,SEEDS,strict=True)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild-native', action='store_true')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if not __debug__ or args.threads < 1:
        parser.error('enable assertions and use a positive thread count')
    raw = (HERE/'native-minors.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == 'a3fc9b252da812eb8272130f227028c1db710700bf653000d9cf4eab2207e048'
    blocks = json.loads(gzip.decompress(raw))
    raw_sources = (HERE/'source-minors.json.gz').read_bytes()
    assert hashlib.sha256(raw_sources).hexdigest() == '2f95a5ad26786676993d0fbfe6ff29094b3b51a188995d5f3837f16f1424c65b'
    sources = json.loads(gzip.decompress(raw_sources))
    summary = check(blocks)
    summary.update(check_source(sources))
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        native_product = folder/'bmn.txt'
        result = product.verify(native_product)
        assert json.loads(json.dumps(result)) == json.loads((HERE/'product.json').read_bytes())
        if args.rebuild_native:
            rebuild(folder,native_product,blocks,sources,args.threads)
    print(json.dumps({**summary, 'terms':result['terms_after_substitutions'],
                      'return_coefficient':result['return_coefficient']}, indent=2))
