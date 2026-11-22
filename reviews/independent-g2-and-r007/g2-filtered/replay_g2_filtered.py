"""Replay the two finite-field correction equations and regenerate the underlying jet maps."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from flint import nmod_mat

from g2_first_primitive import read

ROOT = Path(__file__).resolve().parent
BALANCED = ((2147483647, 731099, 48031520), (1000033, 902117, 439099))


def field(a, p):
    return nmod_mat(a.tolist(), p)


def check(p, seed, minor):
    with np.load(ROOT / f'coefficients-{p}.npz', allow_pickle=False) as v:
        def m(key):
            return field(v[key], p)
        assert list(v['spec']) == [p, seed]
        first, second = m('first'), m('second')
        first_matrix = m('first_matrix')
        assert first.nrows() == 944 and np.count_nonzero(v['first']) == 944
        selected = v['first_matrix'][v['minor_rows']]
        assert int(field(selected, p).det()) == minor
        assert first_matrix * m('first_correction') == -first
        assert second.nrows() == 115 and np.count_nonzero(v['second']) == 115
        before, after, ev = m('before'), m('after'), m('eval')
        assert (before.nrows(), before.rank(), after.rank()) == (115, 101, 13)
        assert after * ev.solve(before) == nmod_mat(13, before.ncols(), p)
        assert before * m('second_correction') == -second
        dual = m('dual')
        assert dual * before == nmod_mat(1, before.ncols(), p)
        assert dual * second == nmod_mat(1, 1, p)
        assert dual * ev * m('cycle').transpose() != nmod_mat(1, 1, p)
        assert after * m('cycle').transpose() == nmod_mat(13, 1, p)
        assert v['homogeneous'].shape == (1, 46) and not np.any(v['homogeneous'])
        assert np.count_nonzero(v['first_correction']) == 796
        assert np.count_nonzero(v['second_correction']) == 101
    print(f'{p}: 944 first and 115 second coordinates cancel; 46 homogeneous pairings vanish.', flush=True)


def rebuild(repo, scratch, workers, primes):
    from build_g2_full_chain import build as classical_build
    from build_g2_filtered import build as quantum_build
    native = scratch / 'native'
    native.mkdir()
    classical_build(repo, native)
    quantum_build(repo, native, ROOT)
    for name, source in [('sector','g2_filtered_sector.cpp'), ('page','g2_page_basis.cpp'),
                         ('products','g2_product_span.cpp')]:
        subprocess.run(['c++','-O3','-std=c++17','-pthread','-I'+str(native),
                        ROOT/source,'-o',native/name],check=True)

    def execute(command, *args, weight='5 2 2 1 1'):
        env = dict(os.environ, G2_SOURCE_WEIGHT=weight, Q1_MUTATE='0', Q1_CMUT='0')
        result = subprocess.run([str(command),*map(str,args)],env=env,check=True,
                                stdout=subprocess.PIPE,text=True)
        return result.stdout.strip()

    def script(name, *args):
        return execute(sys.executable, ROOT/name, *args)

    def compare_packet(directory, p):
        with np.load(ROOT/f'coefficients-{p}.npz', allow_pickle=False) as ref:
            def get(name):
                q, ary = read(directory/(name+'.u32'))
                assert q == p
                return ary
            piv = np.flatnonzero(get('primitive_1_8')[0])
            assert np.array_equal(get('in_1_9')[:,piv],ref['first_matrix'])
            assert np.array_equal(get('primitive_1_8')[:,piv].T,ref['first_correction'])
            for stored, generated in dict(first='acted_cycle_0_9',second='acted_primitive_1_8',
                    before='in_2_8',after='in_2_9',eval='eval_2_8',cycle='cycle_2_8',
                    dual='dual_2_8',homogeneous='paired_cycle_1_8',source='cycle_0_9').items():
                assert np.array_equal(ref[stored],get(generated)), (p,stored)
            assert np.array_equal(ref['second_correction'],get('primitive_2_7').T)

    for p, seed, _ in BALANCED:
        if p not in primes:
            continue
        folder=scratch/str(p)
        for grade, lo, hi, trim in [(0,8,10,1),(1,7,10,0),(2,7,9,0)]:
            execute(native/'sector',folder,grade,lo,hi,p,seed,workers,trim)
        for grade, length, sizes in [(0,9,(2872,1825,1046,1)),(1,8,(2124,1282,796,46)),
                                     (1,9,(944,796,148,0)),(2,8,(115,101,13,1))]:
            row=json.loads(execute(native/'page',folder,grade,length,workers))
            assert tuple(row[k] for k in ('dimension','incoming','outgoing','cohomology')) == sizes
        for grade,length,vec,dual in [(1,8,'cycle',1),(0,9,'cycle',0)]:
            execute(native/'g2_filter_q1',folder,grade,length,vec,seed,p,workers,dual)
        script('g2_first_primitive.py',folder)
        execute(native/'g2_filter_q1',folder,1,8,'primitive',seed,p,workers,0)
        script('g2_first_primitive.py',folder)
        compare_packet(folder,p)

    p=2147483647
    for weight,grade,seed,ranks in [('5 2 2 1 1',-1,731099,(733,307,426)),
                                    ('5 2 2 2 0',0,120781,(2986,1926,1060)),
                                    ('3 3 2 3 0',0,703019,(4345,3098,1247))]:
        directory=scratch/('classical-'+weight.replace(' ',''))
        execute(native/'sector',directory,grade,8,10,p,seed,workers,1,weight=weight)
        row=json.loads(execute(native/'page',directory,grade,9,workers,weight=weight))
        assert tuple(row[k] for k in ('dimension','incoming','outgoing'))==ranks
        assert ranks[0]==sum(ranks[1:])
    directory=scratch/'outgoing'
    execute(native/'sector',directory,0,10,11,p,731099,workers,1)
    _,a=read(directory/'in_0_11.u32')
    assert a.shape==(176,1222) and field(a,p).rank()==176
    for tag,w,seed,expected in [('522',(5,2,2,1,1),691447,(1331,1390)),
                                 ('241',(4,1,1,2,2),448117,(1282,1328))]:
        words=scratch/(tag+'.dat')
        script('g2_product_words.py','--repo',repo,'--recipes',ROOT/f'g2-{tag}-recipes.json.gz',
               '--out',words)
        row=json.loads(execute(native/'products',*w,8,seed,workers,p,words))
        assert (row['boundary'],row['augmented'])==expected
    print('Regenerated the classical ranks, product spans, and specified finite-field corrections.',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true')
    parser.add_argument('--repo',type=Path)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--prime',type=int,choices=[x[0] for x in BALANCED],action='append')
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        raise ValueError('enable assertions and select a positive worker count')
    for p,seed,minor in BALANCED:
        check(p,seed,minor)
    if args.rebuild:
        if not args.repo:
            parser.error('--rebuild requires --repo')
        with tempfile.TemporaryDirectory() as scratch:
            rebuild(args.repo.resolve(),Path(scratch),args.workers,args.prime or [x[0] for x in BALANCED])
