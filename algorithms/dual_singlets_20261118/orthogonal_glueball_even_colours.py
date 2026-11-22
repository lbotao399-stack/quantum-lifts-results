"""Test polarized Pfaffian relations for a second differential into an orthogonal glueball."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product
from multiprocessing import get_context
from random import Random

from orthogonal_glueball_relation_span import determinant, relation_dimension, weights
from orthogonal_glueball_second_pair import second_dual
from orthogonal_glueball_second_sources import ranks
from terminal_commutators import words


@lru_cache(None)
def path_shapes(colours, power):
    excess = 2*power + 1 - colours
    # A path with k letters uses k-1 contractions outside the Pfaffian.
    available = sorted({min(p, p[::-1]) for k in range(2, excess+2)
                        for chars in product('mps', repeat=k) for p in [''.join(chars)]})
    shapes = []
    for count in range(1, excess+1):
        for paths in combinations_with_replacement(available, count):
            if sum(len(p)-1 for p in paths) > excess:
                continue
            letters = ''.join(paths)
            if letters.count('s') <= 6 and all(0 <= power-2-letters.count(s) <= 3 for s in 'pm'):
                shapes.append(paths)
    return tuple(shapes)


def polarize(index, colours, power):
    shapes = path_shapes(colours,power)
    paths = shapes[index % len(shapes)]
    rng = Random(427021 + 10007*colours + 104729*index)
    matrices = [weights(rng,4) for _ in paths]
    if (index//len(shapes)) % 3 == 1:
        chosen = rng.randrange(len(paths))
        if len(paths[chosen]) % 2:
            matrices[chosen] = [0,0,0,1]
        else:
            matrices[chosen] = weights(rng,3)+[0]
    return dict(paths=paths, matrices=matrices,
                single=[rng.sample(range(3),power-2-''.join(paths).count(s)) for s in 'pm'],
                scalar_directions=[[weights(rng,3) for _ in range(p.count('s'))] for p in paths])


def vanishes_on_planes(spec):
    for path, auxiliary in zip(spec['paths'],spec['matrices'],strict=True):
        if path.count('m') >= 2 or path.count('p') >= 2:
            return True
        if (len(path)%2==0 and not auxiliary[3]) or (len(path)%2==1 and not any(auxiliary[:3])):
            return True
    return False


def evaluate(task):
    index, spec, colours, power, scores = task
    relation = determinant(power,spec)
    assert relation.pf_degree > colours
    basis = words((power-2,power-2,2,2,2))
    zero = vanishes_on_planes(spec)
    cartan = [str(relation.cartan(w)) if not zero and len(w) > colours//2 else '0' for w in basis]
    score = sum(c*relation.coefficient(w) for w,c in scores)
    relation.clear()
    return dict(index=index, spec=spec, pfaffian_degree=relation.pf_degree,
                forced_zero_on_planes=zero, cartan=cartan, second_score=str(score))


def search(colours, power, workers, attempts, save):
    lower = [relation_dimension(colours,power,k,False) for k in (0,1)]
    if lower != [0,0]:
        raise ValueError('This search requires the full earlier weights to be relation-free')
    dimension = relation_dimension(colours,power,2)
    dual = second_dual(colours,power)
    scores = [(tuple(map(tuple,w)),Fraction(v)) for w,v in dual['scores']]
    rows = []
    record = dict(colours=colours,power=power,lower_full_weight=lower,flavour_relations=dimension,
                  path_shapes=len(path_shapes(colours,power)),dual=dual,rows=rows)
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for start in range(0,attempts,workers):
            tasks = [(i,polarize(i,colours,power),colours,power,scores)
                     for i in range(start,min(start+workers,attempts))]
            rows.extend(pool.map(evaluate,tasks))
            record.update(ranks(rows))
            save(record)
            if 'witness' in record:
                break
    return record
