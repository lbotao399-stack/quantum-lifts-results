"""Test the second pairing on a character-complete finite-colour relation space."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod
from random import Random

from b2_schur_complex import centralizer, character, partitions
from mixed_incoming import so_canon
from orthogonal_glueball_cartan_basis import matrix
from orthogonal_glueball_first_pair import GHOSTS, PHI
from orthogonal_glueball_second_sources import Relation
from orthogonal_glueball_sources import missing_adjoint_trace
from orthogonal_glueball_third_pair import encode, targets


def relation_dimension(colours, power, grade, flavour_singlet=True):
    """Dimension lost on specialization, with the two dotted weights fixed."""
    ghosts = power - grade
    if flavour_singlet:
        scalars = [(z, Fraction(character((grade,)*3, z), centralizer(z)))
                   for z in partitions(3*grade)]
    else:
        scalars = [(tuple(sorted(a+b+c)), Fraction(1, prod(map(centralizer, (a,b,c)))))
                   for a,b,c in product(partitions(grade), repeat=3)]
    count = Fraction()
    for x,y in product(partitions(ghosts), repeat=2):
        weight = Fraction((-1)**(len(x)+len(y)), centralizer(x)*centralizer(y))
        for z, multiplicity in scalars:
            if multiplicity:
                count += weight*multiplicity*missing_adjoint_trace(tuple(sorted(x+y+z)), colours)
    assert count.denominator == 1 and count >= 0
    return int(count)


def path_shapes(power, colours):
    ghosts = power-2
    shapes = []
    for lengths in ((2,), (3,), (2,2)):
        if 2*power+2-sum(k-1 for k in lengths) <= colours:
            continue
        shapes.extend(paths for paths in product(*(product('mps', repeat=k) for k in lengths))
                      if all(0 <= ghosts-sum(path.count(s) for path in paths) <= 3 for s in 'pm'))
    return shapes


def weights(rng, n):
    row = [rng.randrange(-2,3) for _ in range(n)]
    if not any(row):
        row[rng.randrange(n)] = 1
    return row


def specification(index, power, colours):
    paths = path_shapes(power, colours)
    rng = Random(7902301 + 104729*index + colours)
    paths = paths[index % len(paths)]
    return dict(paths=paths,
                single=[rng.sample(range(3), power-2-sum(path.count(s) for path in paths))
                        for s in 'pm'],
                matrices=[weights(rng, 4) for path in paths],
                scalar_directions=[[weights(rng,3) for _ in range(path.count('s'))] for path in paths])


def determinant(power, specification):
    """Insert W+W^t or W-W^t in each marked path of the polarized Pfaffian."""
    plus, minus = specification['single']
    relation = Relation(power, dict(pattern=(), indices=(), plus=plus, minus=minus))
    for k,(path,matrices,directions) in enumerate(zip(specification['paths'],
                                                      specification['matrices'],
                                                      specification['scalar_directions'])):
        odd = len(path)-path.count('s')
        reverse = (-1)**(len(path)+odd*(odd-1)//2)
        bit = 1 << (6+k) if odd % 2 else 0
        relation.mask |= bit
        relation.used |= 1 << k
        relation.pf_degree -= len(path)-1
        options = []
        axes = iter(directions)
        for slot in path:
            if slot == 's':
                axis = next(axes)
                options.append([(i,1 << (3*(3*i+j)),axis[j])
                                for i,j in product(range(3), repeat=2) if axis[j]])
            else:
                options.append([(3 if slot == 'p' else 4, 0, 1)])
        entries = Counter()
        for assignment in product(*options):
            word = tuple(i for i,_,_ in assignment)
            degree = sum(n for _,n,_ in assignment)
            factor = prod(c for _,_,c in assignment)
            for b,m in enumerate(matrices):
                # J times an auxiliary traceless matrix is symmetric.
                # It therefore multiplies an antisymmetric colour path.
                entries[word,b,degree] += factor*m
                entries[word[::-1],b,degree] += factor*m*reverse*(1 if b==3 else -1)
        for (word,b,degree), c in entries.items():
            if c:
                relation.catalog[word].append((b,bit,degree,1<<k,c))
    return relation


def probes(power, scores, count):
    rng = Random(38119+power)
    probes = []
    initial = [tuple(map(tuple,w)) for w,_ in scores]
    rng.shuffle(initial)
    probes.extend(initial[:24])
    n = 2*power+2
    shapes = [s for s in partitions(n,2) if 3 <= len(s) and max(s) <= 6]
    samples = list(range(3))*2 + [3,4]*(power-2)
    while len(probes) < count:
        shape = rng.choice(shapes)
        rng.shuffle(samples)
        traces, start = [], 0
        for size in shape:
            traces.append(tuple((PHI+GHOSTS)[i] for i in samples[start:start+size]))
            start += size
        sign,w = so_canon(tuple(traces))
        if sign and encode(w) not in probes:
            probes.append(encode(w))
    return probes


def quotient_functional(case):
    """Subtract a functional of the commuting quotient from the ordered second pairing."""
    from flint import fmpq_mat
    dual = case['dual']
    power = dual['power']
    coords = fmpq_mat([r['cartan'] for r in case['rows']])
    values = fmpq_mat([[r['second_score']] for r in case['rows']])
    def pivots(m):
        red,rank = m.rref()
        return [next(j for j in range(red.ncols()) if red[i,j]) for i in range(rank)]
    rows,cols = pivots(coords.transpose()),pivots(coords)
    minor = fmpq_mat([[coords[i,j] for j in cols] for i in rows])
    functional = minor.solve(fmpq_mat([[values[i,0]] for i in rows]))
    assert all(sum(coords[i,j]*functional[k,0] for k,j in enumerate(cols))==values[i,0]
               for i in range(coords.nrows()))
    trace_to_sites = matrix(power)
    basis, ordered, scores = targets(power,dual)
    in_traces = [sum(trace_to_sites[i,j]*functional[k,0] for k,j in enumerate(cols))
                 for i in range(len(basis))]
    difference = Counter({encode(w):c for w,c in scores.items()})
    for w,(i,s) in ordered.items():
        if in_traces[i]:
            difference[encode(w)] -= s*Fraction(str(in_traces[i]))
    return {w:c for w,c in difference.items() if c}, dict(columns=cols, weights=[str(k[0]) for k in functional.tolist()])


def evaluate_relation(spec, power, colours, coordinates, functional):
    relation = determinant(power,spec)
    assert relation.pf_degree > colours
    values = [relation.coefficient(w) for w in coordinates]
    residual = sum(n*relation.coefficient(w) for w,n in functional.items())
    relation.clear()
    return [str(x) for x in values], str(residual)


def run_sample(i):
    spec = specification(i,POWER,COLOURS)
    row,residual = evaluate_relation(spec,POWER,COLOURS,PROBES,FUNCTIONAL)
    return i,spec,row,residual


def run(case, output, workers, attempts, probe_count):
    import json
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    from flint import fmpq_mat
    global POWER,COLOURS,PROBES,FUNCTIONAL
    dual = case['dual']
    POWER,COLOURS = dual['power'],dual['colours']
    PROBES = probes(POWER,dual['scores'],probe_count)
    FUNCTIONAL,summary = quotient_functional(case)
    dimension = relation_dimension(COLOURS,POWER,2)
    lower = [relation_dimension(COLOURS,POWER,k,False) for k in (0,1)]
    assert lower == [0,0]
    selected,failed = [],[]
    record = dict(colours=COLOURS,power=POWER,lower_full_weight_relations=lower,
                  grade_two_flavour_singlet_relations=dimension,site_functional=summary,
                  probes=PROBES,independent=selected,nonzero_residuals=failed)
    rows = []
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for start in range(0,attempts,workers):
            for i,spec,row,residual in pool.map(run_sample,range(start,min(attempts,start+workers))):
                if Fraction(residual):
                    failed.append(dict(index=i,residual=residual))
                if len(rows)<dimension and fmpq_mat(rows+[row]).rank()>len(rows):
                    rows.append(row)
                    selected.append(dict(index=i,spec=spec,probe_values=row,residual=residual))
            record.update(examined=min(attempts,start+workers),rank=len(rows))
            if rows:
                red,rank=fmpq_mat(rows).rref()
                columns=[next(j for j in range(red.ncols()) if red[k,j]) for k in range(rank)]
                record['minor_columns']=columns
                record['minor_determinant']=str(fmpq_mat([[r[j] for j in columns] for r in rows]).det())
            tmp=output.with_suffix('.part')
            tmp.write_text(json.dumps(record,separators=(',',':'))+'\n')
            tmp.replace(output)
            print(COLOURS,record['examined'],len(rows),'of',dimension,'nonzero residuals',len(failed),flush=True)
            if len(rows)==dimension:
                break
    return record


if __name__=='__main__':
    import argparse
    import gzip
    import json
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('second',type=Path)
    parser.add_argument('--colours',type=int,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--attempts',type=int,default=640)
    parser.add_argument('--probes',type=int,default=100)
    args=parser.parse_args()
    if not __debug__ or min(args.workers,args.attempts,args.probes)<1:
        parser.error('positive work bounds and assertions are required')
    case=json.loads(gzip.decompress(args.second.read_bytes()))['cases'][str(args.colours)]
    run(case,args.out,args.workers,args.attempts,args.probes)
