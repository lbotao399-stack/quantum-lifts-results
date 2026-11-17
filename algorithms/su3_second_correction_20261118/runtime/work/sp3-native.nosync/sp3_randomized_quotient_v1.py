"""Sample all trace partitions; saturate Q0 boundaries plus closed products.

Every sampled row is the reduction of a rational extremal projection. Rank
saturation is a rational exclusion; an unsaturated run is only a bound.
"""
import sp_environment
import sp_census as census
import hw_extremal as exact_projection
import native_bmn_projection_v1 as projection
import native_bmn_derivation_rank_v1 as derivation
from sp_projected_basis import evaluate as python_evaluate
from extend_rows_free_columns_v1 import extend_rows
from sp_graviton_singlets import product_values
from fortuity.bmn.words import contents_of_sector, graviton_generators
from protected import graph_expr, expr_sec, monomials
from hw_products import Raising
from projected_checkpoint import atomic, json_bytes
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from collections import Counter, defaultdict
from fractions import Fraction
from math import ceil, comb
from pathlib import Path
import argparse, hashlib, json, multiprocessing, os, pickle, random, time
import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
CONTEXT = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def partitions(total, minimum=2):
    if total == 0:
        yield ()
    for first in range(minimum, total+1):
        for rest in partitions(total-first, first):
            yield (first,)+rest


def setup(weights, prime):
    global CONTEXT
    census.configure(3, prime)
    frames, dimensions, choices = {}, {}, {}
    for spin in range(4):
        sector = (spin,)+tuple(weights)
        dims = {int(d):int(v) for d,v in census.base.exact_hw(sector,0).items() if v}
        dimensions[spin] = dims
        frames[spin] = {}
        for degree, dimension in sorted(dims.items()):
            k = degree+2 if degree else 0
            frames[spin][degree] = census.base.Frame(degree, ceil((dimension+16)/comb(k,degree)),
                2026111701+spin*100003+degree*7919, prime, k=k)
    for spin in (1,2,3):
        contents = list(contents_of_sector((spin,)+tuple(weights)))
        choices[spin] = []
        for content in contents:
            letters = tuple(i for i,n in enumerate(content) for _ in range(n))
            shapes = tuple(partitions(len(letters)))
            choices[spin].append((letters, shapes))
    CONTEXT = dict(weights=tuple(weights), prime=prime, frames=frames,
                   dimensions=dimensions, choices=choices)


def batch(task):
    spin, seed, count = task
    data = CONTEXT; prime = data['prime']; rng = random.Random(seed)
    source, recipes, shapes_used, contents_used = [], [], Counter(), Counter()
    attempts = 0; started = time.monotonic(); projection_seconds = 0
    while len(source) < count and attempts < count*80:
        index = (attempts+seed) % len(data['choices'][spin]); attempts += 1
        letters, shapes = data['choices'][spin][index]
        letters = list(letters); rng.shuffle(letters)
        shape = rng.choice(shapes); offset = 0; traces = []
        for length in shape:
            traces.append(tuple(letters[offset:offset+length])); offset += length
        sign, word = census.base.canon((), tuple(traces))
        if not sign:
            continue
        begin = time.monotonic()
        expr = projection.project({word:1}, (spin,)+data['weights'], prime)
        projection_seconds += time.monotonic()-begin
        if not expr:
            continue
        source.append(expr); recipes.append(word)
        shapes_used[str(shape)] += 1; contents_used[index] += 1
    begin = time.monotonic()
    blocks = [derivation.evaluate(source,frame) for frame in data['frames'][spin-1].values()]
    values = np.concatenate(blocks, axis=1)
    return dict(spin=spin, recipes=recipes, values=values, attempts=attempts,
                trace_partitions=dict(shapes_used), content_patterns=dict(contents_used),
                projected_terms=sum(map(len,source)), projection_seconds=projection_seconds,
                evaluation_seconds=time.monotonic()-begin, seconds=time.monotonic()-started)


def products(sector, output):
    prime = CONTEXT['prime']; generators = []
    for name, expr in graviton_generators(6,1):
        expr = graph_expr(expr)
        if not expr:
            continue
        assert not census.base.action_expr(expr)
        odd = {sum(census.base.count(word)[3:6]) % 2 for word in expr}
        assert len(odd) == 1
        generators.append(dict(name=name, expr=expr, sector=expr_sec(expr),
                               odd=odd.pop(), reflection=0, family='trace'))
    words = monomials(generators,sector,0); raising = Raising(generators,prime)
    images = [{(root,w):c for root in (1,2) for w,c in raising.product(recipe,root).items()}
              for recipe in words]
    kernel, rank = census.base.formal_closed_kernel(images,prime)
    for row in kernel:
        residual = defaultdict(int)
        for index in np.flatnonzero(row):
            for word,c in images[index].items():
                residual[word] += int(row[index])*c
        assert not any(c % prime for c in residual.values())
    coefficients = csr_matrix(kernel)
    assert max(np.diff(coefficients.indptr), default=0)*(prime-1)**2 < 2**63
    values = np.concatenate([np.asarray(coefficients @ product_values(generators,words,frame,degree)) % prime
        for degree,frame in CONTEXT['frames'][sector[0]].items()],axis=1)
    payload = dict(generators=generators,products=words,kernel=kernel,values=values)
    path = output/('products_%d.pkl'%sector[0]); atomic(path,pickle.dumps(payload,protocol=5))
    return values, dict(sector=sector, generators=len(generators), weight_products=len(words),
        highest_weight_product_rows=len(kernel), raising_rank=rank, sha256=sha(path))


def control(output):
    records = []; prime = CONTEXT['prime']
    for spin in (1,2,3):
        saved = batch((spin,2026111719+spin,2)); expected = []
        for word in saved['recipes']:
            exact = exact_projection.project({word:1},(spin,)+CONTEXT['weights'])
            reduced = {w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime
                       for w,c in exact.items() if c}
            assert reduced == projection.project({word:1},(spin,)+CONTEXT['weights'],prime)
            assert not exact_projection.action(exact,0,1)
            assert not exact_projection.action(exact,1,2)
            image = census.base.action_expr(exact)
            assert not census.base.action_expr(image)
            image = {w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime
                     for w,c in image.items() if c}
            blocks = []
            for degree,frame in CONTEXT['frames'][spin-1].items():
                part = {w:c for w,c in image.items() if sum(census.base.count(w)[3:6]) == degree}
                blocks.append(python_evaluate([part],frame)[0])
            expected.append(np.concatenate(blocks))
        expected = np.asarray(expected)
        assert np.array_equal(expected,saved['values'])
        records.append(dict(spin=spin, rows=len(expected), independent_coordinates=int(expected.size),
                            trace_partitions=saved['trace_partitions']))
    report = dict(status='passed_random_partition_rational_projection_and_Q0_controls',
                  weights=CONTEXT['weights'],prime=prime,records=records,driver_sha256=sha(__file__))
    atomic(output/'control.json',json_bytes(report)); print(json.dumps(report),flush=True)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--weights',default='5,4,4')
    parser.add_argument('--prime',type=int,default=1000003); parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=36); parser.add_argument('--seconds',type=int,default=1200)
    parser.add_argument('--control',action='store_true'); args = parser.parse_args()
    start = time.monotonic(); weights = tuple(map(int,args.weights.split(',')))
    args.output.mkdir(parents=True,exist_ok=True); setup(weights,args.prime)
    if args.control:
        return control(args.output)
    assert os.environ.get('SLURM_JOB_ID')
    assert not (args.output/'manifest.json').exists()
    families = {}; joint = {}; boundary = {}; selected = {}; completed = Counter(); profile = Counter()
    for spin in (1,2,3):
        columns = sum(frame.ncols for frame in CONTEXT['frames'][spin-1].values())
        boundary[spin] = (np.zeros((0,columns),dtype=np.int64),[])
        selected[spin] = []
        if spin > 1:
            values, families[spin-1] = products((spin-1,)+weights,args.output)
            rows,pivots,_ = extend_rows(np.zeros((0,columns),dtype=np.int64),[],values,args.prime)
            joint[spin] = (rows,pivots)
    metadata = dict(group='Sp(3)',weights=weights,prime=args.prime,dimensions=CONTEXT['dimensions'],
        frames={s:{d:dict(points=f.P,k=f.k,columns=f.ncols) for d,f in fs.items()} for s,fs in CONTEXT['frames'].items()},
        families=families,driver_sha256=sha(__file__),sampling='All partitions into trace lengths at least two; seeded independent shuffles',
        rational_bridge='Apply the unchanged rational extremal projector to each integer source and to integer lifts of product-kernel rows. Projector denominators are below the pinned prime; exact Q0 closure is preserved.')
    atomic(args.output/'manifest.json',json_bytes(metadata))
    serial = 0; pending = {}; partitions_used = Counter()
    def bounds():
        return {s:sum(CONTEXT['dimensions'][s].values())-len(joint[s+1][0])-len(boundary[s][0]) for s in (1,2)}
    def save(status):
        limits = bounds(); assert min(limits.values()) >= 0
        atomic(args.output/'checkpoint.pkl',pickle.dumps(dict(boundary=boundary,joint=joint,
            selected=selected,next_serial=serial),protocol=5))
        report = dict(status=status,group='Sp(3)',weights=weights,prime=args.prime,
            outgoing_rank_lower_bounds={s:len(v[0]) for s,v in boundary.items()},
            closed_boundary_plus_product_ranks={s-1:len(v[0]) for s,v in joint.items()},
            fortuitous_quotient_upper_bounds=limits,rational_excluded_sectors=[(s,)+weights for s,b in limits.items() if not b],
            completed_batches=dict(completed),profile=dict(profile),trace_partitions=dict(partitions_used),
            seconds=time.monotonic()-start,checkpoint_sha256=sha(args.output/'checkpoint.pkl'),
            manifest_sha256=sha(args.output/'manifest.json'),job_id=os.environ['SLURM_JOB_ID'],
            node=os.environ.get('SLURMD_NODENAME'),d1_increment=0,new_generator_increment=0,
            d2_increment=0,d3_increment=0,nonperturbative_increment=0,
            scope='Saturation excludes additional classical classes modulo the specified universally closed products. Positive bounds do not count classes. No quantum survival is inferred.')
        atomic(args.output/'result.json',json_bytes(report)); print(json.dumps(report),flush=True)
    pool = ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn'),
                               initializer=setup,initargs=(weights,args.prime))
    def submit(spin):
        nonlocal serial
        pending[pool.submit(batch,(spin,2026111723+7919*serial,12))] = serial; serial += 1
    failed = False
    try:
        for i in range(args.workers):
            submit(i%3+1)
        last = time.monotonic()
        while pending:
            done,_ = wait(pending,timeout=2,return_when=FIRST_COMPLETED)
            for future in done:
                task = pending.pop(future); result = future.result(); spin = result['spin']; completed[spin] += 1
                atomic(args.output/('batch_%06d.pkl'%task),pickle.dumps(result,protocol=5))
                begin = time.monotonic(); rows,pivots = boundary[spin]
                rows,pivots,indices = extend_rows(rows,pivots,result['values'],args.prime)
                boundary[spin] = rows,pivots
                for index in indices:
                    selected[spin].append(dict(word=result['recipes'][index],values=result['values'][index]))
                if spin > 1:
                    rows,pivots = joint[spin]
                    rows,pivots,_ = extend_rows(rows,pivots,result['values'][indices],args.prime)
                    joint[spin] = rows,pivots
                profile['elimination_seconds'] += time.monotonic()-begin
                for key in ('projection_seconds','evaluation_seconds','projected_terms'):
                    profile[key] += result[key]
                partitions_used.update(result['trace_partitions'])
                limits = bounds(); assert min(limits.values()) >= 0
                relevant = limits[1] if spin==1 else limits[2] if spin==3 else max(limits.values())
                if relevant and time.monotonic()-start < args.seconds:
                    submit(spin)
            if time.monotonic()-last > 30:
                save('running'); last = time.monotonic()
            if not max(bounds().values()) or time.monotonic()-start >= args.seconds:
                for future in pending: future.cancel()
                break
    except BaseException:
        failed = True
        for future in pending: future.cancel()
        raise
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
        save('failed_checkpointed' if failed else 'rational_quotient_saturated' if not max(bounds().values()) else 'bounded_rank_checkpointed')


if __name__ == '__main__':
    main()
