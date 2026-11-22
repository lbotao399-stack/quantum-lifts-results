"""Prepare and evaluate the corrected F4 endpoints on the four norm-dressed sources."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
from itertools import chain
import json
from pathlib import Path
import subprocess
from zipfile import ZipFile

from endpoint_metric import branches, read
from verify import ARCHIVE, HERE, contents
from ward_boundary import graph, residue


def prepare(z, folder, unfolded=False):
    colour = read(z,'collected-O0/color-basis-native.json')
    rows = branches(colour, read(z,'independent-audit/target-normalization.json'))
    separate = list(rows)
    orbit = {tuple(k):(tuple(f),Fraction(s)) for k,(f,s) in
             read(z,'exact-score/eta0-folded-exact.json')['orbit']}
    if not unfolded:
        shortened = []
        for metric, row in separate:
            collected = Counter()
            for k,c in row.items():
                target,scale = orbit.get(k,(k,1))
                collected[target] += scale*c
            shortened.append((metric,collected))
        separate=shortened
    coordinates = read(z,'cochains/inputs/exact-coordinate-rows.json')['rows']
    first = {tuple(k) for row in coordinates for k,_ in row}
    rest = set().union(*(row.keys() for _,row in separate))-first
    keys = sorted(first)+sorted(rest)
    for src in ('full-score/legacy.cpp', 'full-score/native.in', 'full-score/trace_seed_native.in',
                'exact-score/exact_anchor.cpp'):
        (folder/Path(src).name).write_bytes(contents(z,src))
    data = contents(z,'exact-score/eta0-all-exact.in').decode().splitlines()
    metric = data[:1+int(data[0])]
    (folder/'queries.in').write_text('\n'.join(metric+[str(len(keys))]+
                                     [' '.join(map(str,k)) for k in keys])+'\n')
    primes = (1000003, 1000033)
    ends = [dict(graph(z,f'graph-{p}.in',p)) for p in primes]
    assert all(set(e).issubset(keys) for e in ends)
    weights = [[str(row.get(k,0)) for _,row in separate]+[e.get(k,0) for e in ends] for k in keys]
    (folder/'weights.json').write_text(json.dumps(weights,separators=(',',':')))
    key_index = {k:i for i,k in enumerate(keys)}
    restrictions = [[[key_index[tuple(k)],c] for k,c in row] for row in coordinates]
    (folder/'coordinate-rows.json').write_text(json.dumps(restrictions,separators=(',',':')))
    banks = ['4']
    for source in read(z,'full-score/closed-w1-recipes.json'):
        actual = []
        for traces,c in source['Psi']:
            coefficient = 36*Fraction(c)
            assert coefficient.denominator == 1
            terms = list(chain.from_iterable([len(w),*w] for w in traces))
            actual.append(' '.join(map(str,[int(coefficient),len(traces),*terms])))
        banks.extend([str(len(actual)),*actual])
    (folder/'recipes.in').write_text('\n'.join(banks)+'\n')
    # Preserve the native field and sign operations. Evaluate the four input
    # coefficient vectors together while sharing every colour trace and norm.
    cpp=(folder/'exact_anchor.cpp').read_text()
    for old,new in [('array<Z,2>','array<Z,4>'), ('bank[2]','bank[4]'),
                    ('assert(n==2);','assert(n==4);'), ('int measure=1;','constexpr int measure=4;'),
                    ('Z v[2]','Z v[4]'), ('Z z[2]','Z z[4]'), ('Z t[2]','Z t[4]'),
                    ('cout<<"Q "<<i<<" "<<dec(a[0])<<endl;',
                     'cout<<"Q "<<i; for(const auto v:a) cout<<" "<<dec(v); cout<<endl;')]:
        assert old in cpp
        cpp=cpp.replace(old,new)
    (folder/'four_anchor.cpp').write_text(cpp)
    return len(keys)


def contract(folder, start, stop):
    output = subprocess.check_output([str((folder/'anchor').resolve()), 'native.in',
                                      'recipes.in','queries.in', str(start),str(stop),'0'],
                                     text=True, cwd=folder)
    values = []
    for index, line in enumerate(output.splitlines(), start):
        tag,j,*n = line.split()
        assert tag=='Q' and index == int(j) and len(n)==4
        values.append(list(map(int,n)))
    assert len(values)==stop-start
    return values


def result(folder, values):
    raw = json.loads((folder/'weights.json').read_bytes())
    weights = [[Fraction(v) for v in r[:2]]+r[2:] for r in raw]
    rows = json.loads((folder/'coordinate-rows.json').read_bytes())
    query = (folder/'queries.in').read_text().splitlines()
    offset = int(query[0])+1
    keys = [tuple(map(int, line.split())) for line in query[offset+1:]]
    samples = json.loads((HERE/'A8-5693880.json').read_text())['selected_colour_keys_and_integer_numerators']
    by_key = dict(zip(keys, values, strict=True))
    assert len(by_key) == len(weights) == int(query[offset])
    assert all(by_key[tuple(k)] == n for k,n in samples)
    primes = (1000003,1000033)
    results = []
    for i,n in enumerate(zip(*values,strict=True)):
        # The archive's positive raw anchor is 6 times this test input.
        denom=147456*6
        v0=sum(c*r[0] for c,r in zip(n,weights,strict=True))/denom
        v1=sum(c*r[1] for c,r in zip(n,weights,strict=True))/denom
        corrected=-4*v0+2*v1
        moduli=[]
        for j,p in enumerate(primes):
            k=residue(Fraction(sum(c*r[2+j] for c,r in zip(n,weights,strict=True)),denom),p)
            moduli.append(dict(prime=p,etaK=k,corrected_sum=(residue(corrected,p)+k)%p))
        coord=[sum(n[idx]*Fraction(c) for idx,c in row)/denom for row in rows]
        results.append(dict(source=i,branch_values=[str(v0),str(v1)],eta0=str(corrected),
                            reductions=moduli,coordinates=list(map(str,coord)),
                            raw_numerators_sha256=sha256(json.dumps(n,separators=(',',':')).encode()).hexdigest()))
    return dict(archive_sha256=ARCHIVE,raw_components=len(weights),sources=results)


def evaluate(folder, workers, compiler='g++'):
    subprocess.run([compiler,'-std=c++20','-O3','-ftrapv','four_anchor.cpp','-o','anchor'],
                   check=True,cwd=folder)
    count=len(json.loads((folder/'weights.json').read_bytes()))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        chunks=[pool.submit(contract,folder,start,min(start+40,count)) for start in range(0,count,40)]
        values=[row for part in chunks for row in part.result()]
    return result(folder, values)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('folder',type=Path)
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--unfolded',action='store_true')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--compiler',default='g++')
    args=parser.parse_args()
    if not __debug__ or args.workers<1 or sha256(args.archive.read_bytes()).hexdigest()!=ARCHIVE:
        parser.error('use the pinned archive, positive worker count and enabled Python assertions')
    args.folder.mkdir(parents=True,exist_ok=True)
    with ZipFile(args.archive) as z:
        count=prepare(z,args.folder,args.unfolded)
    if args.run:
        print(json.dumps(evaluate(args.folder,args.workers,args.compiler),indent=2))
    else:
        print(count,'raw colour components per source')
