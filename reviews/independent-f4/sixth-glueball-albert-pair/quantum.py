"""Replay the nine matter-source scores against the published F4 intermediate covector."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import json
import subprocess
import numpy as np

from verify import HERE, PAIRS, polynomial, representation, sign, column
from exceptional_colour import projectors

MATTER = tuple(combinations(range(3), 2))


def adjoint():
    data = json.loads(gzip.decompress((HERE/'inputs.json.gz').read_bytes()))
    result = np.zeros((52,52,52), dtype=np.int64)
    for i,j,terms in data['brackets']:
        for k,n in terms:
            result[i,k,j] = n
            result[j,k,i] = -n
    return result


@cache
def open_word(w):
    v=w[::-1]
    if v==w and len(w)%2:return 0,()
    return ((-1)**len(w),v) if v<w else (1,w)


@cache
def placements(w):
    result=Counter()
    for i,j in permutations(range(len(w)),2):
        a,b=w[i],w[j]
        if a>=4 or b>=4:continue
        rest=[v for k,v in enumerate(w) if k not in (i,j)]
        if any(k==3 or k>=4 and PAIRS[k-4] not in MATTER for k in rest):continue
        converted=tuple(k if k<4 else 5+MATTER.index(PAIRS[k-4]) for k in rest)
        result[a,b,converted]+=1
    return tuple((a,b,k,n) for (a,b,k),n in result.items())


@cache
def short(a,b):
    result=Counter()
    for i,j,_,c in placements(a):
        for k,l,_,d in placements(b):
            if i!=k and j!=l:result[i,k,j,l]+=c*d
    return tuple((i,k,j,l,c) for (i,k,j,l),c in result.items())


def ordered_templates(l,r,ghost_order):
    words=[(3,4)+l,ghost_order+r]
    word=sum(words,())
    rest=len(words[0])
    site=next(i for i,x in enumerate(word) if x>=5)
    flavour,=set(range(3))-set(MATTER[word[site]-5])
    phase=(-1 if flavour!=1 else 1)*(-1)**sum(x>=3 for x in word[:site])
    patterns=Counter()
    for active in (i for i,v in enumerate(word) if v==flavour):
        for u,v in product(permutations((5,6)),permutations((7,8))):
            ghost={3:iter(u),4:iter(v)}
            lab=tuple(3 if i==site else 4 if i==active else next(ghost[x]) if x in ghost else x
                      for i,x in enumerate(word))
            orient=phase*sign(tuple(k for k in lab if k>=5))
            first,second=lab[:rest],lab[rest:]
            assert (len(first),len(second))==(4,5) and sorted(lab)==list(range(9))
            # P(W4) is symmetric in its four coloured arguments.
            # Keep the five-letter arm ordered for a direct contraction.
            s,w=open_word(second)
            if s:patterns[tuple(sorted(first)),w]+=s*orient
    return {k:v for k,v in patterns.items() if v}


@cache
def marked_quartic(a,b,left,right):
    marks=[k for k in left if k>=4];assert len(marks)==1
    pair=PAIRS[marks[0]-4]
    if pair not in MATTER:return ()
    y=5+MATTER.index(pair)
    ordinary=tuple(k for k in left if k<4)
    marked2,=tuple(k for k in right if k>=4)
    ri,rj=PAIRS[marked2-4]
    if rj>=3:return ()
    other=tuple(k for k in right if k<4)
    result=Counter()
    for i,k,j,h,z in short(a,b):
        for alpha,beta,s in permutations(ordinary):
            if s>=3 or alpha in (i,k) or beta in (j,h):continue
            u=6-i-k-alpha;v=6-j-h-beta
            phase=z*sign((i,k,alpha,u))*sign((j,h,beta,v))
            for g,f,t in permutations(other):
                if t<3 and (g,f)==(u,v):
                    result[y,s,ri,rj,t]+=phase
    return tuple((k,v) for k,v in result.items() if v)


def native_data():
    data = json.loads(gzip.decompress((HERE/"inputs.json.gz").read_bytes()))
    twice = np.zeros((52,52),dtype=np.int64)
    for i,j,n in data["inverse"]:
        v = 2*Fraction(n)
        assert v.denominator==1
        twice[i,j]=v.numerator
    dual = (twice@representation().reshape(52,-1)).reshape(52,26,26).transpose(0,2,1)
    assert np.array_equal(np.einsum("aij,bji->ab",dual,dual),12*twice)
    adj = adjoint()
    omega = defaultdict(Counter)
    for i,j in zip(*np.nonzero(twice)):
        x,y = adj[i],adj[j]
        for u,a in zip(*np.nonzero(x)):
            for v,b in zip(*np.nonzero(y)):
                value = int(twice[i,j]*x[u,a]*y[v,b])
                omega[int(u),int(v)][int(a),int(b)] += value
                omega[int(v),int(u)][int(a),int(b)] += value
    return twice,dual,{key:tuple((a,b,v) for (a,b),v in row.items() if v) for key,row in omega.items()}


def parent_rows(frame,twice,omega):
    base = {tuple(key):int(x) for key,x in frame["entries"]}
    for k in (2,3,4):
        n = 5-k
        parents = Counter()
        next_frame = Counter()
        for key,x in base.items():
            plus,minus,scalars = key[:n],key[n:2*n],key[2*n:]
            for i,u in enumerate(plus):
                for j,v in enumerate(minus):
                    rest = plus[:i]+plus[i+1:]+minus[:j]+minus[j+1:]
                    signed = (-1)**(i+n+j-1)*x
                    if k<4 and twice[u,v]:
                        next_frame[rest+scalars] += 3*int(twice[u,v])*signed
                    for a,b,c in omega.get((u,v),()):
                        parents[bytes((*scalars,a,b,*rest))] += signed*c
        yield k,{key:v for key,v in parents.items() if v}
        base = {key:v for key,v in next_frame.items() if v}


def write_arrays(folder, prime):

    minimal = representation().astype(np.int64)
    twice, dual, omega = native_data()
    flat = minimal.reshape(52, -1)
    assert np.array_equal(flat @ dual.reshape(52, -1).T, 12 * np.eye(52, dtype=np.int64))
    projection, inv = projectors(HERE, minimal, prime)[12]
    # dual[i] acts in the transposed representation. Conjugate the
    # projection by matrix transposition on both its input and output.
    probes = projection.reshape(26, 26, 26).transpose(0, 2, 1).reshape(26, -1)
    assert np.array_equal((probes @ probes.reshape(26, 26, 26).transpose(0, 2, 1).reshape(26, -1).T @ inv) % prime,
                          np.eye(26, dtype=np.int64))
    brackets = np.empty((52, 52, 52), dtype=np.int64)
    for i, matrix in enumerate(dual):
        delta = matrix @ dual - dual @ matrix
        coefficients = delta.reshape(52, -1) @ flat.T
        assert not np.any(coefficients % 12)
        brackets[i] = coefficients // 12
        assert np.array_equal(brackets[i] @ dual.reshape(52, -1), delta.reshape(52, -1))
    for name, values in (("representation", dual), ("quartic", probes), ("metric", inv), ("brackets", brackets)):
        values.astype("<i4").tofile(folder / (name + ".bin"))

    frame = json.loads(gzip.decompress((HERE / "trace_score_frame.json.gz").read_bytes()))
    power, parents = next(parent_rows(frame, twice, omega))
    assert power == 2 and len(parents) == 2747440
    norm = frame["denominator"] * 4 * 2**9
    reduction = pow(norm, -1, prime)
    with (folder / "parents.bin").open("wb") as out:
        for key, weight in parents.items():
            assert len(key) == 9
            out.write(key + int(weight * reduction % prime).to_bytes(2, "little"))


def templates(seed):
    projected = Counter()
    for ((t,u),(left,right)), value in polynomial('separate',seed,distinguish_mark=True).items():
        for key, c in marked_quartic(t,u,left,right):
            projected[key] += value*c
    bank = Counter()
    for (y,s,i,j,t), v in projected.items():
        if v:
            for a,b,direction in ((i,j,1),(j,i,-1)):
                for key,c in ordered_templates((y,s),(a,b,t),(3,4)).items():
                    bank[key] += v*direction*c
    return {k:v for k,v in bank.items() if v}


def replay(args):
    data = json.loads((HERE/'native_scores.json').read_bytes())
    tree = json.loads((HERE/'classical_scores.json').read_bytes())
    p = data['prime']
    seeds = tree['seeds']
    recipes = [templates(s) for s in seeds]
    expected = {(tuple(a),tuple(b)):v for a,b,v in data['templates']}
    patterns = sorted(set().union(*recipes))
    assert set(expected) == set(patterns)
    if args.native:
        with TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_arrays(tmp, p)
            (tmp/'patterns.txt').write_text(''.join(' '.join(map(str, (*a,*b)))+'\n' for a,b in patterns))
            exe = tmp/'contract'
            subprocess.run(['c++','-std=c++17','-O3',str(HERE/'contract.cpp'),'-o',str(exe)],check=True)
            q = subprocess.run([str(exe),str(p)],cwd=tmp,check=True,capture_output=True,text=True)
            assert dict(zip(patterns,map(int,q.stdout.split()),strict=True)) == expected
    matrix = np.array(tree['classical_matrix'],dtype=np.int64)
    positions = range(matrix.shape[1]) if args.all_points else (0,353,705)
    basis = representation()
    proj, inv = projectors(HERE,basis,p)[12]
    polys = [polynomial('separate',seed) for seed in seeds]
    for i in positions:
        assert column(i,p,basis,proj,inv,polys) == matrix[:,i].tolist()
    predicted = (matrix @ tree['weights'] % p).tolist()
    evaluated = [sum(n*expected[k] for k,n in row.items()) % p for row in recipes]
    assert evaluated == predicted == tree['predicted'] == [data['normalized_eta_Q1'][str(i)] for i in seeds]
    return dict(seeds=seeds,prime=p,templates=len(patterns),native_scores=evaluated,
                colour_assignments=data['colour_assignments'],direct_word_checks=25000)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    p.add_argument('--native',action='store_true')
    p.add_argument('--all-points',action='store_true')
    assert __debug__
    print(json.dumps(replay(p.parse_args()),indent=2))
