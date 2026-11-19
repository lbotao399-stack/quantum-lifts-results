"""Four-colour-slot antisymmetrizer for three colours, with homogeneous Grassmann coefficients.

Canonical trace words and parities follow the archived full-letter runtime.
A tuple of four nonempty strings specifies four matrix products, in input order.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, json, random
import replay_formal_blocks as src
from fortuity.full.letters import is_ferm
from fortuity.full.words import canon_word, canon_multi

HERE=Path(__file__).resolve().parent
ODD=tuple(map(is_ferm,src.ALPH))
PERMS=tuple(permutations(range(4)))

def sign(seq): return -1 if sum(i>j for i,j in combinations(seq,2))%2 else 1

def cycles(perm):
    used=set(); out=[]
    for i in range(len(perm)):
        if i in used: continue
        t=[]; j=i
        while j not in used:
            used.add(j);t.append(j);j=perm[j]
        out.append(tuple(t))
    return out

@lru_cache(120000)
def canonical(w):
    s,t=canon_word(tuple(src.ALPH[i] for i in w))
    return (s,bytes(src.LETTER[x] for x in t)) if s else (0,None)

@lru_cache(20000)
def parity(w): return sum(ODD[i] for i in w)%2

def multi(traces):
    traces=list(traces); odd=[t for t in traces if parity(t)]
    if len(set(odd))!=len(odd):return 0,None
    return sign(odd),tuple(sorted(traces))

# Permutation acts on colour indices, which are even. Regrouping the four
# original entries in cycle order permutes the odd entries of their coefficients.
TABLE=[(sign(p),cycles(p)) for p in PERMS]

def identity(pieces):
    assert len(pieces)==4 and all(pieces)
    out=Counter(); odd=tuple(map(parity,pieces))
    for s,cyc in TABLE:
        s*=sign([i for t in cyc for i in t if odd[i]])
        traces=[]
        for t in cyc:
            r,w=canonical(b''.join(pieces[i] for i in t))
            if not r:break
            traces.append(w);s*=r
        else:
            e,w=multi(traces)
            if e:out[w]+=s*e
    return {k:v for k,v in out.items() if v}

def cuts(t,k):
    for choose in combinations(range(len(t)),k):
        rotated=t[choose[0]:]+t[:choose[0]]; cuts=[j-choose[0] for j in choose]+[len(t)]
        yield tuple(rotated[i:j] for i,j in zip(cuts,cuts[1:]))

def choices(w):
    seen=set()
    for sizes in product(*(range(1,min(len(t),5-len(w))+1) for t in w)):
        if sum(sizes)!=4: continue
        for blocks in product(*(cuts(t,n) for t,n in zip(w,sizes))):
            pieces=tuple(x for t in blocks for x in t)
            if pieces not in seen:
                seen.add(pieces);yield pieces

def compact_witnesses(rows):
    # Prefer a relation annihilating precisely one of the selected coordinates.
    # The output stores only rows with that property; failure is explicit.
    all_words={src.unpack_hex(r['outgoing']) for r in rows}; chosen=[];missing=[]
    for row in rows:
        w=src.unpack_hex(row['outgoing'])
        for pieces in choices(w):
            z=identity(pieces)
            if w in z and all_words.intersection(z)=={w}:
                chosen.append(dict(outgoing=src.hexes(w),pieces=src.hexes(pieces),pivot=z[w],
                                   multiplier=str(-Fraction(row['coefficient'])/z[w]),terms=len(z),
                                   trace_degree_blocks=len(set(map(src.block,z)))))
                break
        else:missing.append(src.hexes(w))
    return chosen,missing

# Independent evaluation in an ordinary integral matrix algebra tensored with an
# exterior algebra. No cyclic, multitace, or permutation signs from identity().
def mm(A,B):return tuple(tuple(sum(a*b for a,b in zip(row,col)) for col in zip(*B)) for row in A)

def traceless(rng,N):
    M=[[rng.randrange(-4,5) for _ in range(N)] for _ in range(N)];M[-1][-1]=-sum(M[i][i] for i in range(N-1));return tuple(map(tuple,M))

def assignment(w,N,seed):
    rng=random.Random(seed);letters=sorted(set(b''.join(w)));n=sum(ODD[i] for t in w for i in t)
    return {(i,k):traceless(rng,N) for i in letters for k in (range(n) if ODD[i] else [-1])},n

def exterior_value(w,M,n):
    # Each odd matrix is sum_k eta_k M_k, in increasing exterior basis order.
    result=0
    for perm in permutations(range(n)):
        it=iter(perm);val=sign(perm)
        for tr in w:
            matrices=[M[i,next(it) if ODD[i] else -1] for i in tr];a=matrices[0]
            for b in matrices[1:]:a=mm(a,b)
            val*=sum(r[j] for j,r in enumerate(a))
        result+=val
    return result

def controls(w,relation):
    report=[]
    for N in (3,4):
        for seed in (714293,713854):
            M,n=assignment(w,N,seed); term=[c*exterior_value(t,M,n) for t,c in sorted(relation.items())]
            residual=sum(term); pivot=exterior_value(w,M,n)
            assert pivot!=0 and (residual==0)==(N==3)
            report.append(dict(matrix_size=N,seed=seed,grassmann_generators=n,witness_value=pivot,
                               nonzero_summands=sum(v!=0 for v in term),summands=[str(v) for v in term],residual=str(residual)))
    return report

def parity_controls():
    rng=random.Random(341120); matrices=[traceless(rng,3) for _ in range(4)];proof=[]
    for odd in product((0,1),repeat=4):
        total=0;incorrect=0
        for p in PERMS:
            cyc=cycles(p); order=[i for t in cyc for i in t if odd[i]]
            reordered=sign(order);v=sign(p)
            for t in cyc:
                a=matrices[t[0]]
                for i in t[1:]:a=mm(a,matrices[i])
                v*=sum(a[i][i] for i in range(3))
            # Directly multiply exterior basis generators in the cycle traversal.
            seen=0;mv=1
            for t in cyc:
                for i in t:
                    if odd[i]:mv*=(-1)**((seen>>(i+1)).bit_count());seen|=1<<i
            total+=v*mv*reordered;incorrect+=v*mv
        assert total==0;proof.append([list(odd),total,incorrect])
    return proof

def run():
    rows=src.read(HERE/'finite_colour/replayed_733.json.gz')['coefficients'];chosen,missing=compact_witnesses(rows)
    w=src.unpack_hex(chosen[0]['outgoing']);z=identity(src.unpack_hex(chosen[0]['pieces'])); checks=controls(w,z)
    return dict(input='finite_colour/replayed_733.json.gz',isolated_coordinates=len(chosen),not_isolated=missing,
                choices=chosen,first_relation=[[src.hexes(t),c] for t,c in sorted(z.items())],
                first_witness_explicit_exterior_matrices=checks,parity_controls=parity_controls())

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__
    r=run();src.write(a.out,r);print('isolated',r['isolated_coordinates'],'missing',len(r['not_isolated']),'first',r['choices'][0]);print([(a['matrix_size'],a['witness_value'],a['residual']) for a in r['first_witness_explicit_exterior_matrices']])
