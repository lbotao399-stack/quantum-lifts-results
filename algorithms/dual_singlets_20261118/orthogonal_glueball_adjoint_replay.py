"""Recompute the third coordinate by pairing the source with quantum-parent columns."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product

from mixed_incoming import so_canon
from fortuity.full.letters import Q_letter, is_ferm
from orthogonal_glueball_first_pair import GHOSTS, PHI, PfaffianCoefficient
from orthogonal_glueball_relation_span import quotient_functional
from orthogonal_glueball_second_sources import Relation
from orthogonal_glueball_third_pair import decode, encode, targets
from rank_four_lowered_replay import Q1d_multi, terminal


@lru_cache(None)
def reconnect(lengths, i, j):
    """All inverse orthogonal Fierz reconnections at the two marked vertices."""
    links = {}
    offset = 0
    for length in lengths:
        for k in range(length):
            left = 2*(offset+k)+1
            right = 2*(offset+(k+1)%length)
            links[left],links[right] = right,left
        offset += length
    paths = set()
    marked_edges = set(product((2*i,2*i+1),(2*j,2*j+1)))
    # An output edge can join both marked letters. Its other reconnected edge
    # can then lie anywhere, including on a separate trace with neither mark.
    for a in (2*i,2*i+1):
        if links[a]//2==j:
            marked_edges.update((a,c) for c,d in links.items() if c<d)
    for a,c in marked_edges:
        b,d = links[a],links[c]
        pairs = [((a,b),(c,d))]
        if len({a,b,c,d})==4:
            pairs.extend([((a,c),(b,d)),((a,d),(b,c))])
        for replacements in pairs:
            edges = links.copy()
            for u,v in replacements:
                edges[u],edges[v] = v,u
            unseen = set(range(offset))
            cycles = []
            while unseen:
                start = 2*min(unseen)
                pos = start
                cycle = []
                while True:
                    cycle.append(pos//2)
                    unseen.remove(pos//2)
                    pos = edges[pos^1]
                    if pos==start:
                        break
                cycles.append(tuple(cycle))
            if all(len(t)>1 for t in cycles):
                paths.add(tuple(cycles))
    return tuple(sorted(paths))


def parents(support, active_even=True):
    result = set()
    for word in support:
        letters = [v for t in word for v in t]
        lengths = tuple(map(len,word))
        for i,p in enumerate(letters):
            if p!=GHOSTS[0]:
                continue
            for j,m in enumerate(letters):
                if m!=GHOSTS[1]:
                    continue
                for layout in reconnect(lengths,i,j):
                    for u,v in ((i,j),(j,i)):
                        active = next(t for t in layout if u in t)
                        if active_even and len(active)%2:
                            continue
                        for k in range(3):
                            labelled = list(letters)
                            labelled[u] = (0,0,tuple(f for f in range(3) if f!=k))
                            labelled[v] = PHI[k]
                            sign,w = so_canon(tuple(tuple(labelled[z] for z in t) for t in layout))
                            if sign:
                                result.add(w)
    return sorted(result)


@lru_cache(50000)
def inverse_reversal(trace):
    """Enumerate the two positions in an odd word merged by the reversal primitive."""
    out = Counter()
    visited = set()
    for oriented in {v[k:]+v[:k] for v in (trace,trace[::-1]) for k in range(len(trace))}:
        pos = next(k for k,v in enumerate(oriented) if len(v[2])==2)
        jet = oriented[pos]
        for coeff,x,y in Q_letter(jet):
            assert x in PHI and y in PHI and abs(coeff)==1
            middle = oriented[:pos]
            tail = oriented[pos+1:]
            for last in range(len(tail)+1):
                end, prefix = tail[:last],tail[last:][::-1]
                original = prefix + (x,) + middle + (y,) + end
                s,w = so_canon((original,))
                if not s or w!=(original,):
                    continue
                i,j = len(prefix),len(prefix)+1+len(middle)
                key = original,i,j
                if key in visited:
                    continue
                visited.add(key)
                # The current bubble ordering is middle,x,y,end,reversed(prefix).
                indices = [*range(i+1,j), i,j,*range(j+1,len(original)),*range(i-1,-1,-1)]
                odd = [z for z in indices if is_ferm(original[z])]
                koszul = sum(z>t for k,z in enumerate(odd) for t in odd[k+1:])
                koszul += sum(map(is_ferm,middle))
                # Parent and spectator traces will be canonicalized together below.
                out[original,oriented] += Fraction(coeff*(-1)**koszul,2)
    return tuple((w,t,c) for (w,t),c in out.items() if c)


def leading_terms(parent):
    active = next(k for k,t in enumerate(parent) if any(len(l[2])==2 for l in t))
    if len(parent[active])%2:
        return {}
    spectators = parent[:active]+parent[active+1:]
    out = Counter()
    for original,merged,c in inverse_reversal(parent[active]):
        s,source = so_canon(spectators+(original,))
        if not s:
            continue
        i = next(k for k,t in enumerate(source) if len(t)%2)
        if source[i]!=original:
            continue
        pre = sum(is_ferm(l) for tr in source[:i] for l in tr)
        sign,term = so_canon(source[:i]+(merged,)+source[i+1:])
        if sign:
            assert term==parent
            out[source] += (-1)**pre*sign*c
    return {w:c for w,c in out.items() if c}


def source_relation():
    rel = Relation(6,dict(pattern=(),indices=(),plus=(0,1,2),minus=(0,1,2)))
    rel.project = PfaffianCoefficient(3,0)
    # Six single ghost markers and nine single scalar insertions in a 22-by-22 Pfaffian.
    rel.pf_degree = 15
    return rel


def interval(bounds):
    from hashlib import sha256
    import json
    start,end = bounds
    relation = source_relation()
    score = Fraction()
    hits = 0
    digest = sha256()
    for v in PARENTS[start:end]:
        adjoint = Fraction()
        for raw,c in Q1d_multi(v,11,'SO_odd').items():
            s,w = so_canon(raw)
            if s:
                adjoint -= s*c*DUAL.get(w,0)
        if adjoint:
            source = sum(c*relation.coefficient(encode(w)) for w,c in leading_terms(v).items())
            if source:
                hits += 1
                score += adjoint*source
                digest.update((json.dumps([v,str(adjoint),str(source)],separators=(',',':'))+'\n').encode())
    relation.clear()
    return [start,end,hits,str(score),digest.hexdigest()]


def run(case,workers):
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    global PARENTS,DUAL
    raw,_ = quotient_functional(case)
    DUAL = {decode(w):n for w,n in raw.items()}
    PARENTS = parents(DUAL)
    print('dual support',len(DUAL),'quantum parents',len(PARENTS),flush=True)
    splits = [(i,min(len(PARENTS),i+80)) for i in range(0,len(PARENTS),80)]
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        blocks = list(pool.map(interval,splits))
    total = sum(Fraction(row[3]) for row in blocks)
    return dict(parents=len(PARENTS),dual_support=len(DUAL),nonzero=sum(r[2] for r in blocks),
                coordinate=str(total),blocks=blocks)


def first_interval(bounds):
    start,end=bounds
    relation=source_relation()
    total=Counter()
    active=0
    for v in PARENTS[start:end]:
        projected=Counter()
        for raw,c in Q1d_multi(v,11,'SO_odd').items():
            s,w=terminal(raw)
            if s:
                projected[POSITIONS[w]]+=s*c
        if any(projected.values()):
            coefficient=sum(c*relation.coefficient(encode(w)) for w,c in leading_terms(v).items())
            if coefficient:
                active+=1
                for i,n in projected.items():
                    total[i]+=n*coefficient
    relation.clear()
    return start,end,active,{i:str(v) for i,v in total.items() if v}


def first_action(case,workers):
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    global PARENTS,POSITIONS
    basis,ordered,_=targets(6,case['dual'])
    POSITIONS={term:i for i,term in enumerate(basis)}
    PARENTS=parents(ordered)
    print('first-action parents',len(PARENTS),flush=True)
    bounds=[(i,min(i+80,len(PARENTS))) for i in range(0,len(PARENTS),80)]
    total=Counter()
    active=0
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for _,_,number,values in pool.map(first_interval,bounds):
            active+=number
            for i,c in values.items():
                total[i]+=Fraction(c)
    return dict(parents=len(PARENTS),nonzero=active,first_trace=[str(total[i]) for i in range(len(basis))])


if __name__=='__main__':
    import argparse
    import gzip
    import json
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('second',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--first',action='store_true')
    args=parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('run with assertions enabled and positive worker count')
    case=json.loads(gzip.decompress(args.second.read_bytes()))['cases']['11']
    output=(first_action if args.first else run)(case,args.workers)
    args.out.write_text(json.dumps(output,separators=(',',':'))+'\n')
    print(output.get('coordinate',output['nonzero']),flush=True)
