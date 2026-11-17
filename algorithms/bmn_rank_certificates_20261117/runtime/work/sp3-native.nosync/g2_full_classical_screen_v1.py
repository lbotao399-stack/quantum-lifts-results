"""Targeted G2 full-letter cohomology bounds with exact character dimensions.

Only the classical differential is computed. Evaluation ranks give lower
bounds over Q; a positive cohomology remainder is not a certified class.
The G2 matrices are constructed as the stabilizer of an explicit 3-form.
No gauge-group trace completeness assumption is used: each basis must
saturate the independently computed Weyl-character dimension.
"""
import sp_environment
import full_extremal as projection
import full_hw_so as full
from fortuity.full.evalg import Frame, nferm
from fortuity.full.letters import is_ferm
from fortuity.full.wordenum import multitraces_with_content
from sp_full_projected_basis import evaluate
from sp_projected_basis import independent_rows
from screen_gaugino_source_dimensions_v1 import content_patterns
from lie_root_dimension_v1 import multiplicity_parts
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import ceil, comb, lcm
from pathlib import Path
from flint import fmpq_mat, nmod_mat
import argparse
import hashlib
import json
import numpy as np
import pickle
import random
import signal
import time


HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


@lru_cache(None)
def dimensions(mu, length, group):
    total = defaultdict(int)
    rho = (2, 1, 0)
    for spin_shift, spin_sign in ((0, 1), (1, -1)):
        for permutation in permutations(range(3)):
            sign = (-1)**sum(permutation[i] > permutation[k]
                             for i in range(3) for k in range(i+1, 3))
            shifted = (mu[0]+spin_shift, mu[1]-spin_shift)+tuple(
                mu[i+2]+rho[i]-rho[permutation[i]] for i in range(3))
            if min(shifted) < 0:
                continue
            for bosons, fermions, count in content_patterns(shifted, length):
                total[sum(fermions)] += spin_sign*sign*count*multiplicity_parts(bosons, fermions, group)
    assert all(value >= 0 for value in total.values()), (mu, length, group, total)
    return {degree: value for degree, value in sorted(total.items()) if value}


def sector(mu, length):
    theta = sum(mu[2:])
    return (2*mu[0]+theta-length, 2*mu[1]+theta-length)+tuple(mu[0]+mu[1]+2*x for x in mu[2:])


@lru_cache(None)
def contents(mu, length):
    if min(mu) < 0 or length < 0 or sum(mu) < length:
        return ()
    maximum = sum(mu)-length+1
    alphabet = []
    for i, j in product(range(mu[0]+1), range(mu[1]+1)):
        for bits in product((0, 1), repeat=3):
            weight = (i, j)+bits
            if 1 <= sum(weight) <= maximum and all(weight[k] <= mu[k] for k in range(5)):
                alphabet.append((weight, (i, j, tuple(k for k, bit in enumerate(bits) if bit))))
    alphabet.sort(reverse=True)
    answer = []
    def visit(start, remaining, left, word):
        if not left:
            if not any(remaining):
                answer.append(tuple(sorted(word)))
            return
        for index in range(start, len(alphabet)):
            weight, letter = alphabet[index]
            count_max = min([left]+[remaining[k]//weight[k] for k in range(5) if weight[k]])
            for count in range(1, count_max+1):
                visit(index+1, tuple(remaining[k]-count*weight[k] for k in range(5)),
                      left-count, word+(letter,)*count)
    visit(0, mu, length, ())
    return tuple(answer)


@lru_cache(None)
def g2_basis():
    phi = np.zeros((7, 7, 7), dtype=np.int64)
    for triple, coefficient in (((0,1,2),1), ((0,3,4),1), ((0,5,6),1),
                                ((1,3,5),1), ((1,4,6),-1), ((2,3,6),-1), ((2,4,5),-1)):
        for order in permutations(range(3)):
            sign = (-1)**sum(order[i] > order[j] for i in range(3) for j in range(i+1,3))
            phi[tuple(triple[i] for i in order)] = coefficient*sign
    antisymmetric = []
    for i, j in combinations(range(7), 2):
        matrix = np.zeros((7,7), dtype=np.int64)
        matrix[i,j], matrix[j,i] = 1, -1
        antisymmetric.append(matrix)
    equations = []
    for i,j,k in combinations(range(7),3):
        equations.append([int(sum(matrix[t,i]*phi[t,j,k]+matrix[t,j]*phi[i,t,k]+
                                  matrix[t,k]*phi[i,j,t] for t in range(7))) for matrix in antisymmetric])
    reduced, rank = fmpq_mat(equations).rref()
    assert rank == 7
    pivots = [next(j for j in range(21) if reduced[i,j]) for i in range(rank)]
    result = []
    for j in range(21):
        if j in pivots:
            continue
        row = [Fraction(0)]*21
        row[j] = Fraction(1)
        for i, pivot in enumerate(pivots):
            row[pivot] = -Fraction(str(reduced[i,j]))
        denominator = lcm(*(x.denominator for x in row))
        result.append(sum((int(x*denominator)*matrix for x,matrix in zip(row,antisymmetric)),
                          np.zeros((7,7),dtype=np.int64)))
    basis = np.array(result, dtype=np.int64)
    assert len(basis) == 14
    flat = fmpq_mat(basis.reshape(14,49).tolist())
    assert flat.rank() == 14
    checks = 0
    for matrix in basis:
        assert not (matrix+matrix.T).any()
        assert not (np.einsum('ti,tjk->ijk',matrix,phi)+
                    np.einsum('tj,itk->ijk',matrix,phi)+
                    np.einsum('tk,ijt->ijk',matrix,phi)).any()
    for x,y in combinations(basis,2):
        bracket = x@y-y@x
        assert fmpq_mat(np.concatenate((basis.reshape(14,49),bracket.reshape(1,49))).tolist()).rank() == 14
        checks += 1
    return basis, dict(stabilizer_constraint_rank=rank, dimension=14,
                      exact_bracket_closure_checks=checks, exact_three_form_checks=14,
                      basis_sha256=hashlib.sha256(basis.tobytes()).hexdigest())


class GaugeFrame(Frame):
    def __init__(self, group, degree, dimension, letters, prime, seed):
        k = degree+2 if degree else 0
        points = ceil((dimension+24)/comb(k,degree))
        if group == 'G2':
            super().__init__(7,points,k,[degree],letters,p=prime,seed=seed,group='SU')
            basis, _ = g2_basis()
            rng = np.random.default_rng(seed)
            for letter in sorted(letters):
                odd = is_ferm(letter)
                coefficients = rng.integers(0,prime,(points,k if odd else 1,14),dtype=np.int64)
                self.L[letter] = (int(odd),np.einsum('...t,tij->...ij',coefficients,basis)%prime)
        elif group == 'Sp2':
            super().__init__(2,points,k,[degree],letters,p=prime,seed=seed,group='Sp')
        else:
            raise NotImplementedError(group)


def reduce_expr(expr, prime):
    return {word: Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime
            for word,c in expr.items() if c}


def build(mu, length, group, prime, output):
    output.mkdir(parents=True,exist_ok=True)
    if (output/'result.json').exists():
        result=json.loads((output/'result.json').read_text())
        assert result['basis_sha256']==digest(output/'basis.pkl')
        return pickle.loads((output/'basis.pkl').read_bytes()),result
    start = time.monotonic()
    dims = dimensions(mu,length,group)
    words = set()
    for content in contents(mu,length):
        for word in multitraces_with_content(content,length):
            sign, canonical = full.so_canon(word)
            if sign:
                words.add(canonical)
    words = sorted(words)
    random.Random(393161).shuffle(words)
    print(json.dumps(dict(stage='enumerated',group=group,mu=mu,length=length,
                          words=len(words),dimensions=dims)),flush=True)
    frames,values,chosen,candidates = {},{},{},{}
    for degree,dimension in dims.items():
        letters=sorted({letter for word in words if nferm(word)==degree for trace in word for letter in trace})
        frames[degree]=GaugeFrame(group,degree,dimension,letters,prime,991031+7919*degree)
        values[degree]=np.zeros((0,frames[degree].ncols),dtype=np.int64)
        chosen[degree],candidates[degree]=[],[]
    def flush(degree):
        if not candidates[degree]:return
        new=evaluate([expr for _,expr in candidates[degree]],frames[degree])
        combined=np.concatenate((values[degree],new))
        indices=independent_rows(combined,prime)
        combined_sources=chosen[degree]+candidates[degree]
        chosen[degree]=[combined_sources[i] for i in indices]
        values[degree]=combined[indices]
        candidates[degree]=[]
        assert len(indices)<=dims[degree]
    examined=0
    for word in words:
        degree=nferm(word)
        if degree not in dims or len(chosen[degree])==dims[degree]:continue
        examined+=1
        expr=projection.project({word:1},sector(mu,length),prime)
        if not expr:continue
        for raising in (1,2,3):assert not full.raise_expr(expr,raising,prime)
        candidates[degree].append((word,expr))
        if len(candidates[degree])>=8:flush(degree)
        if all(len(chosen[d])==dimension for d,dimension in dims.items()):break
    basis=[]
    for degree in dims:
        flush(degree)
        assert len(chosen[degree])==dims[degree],('trace_span_incomplete',group,mu,length,degree,len(chosen[degree]),dims[degree])
        for word,expr in chosen[degree]:
            rational=projection.project({word:1},sector(mu,length))
            assert reduce_expr(rational,prime)==expr
            for row,column in ((0,1),(1,2),('spin',1)):
                assert not projection.action(rational,row,column)
            basis.append(rational)
    (output/'basis.pkl').write_bytes(pickle.dumps(basis,protocol=4))
    result=dict(group=group,mu=mu,sector=sector(mu,length),length=length,prime=prime,
                dimensions=dims,dimension=len(basis),canonical_words=len(words),examined=examined,
                rational_highest_weight_checks=3*len(basis),basis_sha256=digest(output/'basis.pkl'),
                seconds=time.monotonic()-start)
    save(output/'result.json',result)
    print(json.dumps(dict(stage='basis_complete',**result)),flush=True)
    return basis,result


def map_ranks(mu,length,group,prime,output):
    start=time.monotonic()
    bases={};records={}
    for offset in (-1,0,1):
        bases[offset],records[offset]=build(mu,length+offset,group,prime,output/('length%d'%(length+offset)))
    images={offset:[sp_environment.canonical(full.fw.Q_expr(expr)) for expr in bases[offset]]
            for offset in (-1,0)}
    for expr in images[-1]+images[0]:
        assert not sp_environment.canonical(full.fw.Q_expr(expr))
    matrices={}
    for offset in (-1,0):
        target=bases[offset+1]
        dims=dimensions(mu,length+offset+1,group)
        parts=[]
        expressions=target+images[offset]
        for degree,dimension in dims.items():
            filtered=[{word:c for word,c in expr.items() if nferm(word)==degree} for expr in expressions]
            alphabet=sorted({letter for expr in filtered for word in expr for trace in word for letter in trace})
            frame=GaugeFrame(group,degree,dimension,alphabet,prime,433019+7919*degree+1009*offset)
            values=evaluate([reduce_expr(expr,prime) for expr in filtered],frame)
            assert nmod_mat(values[:len(target)].tolist(),prime).rank()==dimension
            assert nmod_mat(values.tolist(),prime).rank()==dimension
            parts.append(values)
        values=np.concatenate(parts,axis=1) if parts else np.zeros((len(expressions),0),dtype=np.int64)
        matrices[offset]=values[len(target):]
        np.save(output/('incoming.npy' if offset==-1 else 'outgoing.npy'),matrices[offset])
    ranks={offset:int(nmod_mat(matrix.tolist(),prime).rank()) if matrix.shape[1] else 0
           for offset,matrix in matrices.items()}
    remainder=len(bases[0])-ranks[-1]-ranks[0]
    assert remainder>=0
    result=dict(status='complete_classical_rank_bound',group=group,mu=mu,length=length,prime=prime,
                dimensions={str(k):len(v) for k,v in bases.items()},incoming_rank_lower_bound=ranks[-1],
                outgoing_rank_lower_bound=ranks[0],cohomology_upper_bound=remainder,
                zero_cohomology_over_Q_certified=remainder==0,
                exact_formal_nilpotence_checks=sum(map(len,images.values())),
                source_sha256={path.name:digest(path) for path in (Path(__file__),HERE/'lie_root_dimension_v1.py',
                    HERE/'screen_gaugino_source_dimensions_v1.py',HERE/'full_extremal.py')},
                evidence_sha256={str(path.relative_to(output)):digest(path) for path in output.rglob('*')
                                 if path.is_file() and path.name!='result.json'},
                seconds=time.monotonic()-start,
                scope='Exact rational cochain dimensions and finite-field rank lower bounds. Positive remainder is not a class or a quantum lift.',
                d1_increment=0,d2_increment=0,d3_increment=0,new_generator_increment=0)
    if group=='G2':result['lie_algebra_controls']=g2_basis()[1]
    save(output/'result.json',result)
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--group',choices=('G2','Sp2'),default='G2')
    parser.add_argument('--order',type=int,choices=(0,1,2,3),required=True)
    parser.add_argument('--prime',type=int,default=1000003)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=int,default=900)
    args=parser.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('G2 classical screen cap')))
    signal.alarm(args.seconds)
    mu=(3-args.order,3-args.order,args.order,args.order,args.order)
    length=6 if args.order==0 else 5+args.order
    args.output.mkdir(parents=True,exist_ok=True)
    map_ranks(mu,length,args.group,args.prime,args.output)
