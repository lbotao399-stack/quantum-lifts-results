"""Test the restricted scalar classical equation in the factors of subgroups of F4."""
import argparse
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

import mixed_incoming as m
from charge_maps import dimensions as classical_dimensions
from g2_full_classical_screen_v1 import dimensions as exceptional_dimensions, g2_basis
from mixed33_exact import read
from rank_four_s5 import decode_rows
from rank_four_sixth import G, IntegralFrame


class IntegralLetters(IntegralFrame):
    def mat_mul_letter(self, arr, degree, letter):
        parity, right = self.L[letter]
        if parity:
            return super().mat_mul_letter(arr, degree, letter)
        return arr @ right, degree


def generators(group):
    if group == "G2":
        return g2_basis()[0]
    if group == "SU3":
        elementary = np.eye(9,dtype=np.int64).reshape(3,3,3,3)
        return np.stack([elementary[i,j] for i in range(3) for j in range(3) if i!=j]
                        + [elementary[0,0]-elementary[1,1], elementary[1,1]-elementary[2,2]])
    if group == "SO7":
        elementary = np.eye(49, dtype=np.int64).reshape(7, 7, 7, 7)
        return np.array([elementary[i,j]-elementary[j,i] for i in range(7) for j in range(i)])
    assert group == "Sp3"
    elementary = np.eye(9, dtype=np.int64).reshape(3, 3, 3, 3)
    zero = np.zeros((3,3), dtype=np.int64)
    out = [np.block([[a,zero],[zero,-a.T]]) for a in elementary.reshape(9,3,3)]
    for i in range(3):
        for j in range(i+1):
            b = elementary[i,j].copy()
            b[j,i] = 1
            out.extend((np.block([[zero,b],[zero,zero]]), np.block([[zero,zero],[b,zero]])))
    return np.stack(out)


def dimension(group, weight, length, degree):
    if group == "G2":
        return exceptional_dimensions(weight, length, group).get(degree, 0)
    return classical_dimensions("Sp" if group == "Sp3" else "SO_odd", weight, length).get(degree, 0)


def integral_rows(rows, group, degree, dim, seed, workers, basis=None, sample_terms=4):
    words = sorted({w for row in rows for w in row})
    assert all(m.nferm(w) == degree for w in words)
    alphabet = sorted({l for w in words for tr in w for l in tr})
    if basis is None:
        basis = generators(group)
    size = basis.shape[1]
    k = degree+2
    points = math.ceil((dim+35)/math.comb(k,degree))
    frame = IntegralLetters(size,points,k,[degree],alphabet,seed=1)
    rng = np.random.default_rng(seed)
    norm = 0
    for letter in alphabet:
        parity, sample = frame.L[letter]
        count = points*(k if parity else 1)
        elements = []
        for _ in range(count):
            ids = rng.choice(len(basis), size=sample_terms, replace=False)
            signs = rng.choice([-1,1], size=sample_terms)
            elements.append(np.einsum("i,ijk->jk",signs,basis[ids]))
        matrices = np.asarray(elements).reshape(sample.shape)
        norm = max(norm, int(abs(matrices).sum(axis=-1).max()))
        frame.L[letter] = parity, matrices
    bound = max(math.factorial(degree)*size**len(w)*norm**sum(map(len,w)) for w in words)
    if bound >= 2**63:
        necklaces = {t for w in words for t in w}
        assert all(math.factorial(degree)*size*norm**len(t) < 2**63 for t in necklaces)
        limits = {t:(int(abs(v).max()), d) for t,(v,d) in frame.necklace_values(necklaces).items()}
        bound = 0
        for word in words:
            prefix, odd = 1, 0
            for t in word:
                maximum, d = limits[t]
                prefix *= math.comb(odd+d,d)*maximum
                odd += d
                bound = max(bound,prefix)
    assert bound < 2**63
    chunks = np.array_split(np.arange(len(words)), min(len(words),workers*2))
    with ThreadPoolExecutor(workers) as pool:
        values = np.concatenate(list(pool.map(
            lambda inds: frame.eval_words([words[i] for i in inds],dtype=np.int64), chunks)))
    assert int(abs(values).max()) <= bound
    # Every term in an exterior coordinate is bounded by the submultiplicative
    # matrix row norm. Apply the singlet coefficients with arbitrary-precision integers.
    positions = {w:i for i,w in enumerate(words)}
    coeff = fmpz_mat(len(rows),len(words))
    for i,row in enumerate(rows):
        for w,c in row.items():
            coeff[i,positions[w]] = c
    evaluated = coeff*fmpz_mat(values.tolist())
    ref = m.Frame(size,points,k,[degree],alphabet,seed=1,p=1000003)
    ref.L = {l:(d,a%ref.p) for l,(d,a) in frame.L.items()}
    tests = np.unique(np.linspace(0,len(words)-1,min(8,len(words)),dtype=int))
    np.testing.assert_array_equal(values[tests]%ref.p, ref.eval_words([words[i] for i in tests],dtype=np.int64))
    return evaluated, dict(seed=seed, points=points, degree=degree, exterior_generators=k,
                           words=len(words), row_norm=norm, word_bound=str(bound), control_words=tests.tolist())


def run(group, power, plan, workers):
    weight = (power,power,1,1,1)
    length = sum(weight)-1
    parent = (read(plan) if plan else m.prepare(weight, length))
    assert tuple(parent["mu"]) == weight and parent["length"] == length
    degree = 2*power-1
    sources = decode_rows(parent,"sources")
    ids = [i for i,r in enumerate(sources) if m.nferm(next(iter(r))) == degree]
    assert all(all(m.nferm(w)==degree for w in sources[i]) for i in ids)
    src = [sources[i] for i in ids]
    images = [m.q(r) for r in src]
    end = decode_rows(m.prepare(weight,length+1),"sources")
    sign, rhs = m.so_canon((G,)*power+(((0,0,(0,)),(0,0,(1,)),(0,0,(2,))),))
    assert sign and not m.q({rhs:sign})
    ds = dimension(group,weight,length,degree)
    dt = dimension(group,weight,length+1,degree+1)
    a, af = integral_rows(src,group,degree,ds,540907+power,workers)
    v, vf = integral_rows(end+images+[{rhs:sign}],group,degree+1,dt,921871+power,workers)
    vv = v.tolist()
    boundary = fmpz_mat(vv[len(end):-1])
    with_rhs = fmpz_mat(vv[len(end):])
    ranks = [a.rank(),fmpz_mat(vv[:len(end)]).rank(),boundary.rank(),with_rhs.rank()]
    assert ranks[0]<=ds and ranks[1]<=dt and ranks[2]<=ranks[3]<=ranks[2]+1
    result = dict(group=group,power=power,weight=weight,source_character=ds,target_character=dt,
                  ranks=ranks,source_frame=af,target_frame=vf,source_ids=ids,
                  rhs=[m.encode(rhs),sign],source_frame_values_sha256=m.digest(m.compact([[int(x) for x in row] for row in a.tolist()])))
    print(group,power,"dimensions",ds,dt,"ranks",ranks,flush=True)
    if ranks[2] < ranks[3] and ranks[0]==ds:
        # Any evaluation of a full independent source set suffices for separation.
        from f4_pure_target import covector
        result["separating_covector"] = covector(with_rhs)
        result["membership"]="excluded over Q"
    elif ranks[2]==ranks[3] and ranks[1]==dt:
        # Evaluation is injective on the entire terminal singlet by its character.
        equations = fmpq_mat(with_rhs).transpose()
        red, rank = equations.rref()
        co = {}
        for i in range(rank):
            j=next(j for j in range(red.ncols()) if red[i,j])
            assert j<len(ids)
            if red[i,red.ncols()-1]:
                co[j]=red[i,red.ncols()-1]
        weighted = fmpq_mat([[co.get(j,0) for j in range(len(ids))]])*fmpq_mat(boundary)
        assert weighted.tolist()[0] == vv[-1]
        result["primitive_in_source_rows"]=[[ids[j],str(c)] for j,c in co.items()]
        result["membership"]="included over Q"
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group",choices=("G2","Sp3","SO7"),required=True)
    parser.add_argument("--power",type=int,required=True)
    parser.add_argument("--plan",type=Path)
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    m.write(args.out,run(args.group,args.power,args.plan,args.workers))
