"""Test a second-page dual on alternating pairs of scalar commutators in F4."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
import gzip
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np
from scipy.sparse import csr_matrix

from dotted_filtration import exterior
from flag_boundaries import HERE, recipe
from single_trace_flag import native26
from source_symmetry import centralizer, character, partitions
from f4_singlet_characters import decomposition


LABELS = ((0,0,0,0), (0,0,0,1), (1,0,0,0), (0,0,1,0), (0,0,0,2))
CHANNELS = (1, 26, 52, 273, 324)
PAIRINGS = ((0,1,2,3), (0,2,1,3), (0,3,1,2))
PIN = '3b30b8ff242977debd90581c12f4a767d7875c6c'
SHA = {1009:'009ee18c0e063efe5bc23fbb09aa32b533a02ec1c5cdb0800311634fadf494a9',
       1013:'80764ebac66165672cf6bbe6c0677b25ca52affb78eac1f9d19ee06ca3c8d044'}


def mod(q, p):
    value = Fraction(q)
    return value.numerator * pow(value.denominator, -1, p) % p


def as_array(matrix):
    return np.asarray(matrix.tolist(), dtype=np.int64)


def trace_pair(a, b):
    return int((a * b.T).sum())


def characters():
    small = {k:[exterior(k).get(w, 0) for w in LABELS] for k in (2,3,4)}
    assert list(small.values()) == [[0,0,1,0,0],[1,0,0,0,1],[0,0,1,0,1]]
    mu = (3,3,3,1)
    labels = (LABELS[-1], (2,0,0,0))
    domain = []
    for w in labels:
        n = sum(Fraction(character(mu, c) * decomposition('F4', c).get(w, 0), centralizer(c))
                for c in partitions(10))
        assert n.denominator == 1
        domain.append(int(n))
    assert domain == [316, 331]
    return {'exterior_End26':small, 'remaining_flag':mu, 'symmetric_colour_multiplicities':domain}


class Channels:
    def __init__(self, p, rep, cubic):
        self.p = p
        self.eye = np.eye(26, dtype=np.int64)
        self.inv26 = pow(26, -1, p)
        form = np.zeros((26,26), dtype=np.int64)
        native = json.loads(gzip.decompress((HERE/'representation26.json.gz').read_bytes()))
        for i,j,q in native['form']:
            form[i,j] = mod(q,p)
        self.form = form
        self.inverse = as_array(nmod_mat(form.tolist(),p).inv())
        self.bases = {}
        for size, matrices in ((26,self.inverse @ cubic % p), (52,rep % p)):
            coordinates = matrices.swapaxes(1,2).reshape(size,-1)
            gram = coordinates @ matrices.reshape(size,-1).T % p
            inverse = as_array(nmod_mat(gram.tolist(),p).inv())
            self.bases[size] = csr_matrix(coordinates), inverse

    def decompose(self, A, B, parity):
        p = self.p
        ab = A @ B % p
        ba = B @ A % p
        # The two input matrices are eigenvectors of transpose in the invariant form.
        symmetric = (ab + parity * ba) * pow(2,-1,p) % p
        skew = (ab - symmetric) % p
        frame = {}
        for size, part in ((26,symmetric),(52,skew)):
            coords, _ = self.bases[size]
            frame[size] = coords @ part.ravel() % p
        return symmetric, skew, int(ab.trace()) % p, frame

    def join(self, first, second):
        p = self.p
        a, b, c, x = first
        aa, bb, cc, y = second
        trivial = c * cc * self.inv26 % p
        joined = {size:int((x[size] @ inv % p) @ y[size]) % p for size, (_,inv) in self.bases.items()}
        return (trivial, joined[26], joined[52], (trace_pair(b,bb)-joined[52]) % p,
                (trace_pair(a,aa)-trivial-joined[26]) % p)


def alternating_columns(x, channels, direct):
    p = channels.p
    pairs = {(i,j):(x[i]@x[j]-x[j]@x[i]) % p for i,j in combinations(range(4),2)}
    # Exterior square of the colour-valued bracket, with both orders retained.
    bracket = np.zeros_like(x[0])
    for perm, phase in zip(((0,1,2,3),(0,2,1,3),(0,3,1,2)), (1,-1,1)):
        i,j,k,l = perm
        a,b = pairs[i,j],pairs[k,l]
        bracket += phase * (a@b + b@a)
    bracket %= p

    def triple(ids):
        i,j,k = ids
        return (x[i]@pairs[j,k] - x[j]@pairs[i,k] + x[k]@pairs[i,j]) % p

    three = triple((0,1,2))
    scalar = int(three.trace()) * channels.inv26 % p
    adjoint = sum((-1)**i * x[i] * int(triple(tuple(j for j in range(4) if j != i)).trace())
                  for i in range(4)) % p
    if direct:
        check = np.zeros_like(bracket)
        for perm in permutations(range(4)):
            phase = (-1)**sum(i>j for a,i in enumerate(perm) for j in perm[a+1:])
            v = x[perm[0]]
            for i in perm[1:]:
                v = v @ x[i] % p
            check += phase * v
        assert np.array_equal(check % p, bracket)
        assert not bracket.trace() % p and not adjoint.trace() % p
        for mat, parity in ((bracket,1),(three,1),(adjoint,-1)):
            assert not np.any((channels.form@mat - parity*mat.T@channels.form) % p)
        assert not np.any((channels.bases[26][0] @ bracket.ravel()) % p)
    return {'marked':(bracket,1), 'adjoint':(adjoint,-1), 'scalar':(scalar*channels.eye,1),
            'three':((three-scalar*channels.eye)%p,1)}


def moment_rows(rep, cubic, n, p):
    channels = Channels(p,rep,cubic)
    rows = []
    for point in range(n):
        weights = np.random.Generator(np.random.PCG64(144336+point)).integers(-3,4,(4,52))
        x = (weights @ rep.reshape(52,-1)).reshape(4,26,26) % p
        tensors = alternating_columns(x,channels,point==0)

        @cache
        def block(i,j):
            u,sg = tensors[i]
            v,tg = tensors[j]
            return channels.decompose(u,v,sg*tg)

        column = []
        for other in ('marked','adjoint'):
            for a,b in (('scalar','scalar'),('scalar','three'),('three','three')):
                factors=('marked',other,a,b)
                for i,j,k,l in PAIRINGS:
                    column.extend(channels.join(block(factors[i],factors[j]),block(factors[k],factors[l])))
        rows.append(column)
    return np.asarray(rows,dtype=np.int64).T


def read_minor(path, p):
    info = json.loads((path/f'selected-{p}.json').read_bytes())
    payload = gzip.decompress((path/f'selected-{p}.u16.gz').read_bytes())
    digest = sha256(payload).hexdigest()
    assert digest == info['residue_u16le_sha256'] == SHA[p]
    assert info['points'] == list(range(703))
    assert Counter(t for t,_ in info['rows']) == {'reference':542,'four':158,'nested':3}
    m = nmod_mat(np.frombuffer(payload,dtype='<u2').reshape(703,703).tolist(),p)
    assert int(m.det()) == info['determinant'] == {1009:1005,1013:734}[p]
    return info,m


def quotient_map(info, p):
    traces=json.loads((HERE/'triple_selected.json').read_bytes())['rows']
    short=json.loads((HERE/'flag_ghost_rows.json').read_bytes())
    residual={r['recipe']:r['coordinates'] for r in short['four_trace_recipes']}
    assert len(residual)==46
    rem = np.zeros((703,7), dtype=np.int64)
    nested = []
    for row,(family,index) in enumerate(info['rows']):
        if family=='reference':
            name,seed=traces[index]
            if name=='short' and len(recipe(seed)[0])==4:
                rem[row,:4] = [mod(q,p) for q in residual[seed]]
        elif family=='nested':
            rem[row,4+len(nested)] = 1
            nested.append(index)
    assert len(nested)==3
    result=nmod_mat(rem.tolist(),p)
    assert result.rank()==7
    return result


def run(path, primes, candidate):
    rep=native26()
    from colour_graph_ghosts import sparse_cubic
    _, lower, _, _ = sparse_cubic()
    cubic=np.zeros((26,26,26), dtype=np.int64)
    for (i,j,k), v in lower.items():
        cubic[i,j,k]=v
    result={'basis_commit':PIN, 'characters':characters(), 'primes':{}}
    for p in primes:
        info,m=read_minor(path,p)
        quotient=quotient_map(info,p)
        # All evaluation rows are on the same 703 integral four-colour points.
        weights=m.solve(quotient)
        assert m*weights==quotient
        forms=moment_rows(rep,cubic,703,p)
        projected=nmod_mat(forms.tolist(),p)*weights
        record={'tested_contractions':len(forms), 'evaluation_rank':nmod_mat(forms.tolist(),p).rank(),
                'projected_rank':projected.rank(), 'projected_rows':as_array(projected).tolist(),
                'min_inactive_dimension':696+projected.rank()}
        if candidate and p==1009:
            dual=json.loads(candidate.read_bytes())
            coeff=np.array(dual['weights'],dtype=np.int64)
            assert dual['prime']==p and len(coeff)>=703 and not np.any(coeff[703:])
            record['candidate_nonzero'] = np.flatnonzero(forms @ coeff[:703] % p).tolist()
        result['primes'][str(p)]=record
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('minor_folder', type=Path)
    parser.add_argument('--prime', type=int, choices=(1009,1013), action='append')
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--out',required=True,type=Path)
    a=parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    result=run(a.minor_folder,a.prime or (1009,1013),a.candidate)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print([(p,r['evaluation_rank'],r['projected_rank'],r.get('candidate_nonzero')) for p,r in result['primes'].items()])
