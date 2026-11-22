"""Evaluate bracket images in the full (4,4,3,3) target of the classical source test."""
from itertools import combinations
import json
from math import comb
import numpy as np
from scipy.sparse import coo_matrix

from colour import items
from mixed_syzygies import Relations
from source_symmetry import multiplicity, singlet_trace
from b2_schur_complex import transpose

PRIME = 1009


class Planes:
    def __init__(self):
        rel = Relations()
        self.c = rel.c
        self.metric = np.array([[self.c.metric.get((i,j), 0) for j in range(52)] for i in range(52)])
        self.pair_map = {}
        for i,j in combinations(range(52),2):
            terms = {}
            for a, y in self.c.br[i][j].items():
                for mask, z in rel.delta[a].items():
                    pair = items(mask)
                    terms[pair] = terms.get(pair, 0)+y*z
            self.pair_map[52*i+j] = [(p,x) for p,x in terms.items() if x]
        self.binom = np.array([[comb(i,k) for i in range(52)] for k in range(5)])
        self.bases = {0:np.empty((1,0), dtype=np.int16)}
        for k in range(1,5):
            b = np.array(list(combinations(range(52),k)), dtype=np.int16)
            self.bases[k] = b[np.argsort(self.index(b))]

    def index(self, rows):
        return sum((self.binom[j+1,rows[:,j]] for j in range(rows.shape[1])),
                   np.zeros(len(rows), dtype=np.int32))

    def bracket_matrix(self, k):
        basis = self.bases[k]
        columns, rows, values = [], [], []
        for i,j in combinations(range(k),2):
            tags = 52*basis[:,i]+basis[:,j]
            order = np.argsort(tags)
            which, first = np.unique(tags[order], return_index=True)
            for tag, lo, hi in zip(which,first,(*first[1:],len(order))):
                source = order[lo:hi]
                for (a,b), value in self.pair_map[int(tag)]:
                    changed = basis[source].copy()
                    changed[:,i], changed[:,j] = a,b
                    inv = np.zeros(len(changed), dtype=bool)
                    for u,v in combinations(range(k),2):
                        inv ^= changed[:,u] > changed[:,v]
                    changed.sort(axis=1)
                    keep = np.all(changed[:,1:] != changed[:,:-1], axis=1)
                    columns.append(source[keep].astype(np.int32))
                    rows.append(self.index(changed[keep]).astype(np.int32))
                    values.append(value*(1-2*inv[keep].astype(np.int32)))
        matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows),np.concatenate(columns))),
                            shape=(len(basis),len(basis)), dtype=np.int64).tocsr()
        matrix.eliminate_zeros()
        return matrix

    def wedge(self, points, degree):
        previous = np.ones((1,len(points)), dtype=np.int64)
        for k in range(1,degree+1):
            basis = self.bases[k]
            out = np.zeros((len(basis),len(points)),dtype=np.int64)
            for j in range(k):
                indices = self.index(np.delete(basis,j,axis=1))
                out += (-1)**(k-1-j)*previous[indices]*points[:,k-1,basis[:,j]].T
            previous = out % PRIME
        return previous

    def moments(self, k, matrix, points, count):
        dual = self.wedge(points@self.metric,k)
        value = self.wedge(points,k)
        answer = []
        for power in range(count):
            if power:
                value = (matrix@value)%PRIME
            answer.append(np.sum(value*dual,axis=0)%PRIME)
        return np.stack(answer,axis=1)

    def on_cubic(self, points):
        points = points.copy()
        for point in points:
            x,y,z = point[:3]
            bracket = np.zeros(52, dtype=np.int64)
            for i,j in combinations(range(52),2):
                for t,v in self.c.br[i][j].items():
                    bracket[t] += (x[i]*y[j]-x[j]*y[i])*v
            linear = (bracket@self.metric)%PRIME
            j = int(np.flatnonzero(linear)[0])
            # This is the reduction of a rational solution with denominator
            # invertible at PRIME, so the rank also bounds the rational rank.
            z[j] = (z[j]-int(linear@z)*pow(int(linear[j]),-1,PRIME))%PRIME
            assert int(linear@z)%PRIME == 0
        return points


def pivots(array):
    basis = {}
    accepted = []
    for i, source in enumerate(array):
        row = source.copy()%PRIME
        while row.any():
            k = int(np.flatnonzero(row)[0])
            if k not in basis:
                basis[k] = row*pow(int(row[k]),-1,PRIME)%PRIME
                accepted.append(i)
                break
            row = (row-row[k]*basis[k])%PRIME
    return accepted


def calculate():
    algebra = Planes()
    rng = np.random.Generator(np.random.PCG64(443336))
    points = rng.integers(-2,3,size=(36,4,52))
    points[0] = 0
    points[0,:,:4] = np.eye(4,dtype=int)
    ordinary = len(points)
    points = np.concatenate((points,algebra.on_cubic(rng.integers(-2,3,size=(32,4,52)))))
    sizes = (3,4)
    tables = {}
    for k,count in zip(sizes,(5,8)):
        matrix = algebra.bracket_matrix(k)
        tables[k] = np.concatenate([algebra.moments(k,matrix,points[i:i+4],count)
                                    for i in range(0,len(points),4)])
        assert np.count_nonzero(tables[k][0]) == 1 and tables[k][0,0]
    products = (tables[4][:,:,None]*tables[3][:,None,:]%PRIME).reshape(len(points),-1)
    divisor = products[ordinary:]
    products = products[:ordinary]
    ranks = [len(pivots(tables[k][:ordinary])) for k in sizes]
    dimensions = [multiplicity((k,k),2*k,singlet_trace) for k in sizes]
    cols = pivots(products.T)
    rows = pivots(products[:,cols])
    assert ranks == dimensions == [4,6] and len(cols) == len(rows) == 24
    assert len(pivots(products[:,1:])) == 23
    div_cols = pivots(divisor.T)
    div_rows = pivots(divisor[:,div_cols])
    assert len(div_cols) == len(div_rows) == 18
    factors = []
    for a,b in ((1,0),(1,1),(2,0),(2,1),(2,2),(3,1),(3,2)):
        s = transpose(tuple(k for k in (a,a,a,b) if k))
        n = sum(s)
        factors.append(dict(shape=s,dimensions=[multiplicity(s,y,singlet_trace) for y in (n-2,n-1,n)]))
    return dict(prime=PRIME,seed=443336,factor_ranks=ranks,columns=cols,rows=rows,
                minor=products[np.ix_(rows,cols)].tolist(),product_rank=len(cols),
                classical_boundary_rank=23,Cartan_probe=products[0].tolist(),
                cubic_divisor_columns=div_cols,cubic_divisor_rows=div_rows,
                cubic_divisor_minor=divisor[np.ix_(div_rows,div_cols)].tolist(),
                factor_cochain_dimensions=factors)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
