"""Check the full rational kernel of the two marked three-exchange families."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
import gzip
from itertools import permutations
import json
from math import gcd
from multiprocessing import get_context
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix, eye, hstack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from seventh_crossings import PAIRS, native26
from integer_section import HERE as SECTION, tensors
from single_trace_flag import __file__ as trace_module

RELATIONS = json.loads((HERE / "relations.json").read_bytes())["relations"]


def row_norm(matrix):
    return int(np.abs(matrix).sum(axis=1).max())


def native_operators():
    representation = json.loads(gzip.decompress(
        (Path(trace_module).parent / "representation26.json.gz").read_bytes()))
    form = np.zeros((26, 26), dtype=np.int64)
    for i, j, value in representation["form"]:
        entry = Fraction(value)
        assert entry.denominator == 1
        form[i, j] = entry.numerator
    dual = 3 * form.copy()
    dual[:2, :2] = [[2, 1], [1, 2]]
    assert np.array_equal(dual @ form, 3*np.eye(26, dtype=np.int64))
    cubic = np.zeros((26,)*3, dtype=np.int64)
    for i, j, k, z in json.loads((SECTION/"cubic.json").read_bytes())["cubic"]:
        cubic[i, j, k] = z
    gen = native26().astype(np.int64)
    assert np.array_equal(form @ gen, -gen.transpose(0, 2, 1) @ form)
    jordan = dual @ cubic
    two = gen[:,None] @ gen[None,:]
    assert not np.any(two.reshape(52**2, -1) @ jordan.transpose(0,2,1).reshape(26,-1).T)

    _, upper, ex, _ = tensors()
    sym = sum(upper.astype(np.int64).transpose((0,*p)) for p in permutations((1,2,3)))
    paired = (np.einsum("ij,kl->ijkl",dual,dual) + np.einsum("ik,jl->ijkl",dual,dual)
              + np.einsum("il,jk->ijkl",dual,dual))
    assert np.array_equal(sym,4*paired)
    exchange = csr_matrix(ex, dtype=np.int64)
    ident = eye(676, dtype=np.int64, format="csr")
    trace = csr_matrix(np.outer(dual.ravel(), form.ravel()))
    swap = ident[np.arange(676).reshape(26,26).T.ravel()]
    sym_free = 39*(ident+swap)-trace
    assert (sym_free@sym_free-78*sym_free).nnz == 0
    eigen_op = ident
    denominator = 1
    for value in (-36,-12,6,84):
        eigen_op = eigen_op @ (exchange-value*ident)
        denominator *= 42-value
    assert (eigen_op@exchange-42*eigen_op).nnz == 0
    assert int(eigen_op.diagonal().sum()) == 52*denominator
    # These 52 independent native matrices have eigenvalue 42 in the raised pair.
    lie = (gen @ dual).reshape(52, 676).T
    assert np.array_equal(exchange @ lie, 42*lie)
    # A nondegenerate trace form on the independent generators checks the dimension.
    assert modular_rank(lie.T @ lie,1009)[0] == 52
    scale = abs(denominator*78)//gcd(denominator,78)
    left23 = scale*ident-(scale//denominator)*eigen_op-(scale//78)*trace
    left23.data //= gcd(*map(int,left23.data))
    return gen, exchange, cubic, 78*ident-trace, sym_free, left23


def modular_rank(data, prime):
    """Row rank, and determinant for full-rank square matrices."""
    matrix = np.asarray(data, dtype=np.int64).copy()%prime
    rank, determinant = 0, 1
    for j in range(matrix.shape[1]):
        pivots = np.flatnonzero(matrix[rank:,j])
        if not len(pivots):
            continue
        i = rank + int(pivots[0])
        if i != rank:
            matrix[[i,rank]] = matrix[[rank,i]]
            determinant = -determinant
        value = int(matrix[rank,j])
        determinant = determinant*value%prime
        mult = matrix[rank+1:,j] * pow(value,-1,prime)%prime
        matrix[rank+1:,j:] = (matrix[rank+1:,j:]-mult[:,None]*matrix[rank,j:])%prime
        rank += 1
        if rank == len(matrix):
            break
    return rank, determinant if rank==len(matrix)==matrix.shape[1] else None


def endpoint_subspaces():
    gen, exchange, cubic, tracefree, sym_free, left23 = native_operators()
    h = gen[:4]
    assert not np.any(h * (1-np.eye(26,dtype=np.int64)))
    weights = np.diagonal(h,axis1=1,axis2=2).T
    all_slots = np.stack(np.unravel_index(np.arange(26**4),(26,)*4)).T
    sector = all_slots[np.all(weights[all_slots].sum(axis=1)==0, axis=1)]
    assert len(sector)==5512
    strides = np.array([26**3,26**2,26,1])
    reverse = np.full(26**4,-1)
    reverse[sector@strides] = np.arange(len(sector))

    def insert(op, pair):
        i,j = pair
        by_input = op.tocsc()
        rows,cols,values = [],[],[]
        for column,t in enumerate(sector):
            index = 26*t[i]+t[j]
            start,stop = by_input.indptr[index:index+2]
            entry = by_input.indices[start:stop]
            other = t@strides-t[i]*strides[i]-t[j]*strides[j]
            r = reverse[other+(entry//26)*strides[i]+(entry%26)*strides[j]]
            assert np.all(r>=0)
            rows.extend(r)
            cols.extend([column]*len(r))
            values.extend(by_input.data[start:stop])
        return csr_matrix((values,(rows,cols)),shape=(len(sector),len(sector)),dtype=np.int64)

    left = insert(tracefree,(0,1)) @ insert(left23,(2,3))
    right = (insert(tracefree,(0,2)) @ insert(sym_free,(1,3))).tocsc()
    cubic_trace = csr_matrix((cubic[tuple(sector[:,1:].T)],(sector[:,0],np.arange(len(sector)))),
                             shape=(26,len(sector)))
    constraints = (cubic_trace @ right).toarray()
    active = constraints[constraints.any(axis=1)]
    assert len(active)==2
    # Two pivot columns parametrise the exact kernel of the remaining cubic contraction.
    choices = {tuple(map(int,v)):i for i,v in enumerate(active.T) if any(v)}
    _, i, j = min((abs(a*d-b*c),i,j) for (a,b),i in choices.items()
                  for (c,d),j in choices.items() if a*d!=b*c)
    (a,b),(c,d) = active[:,i].tolist(),active[:,j].tolist()
    determinant = a*d-b*c
    columns = []
    for k,(u,v) in enumerate(active.T.tolist()):
        if k in (i,j):
            continue
        x,y = d*u-c*v,a*v-b*u
        common = gcd(determinant,x,y)
        col = ((determinant//common)*right[:,k]-(x//common)*right[:,i]-(y//common)*right[:,j])
        if col.nnz:
            col.data //= gcd(*map(int,col.data))
        columns.append(col)
    right = hstack(columns,format="csr")
    assert (cubic_trace@right).nnz==0
    return [insert(exchange,pair) for pair in PAIRS], left, right


@cache
def blocks():
    return endpoint_subspaces()


def check_columns(interval):
    ops, left, right = blocks()
    start, stop = interval
    level = [right[:,start:stop]]
    words = level.copy()
    for _ in range(3):
        level = [op @ previous for previous in level for op in ops]
        words.extend(level)
    bound = max(row_norm(op) for op in ops)**3 * int(np.abs(right.data).max())
    assert bound < 2**63 and len(words)==259
    left_norm = row_norm(left)
    for row in RELATIONS:
        prefix_bound = sum(abs(coefficient) for _,coefficient in row)*bound
        if prefix_bound < 2**63:
            residual = sum((coef*words[i] for i,coef in row),csr_matrix(words[0].shape,dtype=np.int64))
            if not residual.nnz:
                continue
            full_bound = left_norm*int(np.abs(residual.data).max())
            if full_bound < 2**63:
                assert not np.any((left@residual).data)
                continue
        else:
            full_bound = prefix_bound*left_norm
        modulus = 1
        for prime in (1000003,1000033,1000037,1000039,1000081,1000099):
            assert gcd(modulus,prime)==1 and len(row)*prime**2<2**63 and left_norm*prime<2**63
            residual = csr_matrix(words[0].shape,dtype=np.int64)
            for i,coef in row:
                reduced = words[i].copy()
                reduced.data %= prime
                residual = residual+(coef%prime)*reduced
            residual.data %= prime
            assert not np.any((left@residual).data % prime)
            modulus *= prime
            if modulus>2*full_bound:
                break
        assert modulus>2*full_bound
    return stop-start


def verify(workers, points_dir):
    relations = np.zeros((len(RELATIONS),259),dtype=np.int64)
    for i,row in enumerate(RELATIONS):
        for j,coef in row:
            relations[i,j] = coef
    assert modular_rank(relations,1009)[0]==68
    minors = {}
    first = json.loads((points_dir/"third-exchange-1009.json").read_bytes())
    for prime in (1009,1013):
        record = json.loads((points_dir/f"third-exchange-{prime}.json").read_bytes())
        assert record["points"]==first["points"] and record["columns"]==first["columns"]
        with np.load(points_dir/f"third-exchange-{prime}.npz",allow_pickle=False) as data:
            closed, free = data["closed"].astype(np.int64),data["free"].astype(np.int64)
        assert closed.shape==(191,259) and free.shape==(36,259)
        assert not np.any(closed @ (relations % prime).T % prime)
        assert not np.any(free @ (relations % prime).T % prime)
        rank, det = modular_rank(closed[:,record["columns"]],prime)
        assert rank==191 and det==record["determinant"] and det!=0
        minors[prime] = det
    _, _, right = blocks()
    n = right.shape[1]
    intervals = [(k,min(k+28,n)) for k in range(0,n,28)]
    if workers==1:
        checked = sum(map(check_columns,intervals))
    else:
        with ProcessPoolExecutor(workers,mp_context=get_context("fork")) as pool:
            checked = sum(pool.map(check_columns,intervals))
    assert checked==5510
    return dict(zero_weight_dimension=5512, cubic_rows=2, parametrising_columns=checked,
                independent_identities=68, exchange_words=259, marked_recipes=518,
                outgoing_minors=minors, rational_outgoing_rank=191, formal_kernel=327,
                quotient_by_boundaries=0)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--points-dir",type=Path,default=HERE)
    args=parser.parse_args()
    print(json.dumps(verify(args.workers,args.points_dir),indent=2))
