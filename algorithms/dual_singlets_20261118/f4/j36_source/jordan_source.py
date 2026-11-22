"""Test three polarized Albert-algebra tree images at the actual j=36 source weight."""
import argparse
from functools import cache
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import sys
from zipfile import ZipFile

from flint import nmod_mat
import numpy as np

from dressed_tree import colour_samples, columns, evaluation_rows
from mixed_colour_relations import pivots
from verify import ARCHIVE, HERE

sys.path.insert(0,str(HERE.parents[1]))
from f4_colour import basis, jordan_colours


@cache
def column_antisymmetrizer():
    p4, p3 = tuple(permutations(range(4))), tuple(permutations(range(3)))
    perm_sign = {p: (-1)**sum(a>b for j,a in enumerate(p) for b in p[j+1:])
                 for p in p4+p3}
    index = []
    for first,second,third,fourth in product(p4,p4,p3,p3):
        cols = (first,second,third,fourth)
        a = sum(col[0]*4**(3-i) for i,col in enumerate(cols))
        b = sum(col[1]*4**(3-i) for i,col in enumerate(cols))
        u = 4*first[2]+first[3]
        v = 4*second[2]+second[3]
        c, d = third[2], fourth[2]
        sign = perm_sign[first]*perm_sign[second]*perm_sign[third]*perm_sign[fourth]
        index.append((a,b,u,v,c,d,sign))
    assert len(index) == 20736
    return np.array(index,dtype=np.int64).T


def coefficients(colours, jordan, inverse, prime):
    dim = colours.shape[1]
    pairs = (colours[:,None] @ colours[None,:]) % prime
    products = pairs.reshape(16,dim,dim)
    quartics = ((products[:,None] @ products[None,:]) % prime).reshape(256,dim,dim)
    projected = np.einsum('iab,tba->ti',jordan,quartics) % prime
    jordan_pairs = ((projected*inverse % prime) @ projected.T) % prime
    brackets = ((pairs-pairs.swapaxes(0,1)) % prime).reshape(16,dim,dim)
    bracket_pairs = (brackets[:,None] @ brackets[None,:]) % prime
    comm_gram = np.trace(bracket_pairs,axis1=-2,axis2=-1) % prime
    form = np.trace(pairs,axis1=-2,axis2=-1) % prime
    adjacent = np.einsum('uvij,cdji->uvcd',bracket_pairs,pairs) % prime
    mixed = (brackets[:,None] @ colours[None,:]) % prime
    separated = np.einsum('ucij,vdji->uvcd',mixed,mixed) % prime
    a,b,u,v,c,d,sign = column_antisymmetrizer()
    bridge = jordan_pairs[a,b]*sign
    tails = (comm_gram[u,v]*form[c,d] % prime, adjacent[u,v,c,d], separated[u,v,c,d])
    assert dim**2*(prime-1)**2 < 2**63 and len(sign)*prime < 2**63
    return [int(((bridge*tail) % prime).sum() % prime) for tail in tails]


def run(archive_path, prime):
    assert sha256(archive_path.read_bytes()).hexdigest() == ARCHIVE
    derivations, _ = basis()
    jordan = jordan_colours()
    gram = np.einsum('iab,jba->ij',jordan,jordan)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    inverse = np.array([pow(int(v),-1,prime) for v in gram.diagonal()],dtype=np.int64)
    data = json.loads((HERE/'dressed-tree.json').read_bytes())
    free = {row['free'] for row in data['kernel']}
    original = [i for i in range(80) if i not in free]
    with ZipFile(archive_path) as archive:
        old = [row for i,row in enumerate(columns(archive)) if i not in free]
        previous = evaluation_rows(archive,old,prime,70,derivations)
    added = [coefficients(q,jordan,inverse,prime) for q in colour_samples(derivations,prime,70)]
    assert len(original) == 61
    rows = [a+b for a,b in zip(previous,added,strict=True)]
    matrix = nmod_mat(rows,prime)
    indices = pivots(matrix)
    assert indices[:61] == list(range(61))
    chosen = pivots(matrix.transpose())
    determinant = int(nmod_mat([[rows[i][j] for j in indices] for i in chosen],prime).det())
    assert determinant
    return dict(prime=prime, old_parameter_columns=original, old_rank=61,
                augmented_rank_lower_bound=len(indices), extra_pivots=[i-61 for i in indices[61:]],
                minor_rows=chosen, determinant=determinant, new_columns=added)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('--prime',type=int,choices=(1000003,1000033),required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    print(json.dumps(run(args.archive,args.prime),separators=(',',':')))
