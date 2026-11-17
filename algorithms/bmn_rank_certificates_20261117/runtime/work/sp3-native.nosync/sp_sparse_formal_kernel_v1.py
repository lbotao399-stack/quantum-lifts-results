"""Assemble the existing independent formal-kernel blocks directly as CSR."""
from array import array
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
import hw_blocked


def split_kernel_sparse(images, prime):
    count = len(images)
    parent = list(range(count))
    sizes = [1] * count

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left, right):
        left, right = root(left), root(right)
        if left == right:
            return
        if sizes[left] < sizes[right]:
            left, right = right, left
        parent[right] = left
        sizes[left] += sizes[right]

    owners = {}
    for index, expression in enumerate(images):
        for word, coefficient in expression.items():
            if not int(coefficient) % prime:
                continue
            if word in owners:
                union(index, owners[word])
            else:
                owners[word] = index
    groups = defaultdict(list)
    for index in range(count):
        groups[root(index)].append(index)
    del owners
    data = array('q')
    indices = array('q')
    indptr = array('q', [0])
    rank = 0
    rows = 0
    for ids in groups.values():
        kernel, block_rank = hw_blocked.DENSE_KERNEL([images[index] for index in ids], prime)
        rank += block_rank
        for row in kernel:
            nonzero = np.flatnonzero(row)
            for index in nonzero:
                indices.append(ids[int(index)])
                data.append(int(row[index]) % prime)
            indptr.append(len(data))
            rows += 1
        del kernel
    assert rows == count-rank
    matrix = csr_matrix((np.frombuffer(data, dtype=np.int64),
                         np.frombuffer(indices, dtype=np.int64),
                         np.frombuffer(indptr, dtype=np.int64)), shape=(rows, count))
    matrix.sort_indices()
    return matrix, rank
