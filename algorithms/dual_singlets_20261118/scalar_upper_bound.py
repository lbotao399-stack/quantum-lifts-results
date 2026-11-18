"""Exact Cartan witness for the length-12 quotient used by scalar_neighbors.py."""
import json
from itertools import permutations


def compute():
    result = {}
    for group, n in [('SO(7)', 7), ('USp(6)', 6)]:
        matrices = [[[0] * n for _ in range(n)] for _ in range(3)]
        for i, m in enumerate(matrices):
            if n == 7:
                m[2*i][2*i+1], m[2*i+1][2*i] = 1, -1
            else:
                m[i][i], m[i+3][i+3] = 1, -1
        prod = lambda a, b: [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        assert all(prod(a,b)==prod(b,a) for a in matrices for b in matrices)
        gram = [[sum(prod(a,b)[i][i] for i in range(n)) for b in matrices] for a in matrices]
        det = sum((-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1,3)) * gram[0][p[0]] * gram[1][p[1]] * gram[2][p[2]] for p in permutations(range(3)))
        assert det**2 == 64
        result[group] = dict(trace_pairings=gram, determinant=det, length_12_value=det**2)
    return result

if __name__ == '__main__':
    print(json.dumps(compute(), indent=2))
