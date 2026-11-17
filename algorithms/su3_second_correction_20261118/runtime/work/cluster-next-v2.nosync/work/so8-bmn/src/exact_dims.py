"""Exact SO(N) invariant multiplicities on symmetric/exterior adjoint letters.

adjoint = exterior^2 defining.  The SO(N) trivial representations inside
Schur_nu(defining) have nu all even, or N odd rows (determinant type).
Power-sum plethysm p_k[e_2]=(p_k^2-p_(2k))/2 and Murnaghan--Nakayama
give the symmetric-group character of invariant tensors. Symmetrization
of bosonic species and antisymmetrization of fermionic species is exact.
All arithmetic is Python integers; no numerical Haar integration.
"""
from functools import lru_cache
from collections import Counter, defaultdict
from itertools import product
from math import factorial


@lru_cache(None)
def partitions(n, largest=None):
    if n == 0:
        return ((),)
    if largest is None or largest > n:
        largest = n
    return tuple((k,) + tail for k in range(largest, 0, -1)
                 for tail in partitions(n-k, k))


@lru_cache(None)
def character(shape, cycles):
    if not cycles:
        return int(not shape)
    if sum(shape) != sum(cycles):
        return 0
    k = cycles[0]
    ell = len(shape)
    beta = tuple(shape[i] + ell - 1 - i for i in range(ell))
    total = 0
    for i, old in enumerate(beta):
        new = old - k
        if new < 0 or new in beta:
            continue
        sign = (-1)**sum(new < b < old for b in beta)
        moved = sorted(beta[:i] + beta[i+1:] + (new,), reverse=True)
        small = tuple(b - ell + 1 + j for j, b in enumerate(moved))
        small = tuple(v for v in small if v)
        total += sign * character(small, cycles[1:])
    return total


@lru_cache(None)
def invariant_trace(cycles, parity, size=8):
    """Trace of a letter permutation on SO(size)-invariant adjoint tensors."""
    degree = sum(cycles)
    if parity == 0:
        shapes = [tuple(2*v for v in x) for x in partitions(degree) if len(x) <= size]
    elif size % 2 == 0 and degree >= size//2:
        shapes = [tuple(2*(x[i] if i < len(x) else 0)+1 for i in range(size))
                  for x in partitions(degree-size//2) if len(x) <= size]
    else:
        return 0
    terms = {(): 1}
    for k in cycles:
        nxt = defaultdict(int)
        for old, c in terms.items():
            nxt[tuple(sorted(old+(k,k), reverse=True))] += c
            nxt[tuple(sorted(old+(2*k,), reverse=True))] -= c
        terms = {k:v for k,v in nxt.items() if v}
    value = sum(c*sum(character(s, rho) for s in shapes) for rho,c in terms.items())
    den = 2**len(cycles)
    assert value % den == 0, (cycles, parity, value, den)
    return value//den


@lru_cache(None)
def cycle_data(n, fermionic):
    ans = []
    for mu in partitions(n):
        z = 1
        for k, m in Counter(mu).items():
            z *= k**m * factorial(m)
        ans.append((mu, (factorial(n)//z)*((-1)**(n-len(mu)) if fermionic else 1)))
    return tuple(ans)


@lru_cache(None)
def multiplicity(content, parity, size=8):
    # Species permutations are symmetric, so canonicalize their multiplicities.
    bos = sorted(content[:3]+(content[6],), reverse=True)
    fer = sorted(content[3:6], reverse=True)
    return _multiplicity(tuple(bos), tuple(fer), parity, size)


@lru_cache(None)
def _multiplicity(bos, fer, parity, size):
    terms = {(): 1}
    den = 1
    for n, odd in [(n,False) for n in bos]+[(n,True) for n in fer]:
        den *= factorial(n)
        nxt = defaultdict(int)
        for old, c in terms.items():
            for mu, v in cycle_data(n, odd):
                nxt[tuple(sorted(old+mu,reverse=True))] += c*v
        terms = {k:v for k,v in nxt.items() if v}
    value = sum(c*invariant_trace(mu,parity,size) for mu,c in terms.items())
    assert value % den == 0
    value //= den
    assert value >= 0, (bos,fer,parity,value)
    return value


def scalar_ring_check(max_degree=16):
    # O(8) even ring degrees 2,4,6,8; odd part = Pfaffian times this ring.
    c = [1]+[0]*max_degree
    for d in (2,4,6,8):
        for k in range(d,max_degree+1):
            c[k] += c[k-d]
    rows=[]
    for n in range(max_degree+1):
        even=multiplicity((n,0,0,0,0,0,0),0)
        odd=multiplicity((n,0,0,0,0,0,0),1)
        assert even==c[n] and odd==(c[n-4] if n>=4 else 0),(n,even,odd,c[n])
        rows.append((n,even,odd))
    return rows


if __name__ == '__main__':
    import time,json
    t=time.monotonic()
    print(json.dumps({'scalar_ring':scalar_ring_check(12),'seconds':time.monotonic()-t}))
