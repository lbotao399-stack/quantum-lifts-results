"""Integer lowest-weight maps on the rank-one Cartan at fixed Taylor degrees."""

from collections import defaultdict
from functools import lru_cache
from itertools import product

from flint import nmod_mat


def odd(letter):
    return (letter[2].bit_count() + 1) % 2


def weight(letter):
    a, b, mask = letter
    return (a, b, *((mask >> i) & 1 for i in range(3)))


def monomials(t, length):
    target = (t, t, 4-t, 4-t, 4-t)
    letters = [(a, b, mask) for a, b, mask in product(range(t+1), range(t+1), range(8))
               if a+b+mask and all(x <= y for x, y in zip(weight((a,b,mask)), target))]
    weights = [weight(letter) for letter in letters]

    @lru_cache(None)
    def collect(start, slots, remaining):
        if slots == 0:
            return ((),) if not any(remaining) else ()
        if sum(remaining) < slots:
            return ()
        answer = []
        for i in range(start, len(letters)):
            rest = tuple(x-y for x, y in zip(remaining, weights[i]))
            if min(rest) >= 0:
                answer.extend((letters[i],) + others for others in
                              collect(i+odd(letters[i]), slots-1, rest))
        return tuple(answer)

    basis = collect(0, length, target)
    collect.cache_clear()
    return basis


@lru_cache(None)
def letter_images(letter):
    a, b, mask = letter
    degree = a+b+mask.bit_count()
    images = []
    for src in (1, 2):
        if mask & (1 << src) and not mask & (1 << (src-1)):
            images.append((f'E{src}', False, 1, (a, b, mask ^ (3 << (src-1)))))
    if a:
        images.append(('D', False, b+1, (a-1, b+1, mask)))
    for axis, exponent in enumerate((a, b)):
        if exponent == 0:
            continue
        x, y = a-int(axis == 0), b-int(axis == 1)
        if degree > 1:
            images.append((f'K{axis}', False, degree-1, (x, y, mask)))
        for flavour in range(3):
            bit = 1 << flavour
            if not mask & bit:
                sign = (-1)**(mask.bit_count() + (mask & (bit-1)).bit_count())
                images.append((f'T{axis}{flavour}', True, sign, (x, y, mask | bit)))
    if degree > 1:
        for flavour in range(3):
            bit = 1 << flavour
            if mask & bit:
                sign = (-1)**(mask.bit_count() + (mask & (bit-1)).bit_count())
                images.append((f'S{flavour}', True, (degree-1)*sign, (a, b, mask ^ bit)))
    return tuple(images)


def column(monomial):
    image = defaultdict(int)
    prefix_parity = 0
    for i, letter in enumerate(monomial):
        for operator, is_odd, coeff, new in letter_images(letter):
            target = monomial[:i] + (new,) + monomial[i+1:]
            fermions = [l for l in target if odd(l)]
            if len(fermions) != len(set(fermions)):
                continue
            swaps = sum(x > y for j, x in enumerate(fermions) for y in fermions[j+1:])
            sign = (-1)**(swaps + is_odd * prefix_parity)
            image[operator, tuple(sorted(target))] += sign * coeff
        prefix_parity ^= odd(letter)
    return {key: value for key, value in image.items() if value}


def check(t, length, prime=1000003):
    basis = monomials(t, length)
    columns = [column(m) for m in basis]
    labels = sorted(set().union(*(col.keys() for col in columns)))
    label_ids = {key: i for i, key in enumerate(labels)}
    pivots = {}
    chosen = []
    for original in columns:
        work = {label_ids[key]: v % prime for key, v in original.items() if v % prime}
        while work:
            i = min(work)
            coeff = work[i]
            if i not in pivots:
                inv = pow(coeff, -1, prime)
                pivots[i] = {j: v*inv % prime for j, v in work.items()}
                chosen.append(labels[i])
                break
            for j, v in pivots[i].items():
                value = (work.get(j, 0)-coeff*v) % prime
                if value:
                    work[j] = value
                else:
                    work.pop(j, None)
    assert len(chosen) == len(basis)
    minor = nmod_mat([[col.get(key, 0) % prime for col in columns] for key in chosen], prime)
    determinant = int(minor.det())
    assert determinant
    return dict(t=t, length=length, columns=len(basis), rows=len(labels),
                prime=prime, determinant=determinant, kernel=0)


if __name__ == '__main__':
    import json
    for t, length in ((2, 4), (2, 6), (4, 4)):
        print(json.dumps(check(t, length)), flush=True)
