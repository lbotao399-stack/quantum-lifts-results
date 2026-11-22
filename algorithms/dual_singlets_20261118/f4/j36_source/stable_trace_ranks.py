"""Multilinear classical orthogonal traces, through two excess Taylor slots."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations


def odd(block):
    return (1 + block.bit_count()) % 2


def orientation(word):
    if len(word) < 2:
        return (), 0
    pivot = word.index(min(word))
    rotated = word[pivot:] + word[:pivot]
    sign = (-1)**(sum(map(odd, word[:pivot])) * sum(map(odd, word[pivot:])))
    reverse = rotated[:1] + rotated[:0:-1]
    if reverse < rotated:
        count = sum(map(odd, rotated[1:]))
        return reverse, sign * (-1)**(len(word) + count*(count-1)//2)
    return rotated, sign


def partitions(labels, excess):
    if excess == 0:
        yield tuple(1 << i for i in labels)
        return
    for size in range(2, excess+2):
        for group in combinations(labels, size):
            block = sum(1 << i for i in group)
            rest = tuple(i for i in labels if i not in group)
            for parts in partitions(rest, excess-size+1):
                if all(p.bit_count() == 1 or block < p for p in parts):
                    yield (block,) + parts


def basis(degree, excess):
    words = []
    for parts in partitions(tuple(range(degree)), excess):
        if len(parts) < 2:
            continue
        first, *rest = sorted(parts)
        for tail in permutations(rest):
            word = (first,) + tail
            if word <= word[:1] + word[:0:-1]:
                words.append(word)
    return words


@cache
def splits(block):
    """Native transposed divided jets, with all Taylor directions odd and labelled."""
    result = []
    left = (block-1) & block
    while left:
        right = block ^ left
        inversions = sum((right & ((1 << i)-1)).bit_count()
                         for i in range(block.bit_length()) if left & (1 << i))
        exponent = block.bit_count() + (1+left.bit_count())*right.bit_count() + inversions
        result.append((left, right, -((-1)**exponent)))
        left = (left-1) & block
    return tuple(result)


def derivative(word):
    result = Counter()
    sign = 1
    for i, block in enumerate(word):
        for left, right, factor in splits(block):
            key, orientation_sign = orientation(word[:i] + (left, right) + word[i+1:])
            result[key] += sign * orientation_sign * factor
        if odd(block):
            sign = -sign
    return {w: v for w, v in result.items() if v}


def rank_mod_three(columns, rows):
    """A signed pair of bitsets stores the coefficients +1 and -1 modulo three."""
    pivots = {}
    for column in columns:
        positive = negative = 0
        for key, value in column.items():
            residue = value % 3
            if residue == 1:
                positive |= 1 << rows[key]
            elif residue == 2:
                negative |= 1 << rows[key]
        while positive or negative:
            leading = (positive | negative).bit_length()-1
            if leading not in pivots:
                if negative >> leading:
                    positive, negative = negative, positive
                pivots[leading] = positive, negative
                break
            plus, minus = pivots[leading]
            if positive >> leading:
                plus, minus = minus, plus
            own, other = positive | negative, plus | minus
            positive, negative = ((positive & ~other) | (plus & ~own) | (negative & minus),
                                  (negative & ~other) | (minus & ~own) | (positive & plus))
    return len(pivots)


def check(degree):
    endpoints, middle, incoming = [basis(degree, excess) for excess in range(3)]
    assert all(len(b) == len(set(b)) for b in (endpoints, middle, incoming))
    outgoing = {w: derivative(w) for w in middle}
    first_rank = rank_mod_three(outgoing.values(), {w:i for i,w in enumerate(endpoints)})

    def columns():
        for word in incoming:
            differential = derivative(word)
            twice = Counter()
            for mid, value in differential.items():
                for key, factor in outgoing[mid].items():
                    twice[key] += value * factor
            assert not any(twice.values())
            yield differential

    second_rank = rank_mod_three(columns(), {w:i for i,w in enumerate(middle)})
    terminal = len(endpoints)-first_rank
    one_excess = len(middle)-first_rank-second_rank
    if degree % 2:
        assert terminal == 0
    else:
        assert one_excess == 0
    return dict(slots=degree, dimensions=[len(incoming),len(middle),len(endpoints)],
                ranks_mod_3=[second_rank,first_rank],
                upper_bounds=[one_excess,terminal])


if __name__ == '__main__':
    import json
    if not __debug__:
        raise SystemExit('run with assertions enabled')
    for degree in range(1,9):
        print(json.dumps(check(degree)), flush=True)
