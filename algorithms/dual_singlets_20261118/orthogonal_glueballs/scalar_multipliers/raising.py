"""Transfer lower-page exclusions by the strict flavour-to-dotted substitutions."""
from itertools import combinations, product

from verify import A, GHOSTS, target_coordinates
from mixed_incoming import so_canon
from q1d_exact import pair_terms
from rank_four_lowered_replay import terminal
from so8_level24_descendant import check_letter, substitute, substitution


def raw_word(monomial):
    from verify import PHI
    raw = tuple(tuple(x for x,n in zip(GHOSTS+PHI,trace) for _ in range(n)) for trace in monomial)
    sign,w = so_canon(raw)
    assert sign
    return {w:sign}


def collapse(monomial, spins):
    poly = raw_word(monomial)
    for spin in spins:
        poly = substitute(poly,spin,0,True)
    return sum(c*terminal(word)[0] for word,c in poly.items())


def run():
    jets = {n:[(i,j,t) for i in range(n+1) for j in range(n+1-i)
               for t in combinations(range(3),n-i-j)] for n in (1,2)}
    pairs = list(product(jets[1],jets[1]))+list(product(jets[1],jets[2]))+list(product(jets[2],jets[1]))
    assert len(pairs) == 145
    for spin in range(2):
        for x in jets[1]+jets[2]:
            check_letter(x,spin,0,True)
        for x,y in pairs:
            for ix,z in enumerate((x,y)):
                c,v = substitution(z,spin,0,True)
                if c:
                    assert not (pair_terms(v,y) if ix == 0 else pair_terms(x,v))
            for _,left,right in pair_terms(x,y):
                assert not substitution(left,spin,0,True)[0]
                assert not substitution(right,spin,0,True)[0]
    # Each terminal basis trace has the plus ghost before the minus ghost.
    result = {}
    for p in (6,7):
        basis,_ = target_coordinates((p,p,2,0,0))
        values = [collapse(v,(1,0)) for v in basis]
        assert values == [2,1,0]
        mono,_ = target_coordinates((p+1,p,1,0,0))
        assert len(mono) == 1 and collapse(mono[0],(1,)) == -1
        result[p] = values
    print('34 classical letter tests, 290 ordered pair tests, double substitutions',result)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with assertions enabled')
    run()
