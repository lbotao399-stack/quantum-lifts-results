"""Rebuild the [3,1] relation and its preimage from Kim et al., Sec. 4.1 and App. A."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

# Coefficient in the level-34 syzygy and the two rows of its Appendix A contraction matrix.
# Each digit identifies a slot in the order of the tensors printed in that appendix.
TERMS = (
    (1,-914152,'12367','12345'), (2,2468146,'12867','12345'),
    (6,-2724892,'12436','12345'), (8,171002,'12567','12456'),
    (9,-1959484,'12367','12456'), (12,-1139138,'17356','12456'),
    (13,4209083,'18356','12456'), (15,381609,'12346','12345'),
    (16,-610948,'17246','12345'), (17,-4802651,'12846','12345'),
    (18,-1671855,'12835','12345'), (19,-1442031,'16235','12345'),
    (20,1581135,'18235','12345'), (22,283041,'16835','12345'),
    (23,2620424,'12837','12345'), (32,-94806,'12356','23456'),
    (33,106776,'12345','23456'), (34,135828,'17345','23456'),
    (35,-121302,'12436','12356'), (36,111096,'12736','12356'),
    (37,199920,'2','1'), (38,-97104,'4','1'), (40,237048,'4','2'),
)


def free_product(left, right):
    total = Counter()
    for first, a in left.items():
        for second, b in right.items():
            odd = [v for v in first + second if v[0] in ('v2','v3')]
            if len(set(odd)) != len(odd):
                continue
            inversions = sum(x > y for i,x in enumerate(odd) for y in odd[i+1:])
            total[tuple(sorted(first + second))] += (-1)**inversions * a * b
    return {w:c for w,c in total.items() if c}


def use_free_gravitons(dmk):
    def raw(name, *idx):
        return {((name, *idx),): 1}
    def u2(i,j):
        return raw('u2', *sorted((i,j)))
    def u3(i,j,k):
        return raw('u3', *sorted((i,j,k)))
    def v2(i,j):
        return dmk.add((1,raw('v2',i,j)), *[
            (Fraction(-1,3),raw('v2',k,k)) for k in range(3) if i == j])
    def raw3(i,j,k):
        return raw('v3',*sorted((i,j)),k)
    def v3(i,j,k):
        return dmk.add((1,raw3(i,j,k)), *[
            (Fraction(-1,4),raw3(j,h,h)) for h in range(3) if i == k], *[
            (Fraction(-1,4),raw3(i,h,h)) for h in range(3) if j == k])
    dmk.mul_expr = free_product
    for f in (u2,u3,v2,v3):
        setattr(dmk,f.__name__,lru_cache(None)(f))


def relation20(dmk,i):
    terms = []
    for a,b,c in product(range(3),repeat=3):
        terms.append((2,dmk.mul(dmk.v2(a,c),dmk.v2(b,a),dmk.v3(i,c,b))))
        terms.append((-3,dmk.mul(dmk.v2(i,a),dmk.v2(c,b),dmk.v3(a,b,c))))
    return dmk.add(*terms)


def preimage20(dmk,i):
    # Signed trace words of Eq. (A.2), ordered as printed. P,S are the elementary fields.
    P, S, tr = dmk._ph, dmk._ps, dmk.tr
    terms = []
    for a,b,c,d in product(range(3),repeat=4):
        e = dmk.eps(b,c,d)
        if e:
            for x,y,sign in ((i,a,1),(a,i,-1)):
                terms.append((Fraction(e*sign,2),tr(P(x),S(a),P(y),S(b),S(c),S(d))))
                terms.append((Fraction(-5*e,4),tr(P(x),S(b),S(a),P(y),S(c),S(d))))
        e = dmk.eps(i,c,d)
        if not e:
            continue
        pa,pb,sa,sb,sc,sd = P(a),P(b),S(a),S(b),S(c),S(d)
        terms.extend((Fraction(e*n,den),tr(*w)) for n,den,w in (
            (-17,8,(pa,pb,sb,sc,sa,sd)), (-17,8,(pa,pb,sc,sb,sd,sa)),
            (5,2,(pa,pb,sb,sc,sd,sa)), (1,1,(pa,pb,sc,sa,sb,sd)),
            (15,8,(pa,pb,sc,sb,sa,sd)), (3,8,(pa,sa,sb,sc,pb,sd)),
            (3,8,(pa,sc,sa,sb,pb,sd)), (3,8,(pa,sb,sc,sa,pb,sd))))
        terms.extend((Fraction(e*n,den),dmk.mul(tr(*x),tr(*y))) for n,den,x,y in (
            (-7,16,(pa,pb,sa,sb),(sc,sd)), (-2,1,(pa,pb,sc,sd),(sa,sb)),
            (-1,1,(pa,pb,sc),(sa,sb,sd)), (1,4,(pa,pb,sc),(sa,sd,sb)),
            (5,8,(pa,sa,sc),(pb,sd,sb)), (15,8,(pa,sb,sc),(pb,sd,sa))))
    for a,b in product(range(3),repeat=2):
        for x,u in (((P(a),S(a)),1),((S(a),P(a)),-1)):
            for y,v in (((P(b),S(b)),1),((S(b),P(b)),-1)):
                terms.append((Fraction(5*u*v,4),tr(P(i),*x,*y,{6:1})))
    return dmk.add(*terms)


def summands(dmk, *, primitive):
    symbol = 'r' if primitive else 'R'
    t14 = dmk.tensor(symbol+'14_13')
    t16 = dmk.tensor(symbol+'16_23')
    t20 = lru_cache(None)(lambda i: (preimage20 if primitive else relation20)(dmk,i))
    for label,coefficient,upper,lower in TERMS:
        is_twenty = label >= 37
        terms = []
        for values in product(range(3),repeat=len(upper)):
            hi = [0] * (5 if is_twenty else 9)
            lo = [2] * (3 if is_twenty else 7)
            for i,j,k in zip(upper,lower,values,strict=True):
                hi[int(i)] = lo[int(j)] = k
            if is_twenty:
                terms.append((1,dmk.mul(dmk.v2(hi[1],lo[1]),dmk.v3(hi[2],hi[3],lo[2]),t20(hi[4]))))
                continue
            sign = dmk.eps(*hi[1:4])
            if not sign:
                continue
            first = dmk.v2(hi[4],lo[1])
            second = dmk.v2(hi[5],lo[2])
            if label <= 24:
                third = dmk.v3(hi[6],hi[7],lo[3])
                local = t14(hi[8],*lo[4:7])
            else:
                third = dmk.v2(hi[6],lo[3])
                local = t16(hi[7],hi[8],*lo[4:7])
            terms.append((sign,dmk.mul(first,second,third,local)))
        yield label, coefficient * (-1 if primitive and not is_twenty else 1), dmk.add(*terms)


def load_vendor(vendor):
    vendor = vendor.resolve()
    assert sha256((vendor/'fortuity/bmn/dmk.py').read_bytes()).hexdigest() == '535ec58677878ffe5b613ec93bc4539c7e58eeefc8e532b98e791ccff0fc86b8'
    assert sha256((vendor/'fortuity/bmn/words.py').read_bytes()).hexdigest() == 'f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e'
    sys.path.insert(0,str(vendor))
    from fortuity.bmn import dmk, words
    return dmk,words


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vendor',type=Path,help='directory containing the pinned fortuity/bmn/dmk.py')
    parser.add_argument('mode',choices=('free','source'))
    args = parser.parse_args()
    if not __debug__:
        parser.error('enable assertions')
    dmk,words = load_vendor(args.vendor)
    if args.mode == 'free':
        use_free_gravitons(dmk)
    rows = list(summands(dmk,primitive=args.mode == 'source'))
    result = dmk.add(*[(n,row) for _,n,row in rows])
    separate = sum(len(row) for _,_,row in rows)
    if args.mode == 'free':
        print('free summands, separate monomials, surviving monomials',len(rows),separate,len(result))
        if result:
            print('first nonzero coefficient',min(result.items()))
        assert (len(rows),separate,len(result)) == (23,5369,0)
        assert rows[0][2]
    else:
        here = Path(__file__).resolve().parents[1]
        sys.path.insert(0,str(here))
        from verify import even_part, load_source
        expected = load_source(here/'published_O34_input.json.gz',words)
        even = even_part(result,words)
        saved = json.loads(gzip.decompress((here/'O34_C_even.json.gz').read_bytes()))
        saved = {tuple(map(tuple,w)):Fraction(n) for n,w in saved['terms']}
        print('source terms, matching original, even terms, matching even',
              len(result),result == expected,len(even),even == saved)
        assert result == expected and even == saved


if __name__ == '__main__':
    main()
