"""Rational product lower bounds and the two colour-parity Euler characteristics for SU(3)."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import random

import mixed_incoming as m
from b2_schur_complex import transpose
from rank_two_bmn_index import quotient
from rank_two_bmn_product_replay import clean, highest, witness
from sp3_omega_product_replay import odd_jet, odd_slot, ordered, trace_factor
from sp3_omega_products import all_products, factor_polynomial
from two_flavour_hooks import schur_colour_dimensions
from unitary_trace_weights import algebra


@lru_cache(None)
def factor(spec, ordinary_diagonal=False, group="SU3"):
    traces = trace_factor(spec, unitary=True)
    length = sum(spec[:5])
    degrees = {"SU3": (2, 3), "G2": (2, 6)}
    if group not in degrees or length not in degrees[group]:
        raise ValueError("expected a primitive scalar degree for the specified group")
    differential = Counter()
    diagonal = Counter()
    for word, c in traces.items():
        for image, n in m.fw.Q_multi(word).items():
            differential[image] += c*n
        tr, = word
        jets = tuple((a, b, sum(1 << i for i in flavour)) for a, b, flavour in tr)
        sign, canonical = ordered(jets, odd_jet)
        diagonal[canonical] += sign*c
    assert not clean(differential)
    assert algebra.conjugate(traces) == {w: (-1)**length*c for w,c in traces.items()}
    assert clean(diagonal) == factor_polynomial(spec)
    result = Counter()
    for word, n in clean(diagonal).items():
        # The Cartan directions have entries (1,z,z²) and (1,z²,z), z³=1, z!=1.
        # Their defining trace is three exactly when the two counts agree modulo three.
        for colours in product(range(2), repeat=length):
            coefficient = (int(not any(colours)) + int(all(colours)) + (-1)**length
                           if ordinary_diagonal else 3*int((length-2*sum(colours)) % 3 == 0))
            # The G2 representation in dimension seven restricts to 3 + 3* + 1.
            if group == "G2":
                coefficient *= 1 + (-1)**length
            if not coefficient:
                continue
            left = tuple(j for j,k in zip(word,colours) if k == 0)
            right = tuple(j for j,k in zip(word,colours) if k == 1)
            crossings = sum(odd_jet(x)*odd_jet(y) for i,x in enumerate(word) if colours[i]
                            for j,y in enumerate(word[i+1:],i+1) if not colours[j])
            result[left,right] += coefficient*int(n)*(-1)**crossings
    return clean(result)


@lru_cache(maxsize=450)
def labelled(factors, ordinary_diagonal=False, group="SU3"):
    if not factors:
        return {((), ()): 1}
    result = Counter()
    for (left, right), old in labelled(factors[:-1],ordinary_diagonal,group).items():
        for (a,b), n in factor(factors[-1],ordinary_diagonal,group).items():
            sl, lword = ordered(left+a, odd_jet)
            sr, rword = ordered(right+b, odd_jet)
            if sl and sr:
                result[lword,rword] += old*n*sl*sr*(-1)**(odd_slot(right)*odd_slot(a))
    return clean(result)


def products_at_length(weight, length, trials, group="SU3"):
    if group not in ("SU3", "G2"):
        raise ValueError("unsupported colour group")
    if group == "G2" and length % 2:
        return dict(length=length, rational_product_rank=0, inputs=[], minor=None)
    # At even length the central inversion in the G2 Weyl group acts trivially.
    maps, _, highest_dim = quotient(weight,length,colour_group="SU3")
    chosen, pivots = [], {}
    if highest_dim:
        maximum = 6 if group == "G2" else 3
        families = all_products(maximum,(0,0,*weight),length,include_odd_powers=group=="SU3")
        rows = [row for family in families.values() for row in family
                if group != "G2" or all(sum(g[:5]) in (2,6) and factor_polynomial(g) for g in row)]
        random.Random(871173+length).shuffle(rows)
        for row in rows[:trials]:
            image = Counter()
            for two_slots,n in labelled(row,group=group).items():
                sign, folded = ordered(two_slots,odd_slot)
                if sign:
                    for j,v in maps.get(folded, ()):
                        image[j] += sign*n*v
            reduced = {j:n % m.P for j,n in image.items() if n % m.P}
            while reduced:
                j = min(reduced)
                if j not in pivots:
                    inverse = pow(reduced[j],-1,m.P)
                    pivots[j] = {k:v*inverse % m.P for k,v in reduced.items()}
                    chosen.append(row)
                    break
                c = reduced[j]
                for k,v in pivots[j].items():
                    residue = (reduced.get(k,0)-c*v) % m.P
                    if residue:
                        reduced[k] = residue
                    else:
                        reduced.pop(k,None)
            if len(chosen) == highest_dim:
                break
    # The modular search only chooses inputs. Each retained rank is certified in Q.
    minor = witness([highest(labelled(row,group=group),weight) for row in chosen]) if chosen else None
    return dict(length=length, rational_product_rank=len(chosen), inputs=chosen, minor=minor)


def calculate(weight, trials):
    if len(weight) != 3 or sorted(weight, reverse=True) != list(weight) or min(weight) < 0:
        raise ValueError("expected a dominant nonnegative BMN weight")
    rows, euler, protected = [], [0,0], [0,0]
    for y in range(max(weight),sum(weight)+1):
        counts = schur_colour_dimensions(transpose(weight), y, 3)
        for i, dim in enumerate(counts):
            euler[i] += (-1)**y*dim
        row = products_at_length(weight, y, trials)
        row["even_odd_cochains"] = counts
        protected[y%2] += row["rational_product_rank"]
        rows.append(row)
        print("length", y, "independent protected products", row["rational_product_rank"], flush=True)
    lower = [max(0,protected[0]-euler[0]), max(0,protected[1]+euler[1])]
    return dict(group="SU(3)",weight=weight,level=2*sum(weight),sample_limit_per_length=trials,
                rows=rows,euler_with_length_sign_by_C=euler,
                independent_products_by_C=protected,opposite_graviton_parity_lower_by_C=lower,
                opposite_graviton_parity_total_lower=sum(lower))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weight",type=int,nargs=3,required=True)
    p.add_argument("--trials",type=int,default=512)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__ or a.trials < 1:
        raise ValueError("expected enabled assertions and positive trial count")
    data=calculate(tuple(a.weight),a.trials)
    a.out.write_text(json.dumps(data,indent=2)+"\n")
