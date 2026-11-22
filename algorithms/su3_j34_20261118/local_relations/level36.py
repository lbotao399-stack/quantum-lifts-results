"""Use the repaired R20 in the independent [4,1] bosonic seed of Kim et al."""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys

from syzygy import load_vendor, preimage20, relation20, use_free_gravitons

HERE = Path(__file__).resolve().parent


def rebuild(dmk,primitive):
    names = ('12_03','12_11','12_22','14_13','16_23')
    t03,t11,t22,t14,t16 = [dmk.tensor(('r' if primitive else 'R')+s) for s in names]
    t20 = lru_cache(None)(lambda i: (preimage20 if primitive else relation20)(dmk,i))
    terms = []
    sizes = []
    records = json.loads((HERE/'level36_terms.json').read_bytes())
    assert len(records) == len({r[0] for r in records}) == 52
    for label,coefficient,upper,lower in records:
        contracts = []
        for values in product(range(3),repeat=len(upper)):
            hi,lo = [0]*11,[2]*8
            for i,j,k in zip(upper,lower,values,strict=True):
                hi[i] = lo[j] = k
            sign = 1
            if label < 107:
                sign = dmk.eps(*hi[1:4])
                if label <= 2:
                    sign *= dmk.eps(*hi[4:7])
                if not sign:
                    continue
            if label <= 10:
                start = 7 if label <= 2 else 4
                factors = [dmk.v2(hi[start+j],lo[1+j]) for j in range(4)]
                if label <= 2:
                    local = t03(*lo[5:8])
                elif label <= 4:
                    local = t11(hi[8],lo[5])
                else:
                    local = t22(hi[8],hi[9],lo[5],lo[6])
            elif 26 <= label <= 57:
                factors = [dmk.v2(hi[4],lo[1]),dmk.v3(hi[5],hi[6],lo[2]),dmk.v3(hi[7],hi[8],lo[3])]
                local = t14(hi[9],*lo[4:7])
            elif 82 <= label <= 106:
                factors = [dmk.v2(hi[4],lo[1]),dmk.v2(hi[5],lo[2]),dmk.v3(hi[6],hi[7],lo[3])]
                local = t16(hi[8],hi[9],*lo[4:7])
            elif 107 <= label <= 109:
                factors = [dmk.u2(hi[1],hi[2]),dmk.v2(hi[3],lo[1]),dmk.v2(hi[4],lo[2])]
                local = t20(hi[5])
            elif 110 <= label <= 111:
                factors = [dmk.v3(hi[1],hi[2],lo[1]),dmk.v3(hi[3],hi[4],lo[2])]
                local = t20(hi[5])
            else:
                raise ValueError(label)
            contracts.append((sign,dmk.mul(*factors,local)))
        row = dmk.add(*contracts)
        sizes.append(len(row))
        terms.append((coefficient * (-1 if primitive and 26 <= label <= 106 else 1),row))
    return dmk.add(*terms),sum(sizes)


def check_mixed(dmk,words,result):
    sys.path.insert(0,str(HERE.parent))
    from verify import (add_scaled, boundary_basis, even_part, load_source, pairing, polynomial,
                        reduce_row, restrict_block)
    from match_threshold import paper_core
    # Four abelian scalars and one abelian fermion have C=-1 on the SU(2) block.
    projected = dict(dmk.add((1,result),(-1,even_part(result,words))))
    mixed = restrict_block(projected,(4,0,0),words,5)
    assert mixed == restrict_block(result,(4,0,0),words,5)
    record = json.loads((HERE.parent/'certificate.json').read_bytes())['channels'][2]
    ell = {(e,o):Fraction(c) for e,o,c in record['witness']['monomial_coefficients']}
    core = polynomial(paper_core(words))
    source = polynomial(mixed)
    numerator,denominator = pairing(ell,source),pairing(ell,core)
    assert (numerator,denominator) == (19780588800,-20)
    predecessors = words.words_of_sector((6,4,4,4),3)
    basis = boundary_basis(predecessors,words)
    assert all(not pairing(ell,row) for _,row,_ in basis) and len(predecessors) == 197
    multiple = numerator/denominator
    difference = dict(source)
    add_scaled(difference,core,-multiple)
    remaining,used = reduce_row(difference,basis)
    primitive = {predecessors[i]:c for i,c in used.items()}
    assert not remaining and polynomial(words.Q_expr(primitive)) == difference
    assert multiple*35 == -34616030400

    # In multiplying either seed by Tr((phi^1)^2), only two higher abelian
    # coefficients can feed the same SU(2) threshold. Test both before any iteration.
    seed34 = even_part(load_source(HERE.parent/'published_O34_input.json.gz',words),words)
    higher = {}
    for step in (2,4):
        predecessors = words.words_of_sector((6,4-step,4,4),3)
        basis = boundary_basis(predecessors,words)
        for n,seed in ((3,seed34),(4,projected)):
            target = polynomial(restrict_block(seed,(n+step,0,0),words,5))
            remaining,used = reduce_row(target,basis)
            primitive_high = {predecessors[i]:c for i,c in used.items()}
            assert not remaining and polynomial(words.Q_expr(primitive_high)) == target
            higher[f'{n}+{step}'] = dict(monomials=len(target),primitive_terms=len(primitive_high))
    return dict(paper_terms=len(result),odd_terms=len(projected),mixed_terms=len(mixed),
                pairing=str(numerator),core_pairing=str(denominator),mixed_coefficient=str(multiple),
                primitive_terms=len(primitive),higher_scalar_coefficients=higher)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('vendor',type=Path)
    parser.add_argument('mode',choices=('free','source'))
    args = parser.parse_args()
    if not __debug__:
        parser.error('enable assertions')
    dmk,words = load_vendor(args.vendor)
    if args.mode == 'free':
        use_free_gravitons(dmk)
    result,separate = rebuild(dmk,args.mode == 'source')
    print('52 summands, separate monomials, surviving monomials',separate,len(result),flush=True)
    if args.mode == 'free':
        assert separate == 15259 and not result
    else:
        assert separate == 83109 and len(result) == 11686
        print(json.dumps(check_mixed(dmk,words,result),indent=2))
