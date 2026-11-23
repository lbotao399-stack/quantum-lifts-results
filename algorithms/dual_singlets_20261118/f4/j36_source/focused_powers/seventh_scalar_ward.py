"""Six-scalar quotient and homogeneous second Ward parameters for the seventh F4 glueball."""

from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PAGE = HERE.parent / "first_page"
sys.path.insert(0, str(PAGE))
from colour import Colour
from intermediate_rank import RHO, ZERO, exchange_trace, orbit, pure_current
from b2_schur_complex import centralizer, character, partitions
from f4_singlet_characters import CARTAN, decomposition, reflect


def schur_decomposition(shape):
    result = Counter()
    for cycles in partitions(sum(shape)):
        coefficient = Fraction(character(shape, cycles), centralizer(cycles))
        if coefficient:
            for highest, copies in decomposition("F4", cycles).items():
                result[highest] += coefficient * copies
    assert all(x.denominator == 1 and x >= 0 for x in result.values())
    return {w: int(n) for w, n in result.items() if n}


def scalar_image():
    native = Colour()
    coroots = []
    for e, f in zip(native.data["raising"], native.data["lowering"]):
        h = native.br[e][f]
        denominator = sum(native.weights[e][i] * n for i, n in h.items())
        coroots.append(tuple(Fraction(2 * h.get(i, 0), denominator) for i in range(4)))

    def labels(w):
        lab = [sum(x*y for x, y in zip(w, h)) for h in coroots]
        assert all(x.denominator == 1 for x in lab)
        return tuple(map(int, lab))

    assert tuple(labels(native.weights[e]) for e in native.data["raising"]) == CARTAN["F4"]
    weight_ranks = {}
    scalar_source_dimension = 0
    for block in json.loads((PAGE / "six_nonzero_weights.json").read_bytes())["blocks"]:
        highest = labels(block["weight"])
        r = block["classical_rank_lower"]
        assert r + block["known_cycle_rank_lower"] == block["columns"]
        assert min(highest) >= 0 and highest != ZERO and highest not in weight_ranks
        points = {highest}
        todo = [highest]
        while todo:
            w = todo.pop()
            for mirror in enumerate(CARTAN["F4"]):
                image = reflect(w, mirror)
                if image not in points:
                    points.add(image)
                    todo.append(image)
        assert len(points) == block["orbit"] and points.isdisjoint(weight_ranks)
        weight_ranks.update(dict.fromkeys(points, r))
        scalar_source_dimension += len(points) * block["columns"]
    assert scalar_source_dimension == 45503928

    def highest_copies(w):
        return sum(s * weight_ranks.get(tuple(x+1-y for x,y in zip(w,rho)), 0)
                   for rho,s in orbit(RHO).items())

    control = json.loads((PAGE / "six_invariants.json").read_bytes())
    from flint import fmpz_mat
    assert control["invariant_dimensions"] == [0, 3, 4]
    assert int(fmpz_mat(control["minor"]).det()) == 292863838464
    weight_ranks[ZERO] = 3 - highest_copies(ZERO)
    assert highest_copies(ZERO) == 3
    targets = schur_decomposition((2, 2, 2))
    highest = {w: highest_copies(w) for w in targets}
    assert all(0 <= n <= targets[w] for w, n in highest.items())
    assert all(highest_copies(w) == 0 for w in weight_ranks
               if min(w) >= 0 and w not in targets)
    return targets, highest, weight_ranks


def source_dimension(shape, dimension):
    result = prod(Fraction(dimension+j-i, shape[i]-j+sum(row > j for row in shape[i+1:]))
                  for i,row in enumerate(shape) for j in range(row))
    assert result.denominator == 1
    return int(result)


def differentiated_scalar_types():
    # Removing the distinguished scalar from the determinant-square flavour singlet
    # leaves Young shape (2,2,1) on the other five colour slots.
    rest = (2,2,1)
    diagrams = []
    for i in range(len(rest)+1):
        shape = list((*rest,0))
        shape[i] += 1
        if shape == sorted(shape,reverse=True):
            diagram = tuple(v for v in shape if v)
            diagrams.append(dict(shape=diagram,dimension=source_dimension(diagram,52),
                                 total_scalar_derivative=diagram == (2,2,2)))
    assert sum(row["dimension"] for row in diagrams)==52*source_dimension(rest,52)
    return diagrams


def run():
    from intermediate_rank import positive_roots, shifted, tensor_multiplicity
    from dotted_filtration import GRAM2, dot
    sys.path.insert(0, str(HERE.parent / "eight_albert_review"))
    from first_page_tower import Equations

    targets, image, weights = scalar_image()
    quotient = {w:n-image[w] for w,n in targets.items() if n>image[w]}
    def dimension(w):
        roots_dual = [tuple(dot(row,r) for row in GRAM2) for r,_ in positive_roots()]
        result = prod(Fraction(dot(shifted(w),r),dot(RHO,r)) for r in roots_dual)
        assert result.denominator == 1
        return int(result)

    target_dimension = source_dimension((2, 2, 2), 52)
    assert sum(dimension(w)*n for w,n in targets.items()) == target_dimension
    assert sum(dimension(w)*n for w,n in image.items()) == sum(weights.values())
    ghosts = sorted(pure_current(5, 0))
    assert len(ghosts) == 2
    blocks = []
    for u, v in combinations_with_replacement(ghosts, 2):
        def copies(w):
            n = tensor_multiplicity(u, v, w)
            if u != v:
                return n
            n += exchange_trace(u, w)
            assert n >= 0 and n % 2 == 0
            return n // 2
        couplings = {w:copies(w) for w in targets}
        blocks.append(dict(first=u, second=v,
                           before_scalar_brackets=sum(n*couplings[w] for w,n in targets.items()),
                           after_scalar_brackets=sum(n*couplings[w] for w,n in quotient.items())))
    assert quotient[ZERO] == 1
    fifth = Equations(5)
    assert fifth.lam == [741,748] and len(fifth.target_rows) == 4
    from flint import fmpz_mat
    assert fmpz_mat(fifth.target_rows).rank() == 1
    even_parameters = sum(b["after_scalar_brackets"] for b in blocks)
    return dict(scalar_cochain_dimension=target_dimension,
                scalar_split_rank=sum(weights.values()), zero_weight_split_rank=weights[ZERO],
                scalar_quotient_dimension=sum(dimension(w)*n for w,n in quotient.items()),
                scalar_quotient=[[w,n] for w,n in sorted(quotient.items())],
                dotted_exchange_even_blocks=blocks,
                homogeneous_parameters_before_mixed_ghost=even_parameters,
                scalar_singlet_mixed_ghost_rows=fifth.target_rows,
                scalar_singlet_ghost_kernel=fifth.lam,
                homogeneous_upper_bound_after_singlet_mixed_ghost=even_parameters-1,
                differentiated_scalar_colour_types=differentiated_scalar_types())


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
