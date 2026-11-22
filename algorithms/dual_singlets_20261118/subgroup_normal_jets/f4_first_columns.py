"""Evaluate the classical and first supercharges on three level-36 F4 columns."""

import json

from colour_fans import colour_fans, f4
from f4_sixth import verify
from orthogonal_glueball_first_pair import GHOSTS, PHI
from fortuity.full.letters import Q_letter
from q1d_exact import pair_terms, Q1d_multi
from mixed_incoming import so_canon


def calculate():
    pairing = verify()["pairing_on_F4_S6"]
    colour = colour_fans("B4 in F4", f4())["one_loop_colour"]
    tangent, transverse = (int(row["coefficient"]) for row in colour)
    tree, intrinsic, full = [], [], []
    for i, scalar in enumerate(PHI):
        psi = (0, 0, tuple(j for j in range(3) if j != i))
        # All spectator gauginos are tree closed and cannot form a one-loop pair.
        assert not Q_letter(scalar)
        assert all(not Q_letter(g) and not pair_terms(psi,g) and not pair_terms(scalar,g)
                   for g in GHOSTS)
        assert not any(pair_terms(x,y) for x in GHOSTS for y in GHOSTS)
        left, right = [p for j,p in enumerate(PHI) if j != i]
        assert Q_letter(psi) == ((1,left,right),(-1,right,left))
        permutation = [PHI.index(p) for p in (left,right,scalar)]
        tree.append((-1)**sum(x>y for j,x in enumerate(permutation) for y in permutation[j+1:]))

        terms = {(x,y):c for c,x,y in pair_terms(psi,scalar)}
        # The native pair uses [t,P] in its first slot. Both dotted
        # orders are present, with a second minus from the odd exchange.
        factor = int(terms[GHOSTS[::-1]] - terms[GHOSTS])
        intrinsic.append(tangent*factor)
        full.append((tangent+transverse)*factor)
        source = ((psi, scalar),) + (GHOSTS,)*5
        sign, word = so_canon((GHOSTS,)*6)
        assert Q1d_multi(source,9,"SO_odd") == {word: sign*intrinsic[-1]}

    assert tree == [1,-1,1] and full == [36*k for k in tree]
    return dict(tree_coefficients=tree, spin9_first_coefficients=intrinsic,
                F4_first_coefficients=full, first_pairings=[pairing*k for k in full],
                closure_kernel_dimension=2, first_rank_on_that_kernel=0)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    print(json.dumps(calculate(), indent=2))
