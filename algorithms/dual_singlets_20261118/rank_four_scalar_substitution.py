"""Send the three scalar letters to even linear combinations of the two odd gauginos."""
import argparse
import json
from itertools import product
from pathlib import Path

from flint import fmpq_mat, fmpz_mat

from rank_four_terminal import bounds, quotient
from replay_orthogonal import exchanges


def pure_coefficient(word):
    if any(sum(tr) != 2 for tr in word):
        return 0
    letters = [i for tr in word for i,n in enumerate(tr) for _ in range(n)]
    scalars = [i for i,s in enumerate(letters) if s >= 2]
    assert len(scalars) == 6
    total = 0
    for spins in product(range(2),repeat=6):
        new = letters.copy()
        external = []
        for pos,spin in zip(scalars,spins):
            external.append(2*(letters[pos]-2)+spin)
            new[pos] = spin
        if len(set(external)) != 6 or any(new[i]==new[i+1] for i in range(0,12,2)):
            continue
        # Each auxiliary variable crosses one odd letter at every preceding position.
        flips = sum(scalars)+exchanges(external)+sum(new[::2])
        total += (-1)**flips
    return total


def run(symplectic, orthogonal, output):
    case = bounds(3)
    basis, cov = quotient(case["weight"],case["formal"])
    values = [pure_coefficient(w) for w in basis]
    # Recover the singlet coordinates and check the equality on every weight word.
    pivots = []
    rows = cov.tolist()
    for i,row in enumerate(rows):
        if fmpz_mat([rows[j] for j in pivots]+[row]).rank() > len(pivots):
            pivots.append(i)
    sol = fmpq_mat([rows[j] for j in pivots]).solve(fmpq_mat([[values[j]] for j in pivots]))
    assert (fmpq_mat(cov)*sol).entries() == values
    def test(file):
        data = json.loads(file.read_bytes())
        return [sum(c*values[i] for i,c in w["raw"]) for w in data["relations"]]
    sp, so = test(symplectic), test(orthogonal)
    assert not any(so)
    record = dict(weight=case["weight"], singlet_coordinates=[str(c) for c in sol.entries()],
                  weight_coefficients=[(i,v) for i,v in enumerate(values) if v],
                  symplectic_relations=sp, orthogonal_relations=so,
                  detects_cartan_kernel=any(sp), normalized_pure_target="S^6")
    output.write_text(json.dumps(record,separators=(",",":"))+"\n")
    print(record["singlet_coordinates"],sp,so,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("symplectic","orthogonal","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    assert __debug__
    run(**vars(p.parse_args()))
