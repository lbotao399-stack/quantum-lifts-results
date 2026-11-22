"""Verify the local contractions and closed first pairing into a pure glueball power."""
from collections import Counter
from fractions import Fraction
from math import factorial
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from orthogonal_glueball_first_pair import (
    GHOSTS, PHI, PfaffianCoefficient, highest_product, orientations, single_quantum, traces,
)
from orthogonal_glueball_pair_replay import adjoint_coordinate, cycle_pin
from rank_four_lowered_replay import terminal


def coefficient(colours, power):
    if power < 4:
        raise ValueError('the source requires six single ghost paths and three scalar paths')
    return Fraction(3 * 2**(power-4) * (2*colours-3*power-4), 2*factorial(power-4))


def check_local(colours):
    # An active trace contains at most two ghosts, so one double path suffices
    # to recover its entire polynomial before spectator bilinears are multiplied.
    pf = PfaffianCoefficient(1, 1)
    _, bilinear = traces.so_canon((GHOSTS,))
    spectator = {(9,0,0): 1, (20,0,0): Fraction(1,2),
                 (34,0,0): Fraction(1,2), (0,0,1): -2}
    assert pf.trace(bilinear[0]) == cycle_pin(bilinear[0], 1, pf) == spectator
    odd_part = {k: v for k, v in spectator.items() if k[0]}
    assert pf.multiply(pf.multiply(odd_part, odd_part), odd_part) == {(63,0,0): Fraction(3,2)}
    counts = []
    for size in (1, 2):
        subset = bilinear * size
        orientation, wanted = terminal(subset)
        table = Counter()
        active_traces = orientations(PHI + GHOSTS*(size-1))
        for active in active_traces:
            direct = adjoint_coordinate(colours, active, subset)
            assert direct == orientation * single_quantum(colours, active).get(wanted, 0)
            entries = cycle_pin(active, 1, pf)
            assert entries == pf.trace(active)
            for (mask, markers, doubles), value in entries.items():
                table[mask, doubles] += direct * value * pf.weights.get(markers, 0)
        expected = ({(0,0): -12*(colours-2)} if size == 1 else
                    {(9,0): 24, (20,0): 12, (34,0): 12, (0,1): -36})
        assert {k: v for k, v in table.items() if v} == expected
        counts.append(len(active_traces))
    return counts


def verify():
    # One contraction of two colour generators has at most one empty colour
    # loop when a bilinear survives. Two dimensions thus fix its affine term.
    local = {n: check_local(n) for n in (7, 11)}
    rows = []
    for power in range(4, 15):
        for colours in (power+3, power+4):
            raw = highest_product(colours, power, 1)
            actual = (-1)**power * Fraction(raw['coefficient'])
            assert actual == coefficient(colours, power)
            rows.append([colours, power, str(actual)])
    return dict(local_orientations=local, explicit_pairings=rows)


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('run with assertions enabled')
    print(json.dumps(verify()))
