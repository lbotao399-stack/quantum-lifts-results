"""Check the SO(7) target functional and size its incoming correction equations."""
import hashlib
import json
import sys
from itertools import permutations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUNTIME = HERE.parent / "su3_independent_prime_20261118" / "runtime" / "work" / "sp3-native.nosync"
sys.path.insert(0, str(RUNTIME))
from screen_gaugino_source_dimensions_v1 import content_patterns
from exact_dims import _multiplicity


def cochain_dimension(weight, length):
    """Weyl alternation for spin and flavour; the adjoint of SO(7) is exterior square."""
    rho = (2, 1, 0)
    total = 0
    for shift, spin_sign in ((0, 1), (1, -1)):
        for perm in permutations(range(3)):
            inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
            flavour = tuple(weight[2+i] + rho[i] - rho[perm[i]] for i in range(3))
            shifted = (weight[0] + shift, weight[1] - shift) + flavour
            for bosons, fermions, count in content_patterns(shifted, length):
                total += spin_sign * (-1)**inversions * count * _multiplicity(bosons, fermions, 0, 7)
    assert total >= 0
    return total


def check_target():
    digests = {
        "boundaries": "692522aff654bd859c599e7c7bf1fa58f72c783cc367e13891f8f5d7852895ec",
        "SF": "02935df21cba2bf673a8650a2659eaa37f0e3457eb4710c0a2757ccf470e0b39",
        "functional": "92cd8ac99187d25f2378dd724a3d62c5a6ff3936a0aee5226eb64f3348129318",
    }
    prime = 1000003
    with np.load(HERE / "target_frame.npz", allow_pickle=False) as archive:
        assert set(archive.files) == set(digests)
        arrays = {}
        for key, digest in digests.items():
            array = archive[key]
            assert array.dtype == np.dtype("<u4") and np.all(array < prime)
            assert hashlib.sha256(array.tobytes()).hexdigest() == digest
            arrays[key] = array.astype(np.int64)
    boundaries, target, functional = (arrays[k] for k in ("boundaries", "SF", "functional"))
    assert boundaries.shape == (1443, 1870) and target.shape == (1, 1870) and functional.shape == (1870,)
    assert functional.size * (prime - 1)**2 < 2**63
    assert np.all((boundaries @ functional) % prime == 0)
    assert int((target @ functional)[0] % prime) == 1
    return {"prime": prime, "boundary_rows": len(boundaries), "frame_coordinates": len(functional),
            "boundary_entries_checked": int(boundaries.size), "nonzero_boundary_pairings": 0, "SF_pairing": 1}


def run():
    weight = (4, 4, 1, 1, 1)
    target = [cochain_dimension(weight, length) for length in (9, 10, 11)]
    assert target == [4752, 1739, 295]
    incoming = []
    for order in range(1, 5):  # Both derivative weights of the target are four.
        mu = (4-order, 4-order, 1+order, 1+order, 1+order)
        length = 10 + order - 1
        dims = [cochain_dimension(mu, n) for n in (length-1, length, length+1)]
        incoming.append({"order": order, "mu": mu, "length": length, "classical_neighbors": dims})
    return {"group": "SO(7)", "level": 30, "target_mu": weight, "target_length": 10,
            "target_classical_neighbors": target, "modular_target": check_target(), "incoming": incoming}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
