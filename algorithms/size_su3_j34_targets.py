"""Size the full-letter targets of the published (6; 8,5,4) SU(3) BMN sector.

Use the existing untwisted and charge-conjugation-twisted character routines.
Lengths in each row are the domain of Q0, the target of d_r, and its Q0 image
space. All entries are dimensions before taking cohomology.
"""
import json
import sys
from pathlib import Path


runtime = (Path(__file__).resolve().parent / "su3_independent_prime_20261118"
           / "runtime" / "work" / "sp3-native.nosync")
sys.path.insert(0, str(runtime))

from screen_gaugino_source_dimensions_v1 import highest_dimension
from su3_charge_dimensions_v1 import controls, dimensions


def cochains(weight, length):
    even = sum(dimensions(weight, length, 1).values())
    odd = sum(dimensions(weight, length, -1).values())
    full = highest_dimension(weight, length, "SU3")
    assert even + odd == full
    return {"length": length, "C_even": even, "C_odd": odd, "both": full}


def run():
    result = {"charge_conjugation_controls": controls()}
    source = cochains((0, 0, 8, 5, 4), 11)
    assert source["both"] == 26589  # Public unrefined BMN source multiplicity.
    result["j34_source"] = source
    # In the published (8,5,4) column there are 37 multiplets, 34 gravitons,
    # and the P-even minus P-odd count is 31. All gravitons have P=+1.
    parity_odd = (37 - 31) // 2
    assert parity_odd == 37 - 34 == 3
    # P=C*(-1)**Y, where the two nongravitons at spin label 6 have Y=11.
    source_c = (-1) * (-1)**source["length"]
    result["public_column"] = {"multiplets": 37, "gravitons": 34, "P_odd": parity_odd,
                               "source_C": source_c, "scalar_restriction_kernel_dimension": 2 - 1,
                               "source": "davidsd/fortuity-ai@4256f207986029540060fb44dd494746828f65ce, summary/fortuity-summary.tex and data/su3/core_restrictions_j30_j34.json"}
    control = cochains((2, 2, 5, 2, 2), 9)
    assert control["C_odd"] == 63983  # Audited j=30 second-equation target.
    result["j30_second_target_control"] = control
    targets = []
    for order in (1, 2, 3):
        weight = (order, order, 8 - order, 5 - order, 4 - order)
        target_length = 12 - order
        neighbors = [cochains(weight, length) for length in range(target_length - 1, target_length + 2)]
        targets.append({
            "order": order, "mu": list(weight),
            "classical_neighbors": neighbors,
            "target_P": source_c * (-1)**target_length,
            "dense_C_even_incoming_Q0_bytes_at_four_bytes_per_entry":
                4 * neighbors[0]["C_even"] * neighbors[1]["C_even"],
        })
    result["j34_targets"] = targets
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
