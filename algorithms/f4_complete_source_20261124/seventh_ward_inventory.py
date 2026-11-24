"""List every scalar colour type in the complete remaining seventh-power Ward system."""
import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.inputs))
    from seventh_scalar_ward import scalar_image, differentiated_scalar_types
    from intermediate_rank import RHO, exchange_trace, positive_roots, pure_current, shifted, tensor_multiplicity
    from dotted_filtration import GRAM2, dot
    target, image, weights = scalar_image()
    ghosts = sorted(pure_current(5, 0))
    dotted_pairs = list(combinations_with_replacement(ghosts, 2))
    rows = []
    for highest in sorted(target):
        remaining = target[highest] - image[highest]
        if not remaining:
            continue
        dimensions = [Fraction(dot(shifted(highest), tuple(dot(row, r) for row in GRAM2)),
                               dot(RHO, tuple(dot(row, r) for row in GRAM2)))
                      for r, unused in positive_roots()]
        dimension = prod(dimensions)
        assert dimension.denominator == 1
        couplings = []
        for first, second in dotted_pairs:
            multiplicity = tensor_multiplicity(first, second, highest)
            if first == second:
                multiplicity += exchange_trace(first, highest)
                assert multiplicity % 2 == 0
                multiplicity //= 2
            assert multiplicity >= 0
            couplings.append(multiplicity)
        rows.append({"scalar_colour_dynkin_labels": highest,
                     "irreducible_colour_dimension": int(dimension),
                     "terminal_copies": target[highest], "scalar_bracket_image_copies": image[highest],
                     "remaining_scalar_copies": remaining,
                     "invariant_tensor_multiplicities_by_dotted_pair": couplings,
                     "homogeneous_parameters_by_dotted_pair": [remaining * n for n in couplings],
                     "homogeneous_parameters_total": remaining * sum(couplings)})
    totals = [sum(row["homogeneous_parameters_by_dotted_pair"][i] for row in rows) for i in range(3)]
    assert totals == [641, 1614, 1387]
    assert len(rows) == 35
    assert sum(row["irreducible_colour_dimension"] * row["remaining_scalar_copies"] for row in rows) == 88430447
    result = {
        "schema": "f4-seventh-complete-ward-character-inventory-v1",
        "public_commit": "c778cd7c3ceb72ad8082d791b2aee5f0df082164",
        "source_cochain_dimensions_at_lengths_14_15_16": [217297, 62866, 8553],
        "first_correction_cochain_dimension": 23895,
        "first_correction_outgoing_rank": 2959,
        "homogeneous_first_corrections_retained": 20936,
        "dotted_colour_pairs_in_order": dotted_pairs,
        "scalar_colour_types": rows,
        "homogeneous_parameters_after_scalar_equations": sum(totals),
        "scalar_singlet_mixed_ghost_rank": 1,
        "homogeneous_parameter_upper_bound_after_that_equation": sum(totals) - 1,
        "differentiated_scalar_colour_types_before_scalar_bracket_relations": differentiated_scalar_types(),
        "remaining_complete_calculation": "Evaluate the full mixed-ghost equations and both differentiated-scalar colour types (3,2,1) and (2,2,1,1), including every invariant tensor multiplicity listed here and all 20936 first-correction freedoms. Test the complete intersection with the scalar-bracket image. No chosen Albert-word family substitutes for this system.",
        "full_second_map_decided": False,
        "scope": "This complete character inventory refines a public Ward reduction. It supplies no native basis of the listed multiplicity spaces and no rank of the unresolved mixed equations."}
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "scalar_colour_types"}), flush=True)


if __name__ == "__main__":
    main()
