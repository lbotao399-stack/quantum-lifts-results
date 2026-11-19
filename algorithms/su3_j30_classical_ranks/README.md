# SU(3), level 30: the nine-letter target of R-007

From the repository root run `python3 algorithms/su3_j30_classical_ranks/verify_chains.py` (standard library; assertions enabled). The result must equal [`N4_SU3_j30_verified_chains_20261119.json`](../../data/N4_SU3_j30_verified_chains_20261119.json). This verifies the records and their common ordered coordinates. It does not regenerate the large eliminations or multiply the packed matrices again.

The verifier reads the following committed inputs:

- `algorithms/su3_independent_prime_20261118/SOURCE_MANIFEST.json`, the three production profile plans and all 74 ordered source-recipe lists in that package;
- `joint_plan.json`, `joint.py`, `reduce_joint.py` in this folder;
- `data/N4_SU3_j30_{source_order,incoming_rank,joint_rank}_20261119.json`, `data/N4_SU3_j30_independent_residuals_20261118.json`;
- `data/N4_SU3_j30_{original_stages,target_frames}_p1000033.zip` and `data/N4_SU3_j30_solution_receipts_20261118.zip`.

At prime 1,000,033 the first constraint, the coupled constraint and the incoming classical map have ranks 18,003, 50,147 and 32,144. Subtracting gives zero for the incoming first differential on the recorded first-correction source coordinates. The augmented system is consistent on those coordinates for the recorded source. Byte-identical recipes chosen at another prime need not retain full rank after reduction. Finite-field consistency does not prove rational consistency.

The nine-letter cochain dimension is 63,983. The outgoing classical rank is at least 31,230 over the rationals; stacking Cartan evaluation increases this lower bound to 31,838. The incoming rational rank is at least 32,144 and its image lies in both kernels. Thus the rational cohomology dimension is at most 609 and the dimension invisible on the Cartan is at most one. A matching nonzero class or nonzero higher obstruction does not follow from either bound.

[Independent product replay](../replay_su3_j30_packed_residuals.py) · [source order checks](../su3_j30_source_audit/compare.py) · [target-frame minors](../su3_j30_frame_audit/reproject.py) · [three-prime congruences](../../data/N4_SU3_j30_three_prime_reconstruction_20261119.json).
