# Two scoped perturbative survival certificates in undeformed N=4 SYM

This branch records two finite charge windows. They do not add a nonzero
two-loop lift, a nonperturbative claim, or a priority claim.

## Sp(3), j = 28: ω survives perturbatively to all orders

In the [master report, section 2.6](../../Fortuity_Lift_Master.pdf#page=7),
the incoming d₂ into ω was open. Its complete classical source at
`(7,2,2;1,1), Y=10` has rational dimension at most 454. The Cartan
polynomial computation gives 454 independent highest-weight sources
built from protected multigravitons. Using the quantum product
([arXiv:2611.24160v1, Theorem 9.1](https://arxiv.org/pdf/2611.24160v1#page=86)),
all 454 classes have full perturbative completions; they exhaust this
source. Together with the remaining charge and vanishing-cohomology
checks in the report, this yields survival of ω.

- [Argument and field audit](sp3_omega_proof.md)
- [Full Cartan weight rank data](sp3_weight_ranks.json)
- `python3 replay_sp3_cartan.py` recomputes a weight rank and checks the signed Weyl extraction; `--full` recomputes all twelve weight ranks. No dependencies beyond Python's standard library. This replay verifies the Cartan span calculation, not the independent source upper bound from the report.

## G₂, j = 27: a Cartan-null line survives to E₃

The neutral class `D_+^3(S^3)`, with `S = tr(λ_- λ_+)` and dotted highest weight
`(6,3), Y=6`, is detected by a rational functional that vanishes on the
incoming d₁ and corrected d₂ images. The corrected score vanishes on all
46 rational complements to the pure boundary span and on all 12 mixed
source directions. This result needs no assumption of a quantum Cartan
restriction. It does not determine the other target coordinates or d₃.

The [replay index](g2_j27_replay_index.md), [compact verification manifest](g2_j27_exact_manifest.json),
and [mixed scores](g2_j27_mixed_scores.json) are provided here. The index lists
the underlying computation checkpoints by their original names; the full
matrix cache and large coefficient files are not included in this branch.
