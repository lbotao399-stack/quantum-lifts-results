# SU(3) second-correction algorithms, November 18, 2026

This bundle records the ongoing nonabelian, charge-conjugation-odd `j=30` calculation in undeformed N=4 SYM. It contains the second-correction map drivers, their controls, the coupled solver, and the completed negative `Q2` source evaluations. The second coupled equation is not yet solved.

| Space or map | Profiles | Coordinates or couplings |
|---|---:|---:|
| First-correction coordinates | 17 | 28311 |
| First-equation rank, modulo 1000003 | 13 target blocks | 18003 |
| Homogeneous first-correction freedom | | 10308 |
| Second-correction coordinates | 57 | 49328 |
| Second-equation target space | 37 | 63983 |
| Structural Q0 couplings | | 221 |
| Structural Q1 couplings | | 104 |
| Joint variables | | 77639 |

All 94 second-source and second-target profile bases reached their exact character dimensions in production. `basis_inventory.json` verifies the saved seed hashes and counts. It does not independently reconstruct every basis minor entry.

## Coupled algorithm

1. Preserve each original second-target evaluation frame. Add only the source or differential letters absent from that frame.
2. Evaluate all 49328 second-correction coordinates under `Q0`, and all 28311 first-correction coordinates under `Q1`. Keep every structurally allowed profile coupling.
3. Reuse the saved first-equation elimination stages. Move their augmented RHS column from index 28311 to index 77639, retaining their pivot indices and all 10308 homogeneous parameters.
4. Append the second-equation blocks and the already negative `Q2` source. Eliminate in increasing target-dimension order, checkpoint each stage, and back-substitute only after the complete system is assembled.
5. Re-evaluate both original equations on the resulting particular solution. A modular result still requires separate rational certification.

The assembler checks source recipes, source-basis hashes, global row offsets, evaluation-frame equality, complete row coverage, and absence of overlapping blocks. Its dense control includes a system for which fixing the first particular solution falsely produces an obstruction, while retaining its homogeneous freedom solves both equations.

The production Q1 plan has 892 chunks of at most 32 rows. Its existing two-row checkpoints survive interrupted scheduler steps. The solver accepts `--threads` to configure FLINT's internal thread limit; this setting alone is not a measured parallel speedup.

All 804 second-Q0 chunks are complete, covering 49328 source rows. At 2026-11-18 01:24:33 CST, Q1 production had 476 checkpointed rows. The actual saved first stages passed 36006 residual coordinates after moving the RHS column: one saved particular solution and one deterministic homogeneous vector were checked. The remaining free-variable count before imposing the second equation is `77639 - 18003 = 59636 = 49328 + 10308`; see `validation/actual_first_stage_compatibility.json`.

## Reproduce the independent checks

Use Python 3.13 and packages built for that interpreter. Run from this directory:

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python runtime/work/sp3-native.nosync/validate_su3_second_coupled_solve_v1.py
.venv/bin/python runtime/work/sp3-native.nosync/replay_su3_second_rhs_v1.py --evidence rhs_evidence
```

These two commands require no native compilation. The solver control checks a 30-equation, 16-variable system of rank 11, all five homogeneous solutions, an inconsistent control, checkpoint resumption, and rejection of incomplete or incompatible blocks.

The RHS replay verifies 116 input payload hashes and reconstructs all **65406 production coordinates** using the independent Python matrix evaluator. All 23055 source terms lie in target profile 31. Its 480 saved coordinates are nonzero; every other target RHS is zero. This is a replay of the supplied `Q2` image, not a second-obstruction certificate.

To build and check the production map backends:

```sh
.venv/bin/python build_native.py
.venv/bin/python runtime/work/sp3-native.nosync/validate_su3_second_correction_map_v1.py
.venv/bin/python runtime/work/sp3-native.nosync/validate_su3_second_quantum_map_v1.py --basis first_basis
.venv/bin/python runtime/work/sp3-native.nosync/validate_su3_second_rhs_v1.py --source rhs_evidence/quantum_source/Q2.json.gz
```

The Q0 controls cover 906 coordinates at two primes. The Q1 controls cover all 17 source profiles at two primes and 4732 coordinates. The small-frame RHS control covers 902 coordinates. The independent full-frame replay above is stronger for the stored production RHS.

`build_native.py` requires a C++17 compiler. The `.dylib` suffix is used by these drivers on both macOS and Linux; compiled binaries are excluded. `first_basis/` contains the complete first-correction seed recipes needed by the Q1 controls. Generated replay and validation reports are excluded from the frozen source manifest.

## Counts and provenance

Phase increments: `d1 = 0`, audited fortuitous-generator families `= 0`. Cumulative counts: `d1 = 5`, families `= 3`, `d2 = d3 = nonperturbative = 0`. Orders denote supercharge obstructions. The family count does not assert global ring minimality or absolute literature novelty.

The vendored graded trace and matrix-evaluation modules originate from [davidsd/fortuity-ai](https://github.com/davidsd/fortuity-ai), retaining the imported headers and the prior upstream pin `4256f207986029540060fb44dd494746828f65ce`. This bundle adds the profile-resolved quantum maps, preservation of target frames, staged coupled assembly, and independent evidence replay. No new license is assigned to the imported snapshot.

The Q1 production matrices and the coupled second-equation result remain pending. Cluster credentials and account-specific submission scripts are not included.
