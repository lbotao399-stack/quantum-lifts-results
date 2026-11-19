# G2 j=27 neutral Cartan-null class: exact survival to E3

**Proved:** the nonzero class of D_+^3(S^3), in the dotted highest
weight neutral target theta000, n=(6,3), Y=6, survives to E3.
One exact functional ell on the target pairs to 1 with that class
and annihilates all incoming d1 and d2. Neither outgoing differential
can occur from the neutral theta degree. This survival argument does
not use a quantum Cartan-restriction assumption. In particular it does
not claim that the other two classical neutral coordinates have
zero incoming d2. No conclusion at E4 or in other charges.

Witness and replay order (run from project root, one BLAS thread):

| Check | Coefficient or rank witness | Replay |
| --- | --- | --- |
| Target nonzero and target Q0 dual | j27_n63_C6.json, j27_n63_target.json, j27_n63_derivative.json, j27_n63_d1_dual.json | j27_n63_C6.py, j27_n63_target.py, j27_n63_derivative.py, j27_n63_d1_dual.py |
| ell kills d1, exact factor ell Q1 = mu Q0 | j27_n63_d1_factor.json.gz; all 2621 original columns checked | j27_n63_d1_factor.py |
| Actual first-order compensation exists on closed sources | j27_n63_intermediary_control.json (ranks 1530+1133=2663) | j27_n63_intermediary_verify.py |
| Source normalization (classical boundaries and protected products) | j27_n63_k4_reduced.json and two exact head banks | j27_n63_k4_reduced.py; J27_N63_NEUTRAL_C6.md |
| Complete raw corrected score rows | j27_n63_chi4_full_exact.json.gz and j27_n63_chi5_full_exact.json.gz | j27_n63_full_exact_scores.py 4, then 5 |
| Literal incoming boundary span (rank 2536) | j27_n63_pure_incoming_rank.json | j27_n63_pure_incoming_rank.py |
| Pure kernel, remaining 46 independent exact complements | j27_n63_pure_quotient_mod.json and j27_n63_pure_quotient_exact.json | j27_n63_pure_quotient_mod.py, j27_n63_pure_quotient_exact.py |
| Mixed kernel, twelve independent exact closed lifts and scores | j27_n63_k6_allheads_exact.json and j27_n63_chi_allheads_exact.json | j27_n63_k6_allheads_exact.py, j27_n63_chi_allheads_exact.py |
| Assemble exact dimension sandwich and verify coefficient readout | j27_n63_full_rank_verified_exact_zero.json | j27_n63_finalize_pure_quotient.py |

The replay hashes cover exact score and source coefficients, normalization
banks, chain matrices and pivot-selection data, and the compact outcome.
No caches are duplicated. For full construction details see
J27_N63_NEUTRAL_C6.md and J27_N63_TERMINAL_SCOPE.md.
