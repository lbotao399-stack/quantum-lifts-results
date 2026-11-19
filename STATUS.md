# Running computations: daily status

## Four-source SU(3) test and G2 sign correction, November 18, 2026

At SU(3) level 34 a 365-dimensional quotient is exhausted on the 7,000-coordinate frame. The resulting first correction fails on 120 of 120 fresh coordinates. On 64 noncommuting even jets, all 1,485 cubic coordinates vanish upon restriction to kernels annihilating the projected classical images. These tests leave the source’s first obstruction open. [Data, determinant and replay](algorithms/su3_j34_20261118/README.md#four-source-correction-and-closed-jets).

For the $G_2$ incoming candidates the three-fermion coefficient is $-16/3$ in the archive, versus $+16/3$ in the literal printed identification with both flavour orientations positive. The previous conclusion of protection is withdrawn pending an independent physical-current identification. [Classical ranks and both tests](algorithms/G2_protected_incoming_currents_20261118.md).

Increment: **0** first-order pairings and **0** fortuitous-generator families; totals: **6** and **3**. Certified second-order, third-order and nonperturbative lifts: **0, 0, 0**.


## Mixed singlets and the SU(3) level-34 profile-10 span, November 18, 2026

SO(7) and USp(6) have respectively 3 and 2 classical singlets at bidegree (3,3), length 8. The two formal classes survive in both groups. Appending the known first-order source $F$ completes the SO(7) basis. All five terminal singlet dimensions have also been fixed in both groups. [Rational rank argument and trace identities](algorithms/dual_singlets_20261118/README.md#terminal-singlets-and-the-mixed-length-eight-classes).

For the SU(3) first correction, the four sources mapping to target profile 10 have total dimension 69,657. A functional separates the right-hand side from all 2,752 archived recipes on 2,800 coordinates, but fails on additional words from each of the four sources. [Ranks, explicit counterexamples and coupled target supports](algorithms/su3_j34_20261118/README.md#coupled-correction-and-the-profile-10-separator).

New $d_1$ pairings: **0** (cumulative **6**). New fortuitous-generator families: **0** (cumulative **3**). Certified $d_2$, $d_3$ and nonperturbative lifts: **0, 0, 0**.

## SU(3) joint first correction and dual-group scalar singlets, November 18, 2026

At level 34 the 7,211 profile-16 coefficients pass 7,272 production and 432 fresh coordinates in profile 16, but fail on all 84,462 tested coordinates on profiles 13–15. The source's first differential remains open. [Arrays and replay](algorithms/su3_j34_20261118/README.md#fixed-first-correction-across-profiles).

At zero derivative degree, the length-11 and length-12 singlet cohomology dimensions are 0 and 1 for each of SO(7) and USp(6). The four independently replayed minors and the nonzero Cartan polynomial certify these statements over the rationals. [Ranks, integer signs and Cartan witness](algorithms/dual_singlets_20261118/README.md#zero-derivative-degree).

Increment: **0** first-order pairings and **0** generator families; cumulative **6** and **3**. Certified second-order, third-order and nonperturbative lifts: **0, 0, 0**.

## SO(7) level-30 pairing and dual-group singlet comparison, November 18, 2026

The new first differential onto the line of $SF$ adds one highest-weight pairing: the rank increases from 12,525 to 12,526. A8 job `5628225` reconstructed all 156,900,676 entries of the augmented minor and reproduced each of its 212 panel determinants. Its determinant is 942,245 modulo 1,000,003. The target is absent from the second and later pages. [Proof](algorithms/so7_cartan_target_20261118/README.md) · [Replay receipt](algorithms/so7_cartan_target_20261118/incoming_replay.json).

At doubled charges (0,0,8,8,8) and derivative bidegree (4,4), exact rational computation finds one length-eight SO(7) class and none for USp(6). The SO(7) class is $S^4$, already hit at first order. [Complete classical matrices](algorithms/dual_singlets_20261118/README.md). The $G_2$ incoming classical dimensions are certified; their proposed exhaustion by physically protected current products requires a component identification. [Current calculation](algorithms/G2_protected_incoming_currents_20261118.md).

At SU(3) level 34 the first right-hand side contains 517,468 rational terms. In eight profiles of the third target the projected closed quotient vanishes: 1,263 minus 1,206 minus 57 equals zero. The four minors were independently replayed at two primes. No claim is made for the remaining profiles. At level 30 both finite-field solutions now have zero residual in all 100,307 coordinates; rational solvability remains open. [Level-34 calculation](algorithms/su3_j34_20261118/README.md) · [Level-30 receipts](algorithms/verify_su3_j30_solution_receipts.py).

Increment: **1** first-order pairing, **0** generator families. Cumulative: **6** and **3**, respectively. Certified second-order, third-order and nonperturbative lifts: **0, 0, 0**. [Current ledger](data/N4_lift_campaign_counts_20261118.json).

## SU(3) level-34 direct action and level-30 matrix replay, November 18, 2026

Both A8 jobs completed. Job `5622046` computed the level-34 direct $Q_2$ action on eight CPUs in 12 seconds. A local one-process replay reproduced all 140,526 rational coefficients and the zero Cartan polynomial. Job `5622107` used eight CPUs for 385 seconds to multiply the level-30 solution by all 378 saved matrices. Each of its 100,307 modular residuals vanished. [Level-34 computation](algorithms/su3_j34_20261118/README.md) · [Level-30 replay program](algorithms/replay_su3_j30_packed_residuals.py) · [Matrix hashes and residuals](data/N4_SU3_j30_independent_residuals_20261118.json)

This stage adds 0 $d_1$ pairings and 0 fortuitous-generator families. Campaign totals are 5 highest-weight $d_1$ pairings and 3 families modulo the audited products and descendants. Certified $d_2$, $d_3$ and nonperturbative lifts in undeformed $\mathcal N=4$ SYM remain 0, 0 and 0. [Counts and scope](data/N4_lift_campaign_counts_20261118.json)

## Method review and completed cluster calculation, November 18, 2026

At prime 1,000,033, all 24 A8 map shards finished. The joint solve used 64 CPUs for 4,227 seconds. The production residual was zero in all 34,901 first-equation and 65,406 second-equation coordinates. A6 reached the same 50 ordered pivot lists at prime 1,000,003 and stopped at its application time limit before saving its joint solution. Rational certification and second-prime frame injectivity are open. [Saved coordinates and receipt verifier](algorithms/verify_su3_j30_solution_receipts.py).

The $SO(7)$ test can use a target functional pairing to one with $SF$ and to zero with all 1,443 saved boundaries. The $SU(3)$ level-34 parity calculation and the rational $G_2$ partner replay are included with the review. [Method review and action plan](algorithms/N4_lift_method_review_20261118.md) · [Numerical data](data/N4_lift_method_review_20261118.json).

This phase adds zero $d_1$ highest-weight pairings and zero fortuitous-generator families. The cumulative counts are 5 and 3. Verified $d_2$, $d_3$ and nonperturbative lifts in undeformed $\mathcal N=4$ are 0, 0 and 0.

## Earlier BMN phase, November 18, 2026

[Algorithm and complete replay bundle](algorithms/bmn_rank_certificates_20261117/README.md); [numerical table](data/N4_Sp3_G2_BMN_j30_direct_index_20261117.csv).

- Sp(3)=USp(6) and G2 at j=30: 41 direct highest-weight product sectors per group, 12 weight columns. None of the index-minus-certified-product-rank bounds is negative. Cancelling classical classes remain unresolved by this screen.
- Independent replay: 209 payload hashes; 62 nonzero determinants; all 121237 minor entries reconstructed; 7168 additional matrix-evaluation coordinates; 185 exact generator-closure checks; 141 modular raising rows; rational G2 torus certificate.
- The completed nonabelian G2 `(5,5,5)` column has cohomology dimensions `(0,0,0,0,0,4,0,0,0,0)` for spin labels 0 through 9. At spin 5, `3030-1960-1066=4` and the known-product quotient has dimension `1070-1066=4`: known closed products exhaust all four classes. Its replay verifies 13 saved determinants, 27 rational projected recipes, 1308 independent source-probe coordinates, 185 exact generator closures, 97 modular highest-weight product rows, and all 246078 generic product coordinates. Source minor entries beyond the saved probes remain production evaluations.
- Increment in audited d1 pairings: 0; cumulative: 5. Increment in audited fortuitous-generator families: 0; cumulative: 3. Counts are scoped to the current undeformed N=4 campaign; families are counted modulo the audited products and descendants, without a global minimal-ring or absolute literature-novelty claim.
- Certified d2 / d3 / nonperturbative lifts in undeformed N=4: 0 / 0 / 0. Orders denote supercharge obstructions.
- The full nonabelian SU(3) j=30 second obstruction remains under calculation. All 94 second-source and second-target profile bases reached their required dimensions: 49328 and 63983. The independent negative-Q2 replay verifies 116 payload hashes and all 65406 production coordinates. All 23055 input terms lie in profile 31; its 480 coordinates are nonzero. This does not establish an obstruction.
- The [SU3 algorithm bundle](algorithms/su3_second_correction_20261118/README.md) retains all 10308 homogeneous first-correction parameters. Q0 and Q1 maps are complete, with 49328 and 28311 rows. All 37 second-equation elimination stages report consistency modulo 1000003: rank 50147 for 77639 variables. The particular-solution residual replay and rational certification remain pending.
- [Independent-prime algorithms and data](algorithms/su3_independent_prime_20261118/README.md): at 2026-11-18 08:09:54 CST, all 24 shards run on A8 with 32 CPUs each, reusing the complete rational trace recipes at prime 1000033. 322 of 2197 tasks are complete, with no recorded chunk errors. Lossless matrix packing reduces repeated checkpoint-file reads; controls preserve both assembled equations and reject missing or corrupted input. [Numerical snapshot](data/N4_SU3_independent_prime_20261118.json).

## Historical snapshot, November 12, 2026

Updated 2026-11-12 from the cluster harvests (A8 2026-11-12 22:33, A6 2026-11-12 22:34). Ranks are modulo p = 2^31 - 1 unless stated; blocks are labelled as in the master report. Nonzero differentials listed here are raw scan output: not yet rechecked at a second prime and not yet separated into primary classes and descendants. Full records: [data/daily/harvest_2026-11-12.json](data/daily/harvest_2026-11-12.json).

## N=4 SU(4), all letters, levels 28-30 (trace-relation count of fortuitous classes)

| level | shapes | done | positive fortuitous count | open, with an earlier failed attempt |
|---|---|---|---|---|
| 28 | 3 | 3 | none | none |
| 29 | 29 | 15 | none | (3,2,2,1,1,1,1,1), (4,3,3,1,1,1), (4,3,3,2,1), (4,4,2,1,1,1), (5,3,2,1,1,1), (5,3,2,2,1), (5,4,4,1), (5,5,3,1), (6,4,3,1), (6,5,2,1), (7,3,3,1), (7,4,2,1), (8,4,1,1), (9,2,2,1) |
| 30 | 53 | 0 | none | (3,3,3,2,1,1), (4,4,4,1,1) |

(Level 28 counts only the three shapes reopened in this campaign; the other 35 were done before.)

## N=1 abelian theories with a neutral singlet: four-loop search, weights 18-24

- **massive_singlet**: 296 blocks computed, 16 unfinished; d2 tested 2 (nonzero targets 2)
  - d2 rank 2, massive_singlet [20,2,2,-8,[0],[0]]
  - d2 rank 2, massive_singlet [24,2,2,-10,[0],[0]]

## Rerun of blocks left unfinished by earlier tiers

- 91 blocks computed, 94 unfinished; d3 tested 85 (nonzero targets 68), d4 tested 81 (nonzero targets 57)
  - d3 rank 1, singlet [18,2,2,-9,[0],[0,-3]]
  - d3 rank 1, singlet [19,2,1,-10,[0],[0,-2]]
  - d3 rank 1, singlet_cubic_self [16,1,1,-8,[0],[0]]
  - d4 rank 1, singlet_cubic_self [16,1,1,-8,[0],[0]]
  - d2 rank 1, singlet_cubic_self [15,2,1,-7,[0],[0]]
  - d3 rank 2, singlet_cubic_self [15,2,1,-7,[0],[0]]
  - d2 rank 4, n2_quiver2 [10,2,1,-3,[0,0],[1,-1,-1,-1]]
  - d2 rank 4, n2_quiver2 [10,2,1,-3,[0,0],[-1,1,1,0]]
  - d3 rank 2, massive_singlet [22,2,1,-8,[0],[0]]
  - d2 rank 2, n2_quiver2 [11,2,1,-4,[0,0],[1,-1,-1,0]]
  - d2 rank 2, massive_singlet [22,1,1,-9,[0],[0]]
  - d3 rank 1, singlet [18,1,1,-10,[0],[0,0]]
  - ... 4 more

## N=2 SU(2) and SU(2) + flip singlet, weights 15-16 (components)

- **n2su2_w15**: 16 blocks with classes
  - d2 rank 1, n2su2 [15,0,0,-8,0]
  - d2 rank 1, n2su2 [15,1,0,-8,0]
  - d2 rank 1, n2su2 [15,1,1,-8,0]
- **n2su2_w16**: 22 blocks with classes
  - d2 rank 1, n2su2 [16,1,1,-7,1]
- **su2flip_w15**: 19 blocks with classes
  - d2 rank 1, su2flip [15,0,0,-7,3]
  - d2 rank 1, su2flip [15,1,0,-7,3]
- **su2flip_w16**: 20 blocks with classes

## N=2 SU(3), weights 8-12

- **n2su3_w10**: 13 blocks with classes
- **n2su3_w11**: 15 blocks with classes
  - d2 rank 1, n2su3 [11,0,0,-5,-1]
- **n2su3_w12**: 2 blocks with classes
  - d2 rank 1, n2su3 [12,0,0,-6,0]
- **n2su3_w8**: 20 blocks with classes
- **n2su3_w9**: 19 blocks with classes

## N=2 U(1)^2 quiver: attribution of the two-loop lifts at (11;1,1;-4)

- k1: E1 15, ker d1 1, d2 1, carried by spectator products 0, by products and triples 1, with towers 1
- k2: E1 15, ker d1 1, d2 1, carried by spectator products 0, by products and triples 1, with towers 1

## N=4 SU(3) Leigh-Strassler, levels 23-25

- level 23: running
- level 24: running
- level 25: file 11-06 07:24; failed sectors [[8, 0, 2, 1], [8, 0, 3, 0], [11, 0, 1, 0]]; nonzero d_r [[1, [5, 0, 3, 2], 6, 5], [1, [5, 0, 4, 1], 6, 5], [1, [5, 0, 5, 0], 6, 1]]

## Jobs

| cluster | job | name | state | time | cores |
|---|---|---|---|---|---|
| A8 | 5478222_1 | finretry | RUNNING | 4-00:22:42 | 128 |
| A8 | 5478222_0 | finretry | RUNNING | 4-00:25:59 | 128 |
| A8 | 5478224 | finsu2 | RUNNING | 4-00:31:43 | 40 |
