# Running computations: daily status

Updated 2026-11-10 from the cluster harvests (A8 2026-11-10 22:45, A6 2026-11-10 22:43). Ranks are modulo p = 2^31 - 1 unless stated; blocks are labelled as in the master report. Nonzero differentials listed here are raw scan output: not yet rechecked at a second prime and not yet separated into primary classes and descendants. Full records: [data/daily/harvest_2026-11-10.json](data/daily/harvest_2026-11-10.json).

## N=4 SU(4), all letters, levels 28-30 (trace-relation count of fortuitous classes)

| level | shapes | done | positive fortuitous count | open, with an earlier failed attempt |
|---|---|---|---|---|
| 28 | 3 | 2 | none | (5,4,2,1,1) |
| 29 | 29 | 15 | none | (3,2,2,1,1,1,1,1), (4,3,3,1,1,1), (4,3,3,2,1), (4,4,2,1,1,1), (5,3,2,1,1,1), (5,3,2,2,1), (5,4,4,1), (5,5,3,1), (6,4,3,1), (6,5,2,1), (7,3,3,1), (7,4,2,1), (8,4,1,1), (9,2,2,1) |
| 30 | 53 | 0 | none | none |

(Level 28 counts only the three shapes reopened in this campaign; the other 35 were done before.)

## N=1 abelian theories with a neutral singlet: four-loop search, weights 18-24

- **massive_singlet**: 296 blocks computed, 16 unfinished; d2 tested 2 (nonzero targets 2)
  - d2 rank 2, massive_singlet [20,2,2,-8,[0],[0]]
  - d2 rank 2, massive_singlet [24,2,2,-10,[0],[0]]

## Rerun of blocks left unfinished by earlier tiers

- 86 blocks computed, 92 unfinished; d3 tested 80 (nonzero targets 63), d4 tested 77 (nonzero targets 53)
  - d3 rank 1, singlet [18,2,2,-9,[0],[0,-3]]
  - d3 rank 1, singlet [19,2,1,-10,[0],[0,-2]]
  - d3 rank 1, singlet_cubic_self [16,1,1,-8,[0],[0]]
  - d4 rank 1, singlet_cubic_self [16,1,1,-8,[0],[0]]
  - d2 rank 1, singlet_cubic_self [15,2,1,-7,[0],[0]]
  - d3 rank 2, singlet_cubic_self [15,2,1,-7,[0],[0]]
  - d2 rank 4, n2_quiver2 [10,2,1,-3,[0,0],[1,-1,-1,-1]]
  - d2 rank 4, n2_quiver2 [10,2,1,-3,[0,0],[-1,1,1,0]]
  - d2 rank 2, n2_quiver2 [11,2,1,-4,[0,0],[1,-1,-1,0]]
  - d2 rank 2, massive_singlet [22,1,1,-9,[0],[0]]
  - d3 rank 1, singlet [18,1,1,-10,[0],[0,0]]
  - d3 rank 2, singlet_cubic_self [20,1,1,-11,[0],[0]]
  - ... 1 more

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
- **su2flip_w16**: 19 blocks with classes

## N=2 SU(3), weights 8-12

- **n2su3_w10**: 13 blocks with classes
- **n2su3_w11**: 15 blocks with classes
  - d2 rank 1, n2su3 [11,0,0,-5,-1]
- **n2su3_w12**: 2 blocks with classes
  - d2 rank 1, n2su3 [12,0,0,-6,0]
- **n2su3_w8**: 20 blocks with classes
- **n2su3_w9**: 19 blocks with classes

## N=2 U(1)^2 quiver: attribution of the two-loop lifts at (11;1,1;-4)

- running

## N=4 SU(3) Leigh-Strassler, levels 23-25

- level 23: running
- level 24: running
- level 25: file 11-06 07:24; failed sectors [[8, 0, 2, 1], [8, 0, 3, 0], [11, 0, 1, 0]]; nonzero d_r [[1, [5, 0, 3, 2], 6, 5], [1, [5, 0, 4, 1], 6, 5], [1, [5, 0, 5, 0], 6, 1]]

## Jobs

| cluster | job | name | state | time | cores |
|---|---|---|---|---|---|
| A8 | 5497036 | r12_4_5,4,2,1,1 | PENDING | 0:00 | 128 |
| A8 | 5497046 | fk4_30_0 | PENDING | 0:00 | 128 |
| A8 | 5497047 | fk4_30_1 | PENDING | 0:00 | 128 |
| A8 | 5497049 | fk4_30_2 | PENDING | 0:00 | 128 |
| A8 | 5497050 | fk4_30_3 | PENDING | 0:00 | 128 |
| A8 | 5497051 | fk4_30_4 | PENDING | 0:00 | 128 |
| A8 | 5497053 | fk4_30_5 | PENDING | 0:00 | 128 |
| A8 | 5478228 | ls25c | RUNNING | 1-23:31:22 | 128 |
| A8 | 5497044 | fk4_29_3 | RUNNING | 1-05:31:25 | 64 |
| A8 | 5478229 | ls24c | RUNNING | 1-23:29:10 | 128 |
| A8 | 5497043 | fk4_29_2 | RUNNING | 1-05:42:08 | 64 |
| A8 | 5478223_2 | find4 | RUNNING | 1-23:47:43 | 128 |
| A8 | 5478223_1 | find4 | RUNNING | 2-00:00:46 | 128 |
| A8 | 5478223_0 | find4 | RUNNING | 2-00:16:44 | 128 |
| A8 | 5478222_1 | finretry | RUNNING | 2-00:34:31 | 128 |
| A8 | 5478222_0 | finretry | RUNNING | 2-00:37:48 | 128 |
| A8 | 5497040 | fk4_29_0 | RUNNING | 1-06:38:51 | 64 |
| A8 | 5497042 | fk4_29_1 | RUNNING | 1-06:38:51 | 64 |
| A8 | 5478224 | finsu2 | RUNNING | 2-00:43:32 | 40 |
| A8 | 5410760 | quiverB | RUNNING | 2-01:29:23 | 12 |
| A6 | 7557993 | ls25 | PENDING | 0:00 | 64 |
| A6 | 7557992 | ls24 | PENDING | 0:00 | 64 |
| A6 | 7508416 | codex_su2_pilot | PENDING | 0:00 | 4 |
| A6 | 7570961 | ls23 | RUNNING | 8-22:02:07 | 64 |
