# Uniqueness of the full M4 class

All spaces below are gauge and rotation singlets at (n1,n2) = (7,7). The differential Qn lowers Nb by n and raises Nc by one. E1 denotes Q0 cohomology; d1 is induced by Q1.

At Nc = 5, exact source dimensions and rank lower bounds are:

| Nb | Dimension | Incoming Q0 rank | Outgoing Q0 rank | E1 dimension |
|---|---:|---:|---:|---:|
| 0 | 136 | 39 | 97 | 0 |
| 1 | 474 | 221 | 246 | 7 |
| 2 | 358 | 239 | 119 | 0 |
| 3 | 68 | 57 | 10 | 1 |
| 4 | 3 | 2 | 1 | 0 |

Nb >= 5 has zero source dimension. The Nb = 3 class has nonzero Cartan restriction, so the displayed upper bound one is attained. The Nb = 1 ranks are exact: at Nc = 4 the dimension 327 and ranks 106 + 221 saturate; at Nc = 6 the dimension 351 and outgoing rank 105 force the preceding rank 246 to be exact.

The block differential Q0 + Q1 at consecutive Nb values gives the exact induced ranks

$$
\operatorname{rank}d_1^{(2,4)}=464-239-221=4,
\qquad
\operatorname{rank}d_1^{(1,5)}=346-246-97=3.
$$

These finite-field lower bounds saturate the seven-dimensional E1 sector, since d1 squared is zero. Therefore

$$
\dim E_2^{(1,5)}=7-4-3=0.
$$

The only possible incoming d1 source to the Nb = 3 class, at (Nb,Nc) = (4,4), has dimension 14 and ranks 12 + 2, hence zero E1. Higher incoming sources have Nb >= 5 and Nc = 4 and vanish. The complete rational cocycle in `M4_full_exact.pkl` proves that the Nb = 3 class has no outgoing differential. Therefore it is the sole surviving class at Nc = 5, and any two full lifts of the same leading class differ by a full perturbative boundary.

All rank outputs RQ0_*.out and RD_*.out report bad = 0. A nonzero finite-field minor of an actual rational invariant-image matrix proves the rational rank lower bound. Exact character dimensions and the chain relation supply the matching upper bounds; this argument does not infer rational cohomology from an unchecked modular rank equality.
