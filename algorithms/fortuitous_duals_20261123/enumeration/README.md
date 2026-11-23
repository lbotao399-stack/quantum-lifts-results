# Dual-first fortuitous recognition on real Sp(3) data

## Main result / 本轮新增

The known Sp(3)=USp(6) omega class at `j=28`, native weight tuple `(3,3,5,0,0)`,
length 9 now has an explicit **40-term integer functional** on formal
length-10 trace coordinates satisfying, over the rationals,

```
ell Q0(I) = 0 on the COMPLETE 181-dimensional relation space;
ell Q0(omega lift) = 1.
```

Here `I` is the kernel of finite-colour trace specialization in the complete
highest-weight source block, and `ell` is a coefficient functional, not a
physical Wick bra. The native tuple's first two entries are the dotted/Taylor
weights and its last three entries are the flavour weights.
The source is the original 4264-term rational omega from the repository;
its representative and normalization are unchanged. All forty native ordered
trace words and their integer coefficients are in `sp3_intrinsic_dual.json`.
The largest absolute coefficient is 282.

This is an exact fortuitous-obstruction certificate for an existing class.
It is neither a new BPS class nor a new quantum differential. The new
contribution is a complete-relation, small-support detector with an explicit
replay, rather than a nonzero score on one chosen lift without relation control.

The exact pairing, grouped by supported trace-length partitions, is

```
52793/144 - 2989/9 + 157/3 + 1139/36 - 5639/48
= (52793 - 47824 + 7536 + 4556 - 16917)/144
= 144/144
= 1.
```

The respective partitions are `(2,2,2,4)`, `(2,2,3,3)`, `(2,2,6)`,
`(2,3,5)`, and `(2,8)`. These individual contributions are not separate
cohomology invariants. Their sum is protected against every source trace
identity by the full annihilation equation.

## 1. Why this detects the fortuitous quotient

Use the repository's trace-specialization map `pi:C_st -> C_N` and
`I=ker(pi)`. Assume source and predecessor-degree surjectivity and the
chain-map identities. All degrees below are the fixed source grading and
the relevant adjacent grading; no physical normalization is changed.

For a finite closed cochain, choose a stable lift. Its derivative lies in
the next-degree trace identities. A different lift changes that derivative
by `Q0(I)`. A finite boundary has a stable predecessor lift; its derivative
therefore vanishes modulo `Q0(I)` as well. A universally closed class also
has zero derivative. Thus a functional annihilating all `Q0(I)` and pairing
nontrivially with the chosen derivative proves a nonzero fortuitous component.

The exact image of this connecting map is

```
(im Q_st intersect I_next) / Q_st(I).
```

The specialization of `im Q_st` onto `im Q_N` is surjective. Consequently,

```
dim(im Q_st intersect I_next) = rank Q_st - rank Q_N,

dim fortuitous quotient = rank Q_st - rank(Q_st restricted to I) - rank Q_N.
```

No basis of stable cohomology, finite harmonic cohomology or graviton products
is required for this count. **All ranks in this subtraction must be exact**;
subtracting unrelated modular rank lower bounds would not prove the result.

For the actual Sp(3) block,

```
formal source dimension          2870
finite source dimension          2689
complete source-relation dimension 181
rank Q_st                       1198
rank Q_st(I)                     116
rank Q_N                        1081

1198 - 116 - 1081 = 1.
```

The companion calculation certifies all 181 source identities over Q on a
complete 2689-coordinate finite frame. Its three-prime calculation uses
identical integer matrices and a rigorous coefficient bound:

```
CRT modulus = 1000073001431003663
2 * integer evaluation bound = 53644213685846016
1000073001431003663 > 53644213685846016.
```

The finite-frame minor is nonzero modulo 1000003, so the complete frame is
injective over Q. Exact zero evaluations therefore imply genuine identities,
not only zeros at sampled points. The 181 independent identities exhaust
the source kernel by the exact source/finite dimension difference.

`sp3_intrinsic_dual.py` consumes that certified rational basis, checks its
decoded SHA256, and then performs exact integer/Fraction annihilation and
cycle pairings. The finite field is used only to propose the forty entries;
the final proof is over Q.

## 2. Equivalent direct dual equation, without forming I

Let `b` be a linear functional on the finite source cochains and `ell` a
functional on the stable target coordinates. With the existing matrix
conventions, solve

```
ell Q_st = b pi.
```

For a universally closed source, the right-hand side is zero. For a finite
boundary, lift its predecessor and use `Q_st^2=0`; the right-hand side again
vanishes. Thus such a `b` is a cohomology functional annihilating the stable
image. Functionals in `row(Q_N)` vanish on every finite cycle, so they must
still be quotiented out; omitting this quotient would overcount.

The solution space is the row intersection

```
row(Q_st) intersect row(pi),
```

transported to the finite source via the injective pullback of functionals.
Its dimension is

```
rank Q_st + rank pi - rank([Q_st; pi]),
```

where the last matrix stacks the two maps with the same stable source
columns. The fortuitous dual dimension is therefore

```
rank Q_st + rank pi - rank([Q_st; pi]) - rank Q_N.
```

This is the kernel-free row-space version of the connecting-map formula.
It is ordinary finite-dimensional linear algebra; we do not claim a new
homological theorem. Its practical role is to avoid constructing three large
kernel bases merely to enumerate the much smaller fortuitous quotient.

On the actual Sp(3) data, the unchanged archived finite functional has an
explicit factorization with 402 nonzero `ell` coefficients. It is verified
on all 2870 formal source columns, all 2624 formal predecessor rows, and
all 1281 outgoing cycle coordinates. The factorization minor has determinant
452105 modulo 1000003. The exact-rank dimensions are

```
row-intersection dimension       = 1082;
outgoing dual-coboundary dimension = 1081;
fortuitous dual quotient         = 1;
stacked rank                    = 1198 + 2689 - 1082 = 2805.
```

The 402-term factorization is explicitly modular. Its stacked rank is derived
from the complete companion complex, not from a separate dense elimination
of the stacked matrix. The intrinsic 40-term certificate above is stronger:
its relation annihilation and its omega pairing are exact over Q.

## 3. Which trace coordinates can detect the obstruction?

All ranks in this table were recomputed over Q from the certified rational
relation basis, using exact integer matrices after clearing denominators.
The last column adds the original rational omega derivative.

| Trace lengths | Coordinates | Rank Q0(I) | Rank after adding omega derivative |
| --- | ---: | ---: | ---: |
| (2,2,2,4) | 82 | 21 | 21 |
| (2,2,3,3) | 68 | 13 | 13 |
| (2,2,6) | 454 | 101 | 101 |
| (2,3,5) | 392 | 43 | 43 |
| (2,4,4) | 182 | 21 | 21 |
| (2,8) | 1444 | 116 | 117 |
| (3,3,4) | 68 | 13 | 13 |
| (3,7) | 452 | 43 | 43 |
| (4,6) | 458 | 102 | 102 |
| (5,5) | 272 | 54 | 54 |
| (10) | 2148 | 116 | 117 |

A single-trace-only detector exists, as does a detector restricted to `(2,8)`.
The other individual partitions cannot separate this particular derivative
from the projected relation image. This is a concrete support-selection rule
for subsequent detector searches. It is not a statement that the other
partitions vanish, nor a statement of minimum possible support.

## 4. Streaming dual updates: actual SU(3) algorithm control

If columns of `Z` span the annihilator of already admitted rows and `B`
contains new rows, update by

```
new dual columns = Z * kernel(B Z).
```

This equality is exact: a vector `Z c` annihilates the old rows automatically
and annihilates the new rows precisely when `(B Z)c=0`. The residual matrix
has only the current dual dimension as its width. The same update can stream
complete relation images or known universal-cycle images. To certify a final
dimension, every remaining row must be checked or excluded by a completeness
proof; a temporarily zero score is not a stopping certificate.

The actual 1454 by 434 SU(3) projected Q0 matrix gives the same sequence at
both primes 1000003 and 1000033:

| Source profile | Rows | Added rank | Residual width before update | Remaining dual dimension |
| --- | ---: | ---: | ---: | ---: |
| 38 | 791 | 424 | 434 | 10 |
| 39 | 321 | 6 | 10 | 4 |
| 41 | 160 | 0 | 4 | 4 |
| 53 | 182 | 0 | 4 | 4 |

The Gaussian-elimination matrices contain 347872 entries in total, versus
631036 entries for one full matrix. The 1024-dimensional source kernel is
never constructed. All 5816 source/dual pairings and all four actual RHS
pairings are checked at each prime. Four alternating paired timing trials
at each prime all favored the incremental arithmetic; the nanosecond values
are saved verbatim. These timings exclude native differential generation and
do not establish a speed gain for the full physical solver.

This SU(3) test concerns only the projected target31 cokernel. It is not a
fortuitous census or a solution of the coupled second-correction equation.

## Reproduction / 复现

Dependencies are listed in `requirements.txt`. First run the companion
`positive_sector` regeneration and complete rational-relation certificate.
Then, from this directory:

```bash
python3 sp3_dual_factorization.py --input ../positive_sector
python3 sp3_intrinsic_dual.py --input ../positive_sector
```

The first script checks the full modular direct-dual equation. The second
regenerates the forty entries, verifies the complete rational relation basis,
checks the unchanged rational omega, and recomputes every partition rank
over Q. It takes no physical Gram matrix as input.

Optional SU(3) algorithm control, using the saved arrays from the preceding
frontier audit:

```bash
python3 streaming_profile_duals.py --input /absolute/path/to/su3_second
```

The code preserves native word tuples and source order. It introduces no new
trace normalization. Source omega and dual-singlet algorithms are pinned to
`9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`; their shared runtime is pinned to
`082be7070d207e2cdecfd0a91b03d5e014fd303d` in
`lbotao399-stack/quantum-lifts-results`.

### Scientific limits

A detector score becomes a cohomology statement only after finite-colour
closure is established. Candidates must lie in this highest-weight block;
an arbitrary unprojected weight-space expression must first be projected or
separately certified highest weight. A nonzero score on a nonclosed expression
does not certify a BPS state. The forty-term functional detects the known
one-dimensional fortuitous component; it is not automatically a local CFT
observable or a metric-orthogonal projector. The construction offers a
quotient certificate and a way to prioritize candidate enumeration; physical
matrix elements require additional data.
