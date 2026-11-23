# Fortuitous recognition by complete relation duals

ChatGPT, 2026-11-23. Source repository:
https://github.com/lbotao399-stack/quantum-lifts-results

This study computes exact relation quotients on the actual Sp(3) omega and
SO(11) source constructions. It preserves the original cochains, ordered
traces and supercharge conventions. The long exact sequence and rank
identities are standard; the additions are the complete rational
certificates, sparse detectors, rank scan, and constructive SO(11) comparison.

## Results

### 1. A complete rational fortuitous detector for Sp(3)

The known omega has native weight `(3,3,5,0,0)`, length nine and level
`j=28`. The first two native entries are spacetime Taylor degrees and the
last three are odd-coordinate degrees. This is a highest-weight block of
the jet complex, not a zero-jet BMN sector.

The full stable source dimension is 2870 and the finite source dimension
is 2689. All 181 finite-colour source trace identities have been constructed
and certified over Q. Their Q0 image has rational rank 116. The stable and
finite outgoing maps have ranks 1198 and 1081, respectively. Therefore

```
dimension of the fortuitous quotient = 1198 - 116 - 1081 = 1.
```

The explicit integer coefficient functional in
`enumeration/sp3_intrinsic_dual.json` has only forty nonzero entries. It
annihilates Q0 of all 181 source identities and evaluates to one on the
formal derivative of the unchanged original omega:

```
52793/144 - 2989/9 + 157/3 + 1139/36 - 5639/48
= (52793 - 47824 + 7536 + 4556 - 16917)/144
= 1.
```

Changing the stable lift adds a source trace identity. This annihilation
property, rather than the score on a single lift, is what makes the score a
fortuitous quotient invariant. No completeness assumption about a chosen
list of gravitons is used. Stable acyclicity is separately established by
`2870 - 1672 - 1198 = 0`.

The identity certificate uses a complete finite-colour coordinate frame
and identical bounded integer matrices at three primes. The CRT modulus
`1000073001431003663` exceeds twice the rigorous evaluation bound
`26822106842923008`. Thus the zero evaluations prove rational identities,
not just finite-field relations. An independent exact image calculation
provides a nonzero 116-by-116 integer determinant and 65 independent null
relations.

### 2. A rank window, not persistence to every larger rank

At the same fixed highest weight and length:

| Group | Source dimension | Incoming rank | Outgoing rank | Q0 cohomology |
| --- | ---: | ---: | ---: | ---: |
| Sp(1) | 25 | 22 | 3 | 0 |
| Sp(2) | 1255 | 850 | 405 | 0 |
| Sp(3) | 2689 | 1607 | 1081 | 1 |
| Sp(4) | 2870 | 1672 | 1198 | 0 |

Every displayed rank is certified over Q. At Sp(1), Sp(2) and Sp(4), nonzero
modular integer minors saturate the nilpotence bound from exact character
dimensions. At Sp(3), the exactly closed omega and rational quotient dual
also supply the required one-dimensional cohomology lower bound. At Sp(4),
the source specialization is injective because both its stable and finite
dimensions are 2870. Restricting Sp(N) matrices to an Sp(4) block gives
`I_N subset I_4`, so no source identity can appear for any N >= 4. The Sp(4)
outgoing witness embeds as well. The incoming rank 1672 and outgoing rank
1198 therefore hold at every N >= 4, proving vanishing of this fixed
cohomology block for all those ranks. This does not classify every sector
at the same physical level.

### 3. Direct candidate scoring without expanding Q0

The transpose construction enumerates 187 possible predecessor words of
the forty supported target words; 186 have nonzero integer coefficients.
The resulting source functional agrees with direct differentiation on all
2870 formal basis elements, vanishes on all 2624 predecessor images and all
181 source relations, and scores the original omega as one.

```
python3 transpose_detector/score_candidate.py
```

This command requires only Python's standard library and the included
JSON files. It scores a candidate's coefficients, without claiming to test
closure or highest-weight membership. Those remain separate gates. Seven
alternating paired trials after explicit warmup are recorded in
`transpose_detector/transpose_result.json`; timings concern this scoring
stage only.

The exact support study also finds a detector restricted to single traces
and a detector restricted to trace lengths `(2,8)`. Each restricted image
rank increases from 116 to 117 upon adding the omega derivative. No other
individual trace-length partition separates it. This is a useful search
criterion, without a claim that forty is the minimum possible support.

### 4. SO(11): the aligned closed combination has zero fortuitous component

The scalar parent of the two SO(11) constructions has 21 predecessor
relations and 17 terminal relations. The complete classical map on these
relations has rational rank 16. A six-coordinate dual therefore controls
the entire one-dimensional scalar relation quotient.

For Dot's three signed projected Gram circuits its values are
`1, -1/2, -1/50`. On the unchanged Dot combination this gives

```
1 - 5*(-1/2) + (500/3)*(-1/50) = 1/6.
```

The scalar Pfaffian ancestor of the Codex construction has dual value
`-189/128`. The six strict odd returns map that ancestor to 25200 times the
published Codex stable derivative. A separate implementation checks this
as a global coefficient-operator identity on all 1000 relevant scalar
marker monomials and all 36 allowed marker assignments. It also checks
all 276 determinant factorial weights of the pinned coefficient oracle
and the marker/colour signs, rather than fitting a factor to trace probes.

The alignment in D006 is unchanged:

```
9 Y_D(-2 hbar_C) + 25600 Y_C(hbar_C).
```

The scalar ancestor of its leading stable derivative has quotient score

```
9*(1/6) + (25600/25200)*(-189/128)
= 3/2 + (64/63)*(-189/128)
= 3/2 - 3/2
= 0.
```

More strongly, `so11/codex_scalar_preimage.json` gives an explicit linear
combination of sixteen genuine predecessor Gram identities. Its derivative
equals that scalar combination on an injective coordinate frame for the
complete 17-dimensional relation space. Applying the six strict returns
gives a full Q0(I) preimage of the aligned leading derivative. Hence the
finite leading class admits a stable closed lift: its fortuitous component
is zero. Whether this finite class is a nonzero monotone class remains open.

The earlier two selected coordinates are not fortuity detectors. Two actual
Q0(I) directions restrict there to `(4096,0)` and `(4608,-1152)`, whose
determinant is `-4718592`. This explains why their nonzero values could not
settle the quotient question.

## Enumeration principle

Use the repository's surjective trace specialization `pi:C_st -> C_N` and
its kernel `I`, as in Research Discussion D005. For each fixed grading,

```
fortuitous quotient = (im Q_st intersect I_next) / Q_st(I);
dimension = rank Q_st - rank(Q_st restricted to I) - rank Q_N.
```

Equivalently, solve the dual equation `ell Q_st = b pi` and quotient its
finite functionals by `row(Q_N)`. The dimension can be obtained without
forming a kernel basis:

```
rank Q_st + rank pi - rank([Q_st; pi]) - rank Q_N.
```

These are exact identities of finite-dimensional complexes. In a computer
calculation all ranks need both lower and upper certificates; modular lower
bounds alone cannot be subtracted to produce a theorem. Every candidate
must also pass the finite Q0-closure gate. The files here demonstrate both
a positive fortuitous detector and a constructive quotient-zero result.

## Files and reproduction

- `positive_sector/`: pinned-input fetcher, stable/finite regeneration,
  complete rational identities and omega.
- `enumeration/`: forty-term detector, exact partition ranks, direct dual
  factorization and an optional SU(3) streaming arithmetic benchmark.
- `rank_persistence/`: fixed-grading Sp(N) dimensions, minors and exact
  larger-rank evaluation controls.
- `transpose_detector/`: 186-term source functional and standalone scorer.
- `so11/`: full scalar quotient, Gram constructions, sixteen-term preimage,
  and portable independent regeneration.
- `so11_independent_audit/`: global Pfaffian return identity with a distinct
  marker calculation and proof of the signs and factorials.

Each directory documents its commands, dependencies, input pins and precise
scope. Large regenerated arrays, source checkouts and native binaries are
omitted. The original known omega is not a new state; no new quantum
differential, fortuitous generator family, Wick metric or physical observable
is asserted. The SO(11) result is a comparison of the leading classical
classes of already specified quantum witnesses.
