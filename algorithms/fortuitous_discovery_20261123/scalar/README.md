# Two strict SO(15) fortuitous rank births

The selected classical classes are

```
SO(15), (7,7,3;0,0)_16, j=34;
SO(15), (9,5,3;0,0)_16, j=34.
```

Each has fifteen scalar letters, one matter fermion, and one exact
fortuitous quotient direction. Their scalar defects are Q0-boundaries
of finite identities in SO(13), but are not such boundaries in SO(15).
This excludes the lower-rank persistence that invalidated the SO(13)
and SO(15) `(7,5,5)` preliminary candidates.

| Highest weight | SO(15) source relations | SO(15) terminal relations | Exact Q0 rank | Fortuitous quotient |
|---|---:|---:|---:|---:|
| (7,7,3) | 2 | 3 | 2 | 1 |
| (9,5,3) | 3 | 4 | 3 | 1 |

## Full replay

From this folder:

```sh
python verify_rank_birth.py --rebuild-lower
```

This needs Python 3's standard library and a C++17 compiler. It rebuilds
the small SO(15) matrices, actual equivariant primitive, complete source
and terminal character dimensions, the SO(13) rational boundary checks,
and every lower-rank integer entry modulo the full CRT certificate.
No python-flint, NumPy, external repository, or network connection is
required. Without `--rebuild-lower`, it verifies the stored lower-rank
integer matrices and exact rational identities and rebuilds all other
parts.

`RANK_BIRTH_PROOF.md` supplies the complete-frame argument, the explicit
integer reconstruction bound, the comparison of raw and specified
Pfaffians, and the exact rank-support statement. `AUXILIARY_CONVENTIONS.md`
defines the native letters, Pfaffian series, auxiliary coefficient
operators and the flavour-equivariant odd-trace primitive. That file
also documents the SO(13) controls; those controls are not counted here.

## Exact small matrices

For `(7,7,3)`, use these ordinary ordered-trace coefficients:

```
AA.AA.AA.ABC.BB.BB.BB.CC
AA.AA.ABC.ABC.ABC.BB.BB
AACBBC.AB.AB.AB.ABABC
```

The columns are the two complete predecessor images and the specified
scalar Pfaffian:

```
[    -57960,    -7560,    272160 ]
[ -12983040, -1693440,  20321280 ]
[  32336640,  2419200, -52254720 ]
```

```
det = 4236912357801984000,
dual = (1/181440, -1/40642560, 0),
dual(Q0 sources) = (0,0),
dual(Pfaffian) = 272160/181440 - 20321280/40642560
               = 3/2 - 1/2 = 1.
```

For `(9,5,3)`, use

```
AA.AA.AA.AA.ABC.BB.BB.CC
AA.AA.AA.ABC.ABC.ABC.BB
AA.AAAB.AB.AB.ABBCC.AC
AA.AA.AABAAC.ABC.BB.BC
```

The columns are the three complete predecessor images and the specified
scalar Pfaffian:

```
[     86940,    49140,    11340,   -374220 ]
[  17388000,  9828000,  2268000, -24948000 ]
[   1179360,   997920,  -272160,         0 ]
[ -19051200, -4536000, -1088640,  35925120 ]
```

```
det = 13523262578496405504000000,
dual = (-1/249480, 1/49896000, 0, 0),
dual(Q0 sources) = (0,0,0),
dual(Pfaffian) = 374220/249480 - 24948000/49896000
               = 3/2 - 1/2 = 1.
```

The character dimensions and nonsingular displayed matrices prove that
these coefficient frames separate the entire terminal relation spaces.
Their duals annihilate the complete relation image, not a sampled part.

## Complete representatives

The scalar Pfaffian operators in the auxiliary conventions are

```
(7,7,3): D^3 ZZ^2;
(9,5,3): D^3 AA^2 ZZ.
```

All ordinary differentiation factors are retained. The full cochain is
the flavour-equivariant cyclic/reversal averaged odd-trace primitive
applied coefficientwise to the corresponding scalar Pfaffian. Its
universal Q0 derivative is that exact Pfaffian; its derivative vanishes
in fifteen colours. The full coefficient oracle permits direct queries:

```sh
python minimal_scalar_certificates.py --colours 15 --shape 7 7 3 \
  --query-word AA.AA.AA.Aa.BB.BB.BB.CC --native-sorted-theta
# 45360

python minimal_scalar_certificates.py --colours 15 --shape 9 5 3 \
  --query-word AA.AA.AA.AA.Aa.BB.BB.CC --native-sorted-theta
# -62370
```

Here native `a` is the sorted Taylor matter letter with indices `(1,2)`.
The optional query conversion also handles the minus sign between cyclic
`psi_B` and the native sorted `(0,2)` Taylor letter. Q0 is unchanged.

## Exact smaller-rank boundaries

The lower-rank files contain actual full finite identities and every
rational coefficient of the preimages:

| SO(13) shape | Complete terminal frame | Actual Q0 rows in the preimage | Product coordinates | Checked rational coordinates |
|---|---:|---:|---|---:|
| (7,7,3) | 45 | 44 | 0 | 50 |
| (9,5,3) | 90 | 88 | 0,0 | 100 |

The raw lower-rank Pfaffians equal `(-4/3969)` and `(256/1091475)` times
the specified defects respectively, as checked on explicitly expanded
multilinear probes and the unique auxiliary highest-weight line. The
preimages therefore apply to the exact representatives used here.

This proves zero smaller-rank fortuitous component, not zero total
cohomology. By nesting of the odd orthogonal identity ideals, the same
fortuitous components vanish at every odd N<=13. At N=15 they are
nonzero. At odd N>=17, no degree-seventeen relation exists and the
specified primitive is not closed. Their fortuitous rank support among
odd orthogonal groups is exactly SO(15).

## Products, descendants and scope

At length sixteen every factor in a nontrivial product has length at
most fourteen. Its derivative has length at most fifteen, where SO(15)
specialization is injective. A finite closed factor therefore has a
stable closed lift, and every such product is monotone. Both duals kill
the whole stable cohomology image. Proper ordinary BMN raising ancestors
are scalar at level thirty-two, or absent below it, so they are monotone
as well.

These are two independent classical Q0 fortuitous classes in distinct
highest weights. The package asserts no quantum survival, first or
higher differential, or all-loop BPS status. No external publication or
message to another researcher has been made.
