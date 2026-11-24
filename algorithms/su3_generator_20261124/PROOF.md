# Completeness and exact finite-rank test

## 1. Complete native inputs

The weight is `(3,0,0,4,4)` in flavour-then-dotted order. Nonnegative flavour
charges force every letter to contain either no flavour coordinate or only the
first flavour coordinate. Length seven therefore contains exactly three scalar
jets and four fermionic jets. Every derivative distribution of total dotted
degree `(4,4)` is included. Adjacent lengths contain the same three scalar jets
and respectively three or five fermionic jets.

`trace_algebra.py` enumerates every multiset of admissible letters, every
oriented cyclic trace, and every graded product of traces. Tracelessness removes
length-one traces. Graded cyclicity and exchange of odd traces remove exactly
their formal algebraic relations. Charge conjugation is imposed by its direct
coefficient equations. Dotted singlets are the complete integer kernel of the
raising action at weight zero. These finite-dimensional representations are
semisimple, so this kernel is the full dotted-singlet subspace.

The independent character calculation uses Schur-Weyl characters on invariant
colour tensors, subtracts identity-matrix slots by inclusion-exclusion, and
keeps the odd charge-conjugation character. Bosonic letter multiplicities use
symmetric powers; fermionic letter multiplicities use exterior powers. The
dotted-singlet multiplicity is the difference of the two adjacent weight
multiplicities:

| Length | SU(3) weight `(4,4)` | SU(3) weight `(5,3)` | Difference |
|---:|---:|---:|---:|
| 6 | 9,098 | 7,929 | 1,169 |
| 7 | 12,717 | 10,821 | 1,896 |
| 8 | 9,907 | 8,319 | 1,588 |

In the stable trace algebra the corresponding differences are

```text
16270 - 14179 = 2091
33400 - 28409 = 4991
41509 - 34900 = 6609
```

The exact polynomial bases attain these dimensions. Formal singlets have
distinct final nonredundant parity-word coordinates, which gives their
triangular independence witness. The full construction also verifies every
raising equation and every charge-conjugation equation directly.

## 2. Exact finite-colour identities

At each fixed derivative/type multiplicity, the evaluator assigns traceless
integer `3 × 3` matrices to every allowed jet. Each matrix is the sum of two
signed partial permutation or diagonal matrices, so its maximum absolute row
sum is at most two. Fermionic letters take values in an exterior algebra with
one more generator than the fermion degree.

For a trace word of length seven or eight, an exterior coefficient is bounded
by

```text
(fermion degree)! × 2^(word length) × 3^(number of traces).
```

The coefficient sum of each formal singlet gives the corresponding combination
bound. Every bound is checked before signed-64-bit native evaluation. Exact
rational solves then express every formal singlet in the chosen finite basis.
Each chosen evaluation minor is a nonzero reduction of an integer minor.

The independent character dimensions and these minors agree, establishing a
complete injective coordinate frame. Hence an exact zero in every coordinate
is an exact polynomial identity in the native SU(3) space. This implication
uses the complete character-certified space.

There are 26, 18 and 10 frames at lengths six, seven and eight. Their complete
relation dimensions are

```text
2091 - 1169 = 922
4991 - 1896 = 3095
6609 - 1588 = 5021.
```

For each formal row outside the chosen basis, the relation is that row minus
its printed rational linear combination of finite-basis rows. Every such
relation is present in `frames.json.gz`; none is selected by an ansatz.

## 3. Exact differential and ordinary cohomology

The differential is derived from the complete relative superfield action.
`maps.py` first differentiates every source polynomial in the stable trace
algebra. Its coordinates in the complete target singlet basis are determined
by triangular coefficient elimination. The entire resulting oriented trace
polynomial is compared with the direct derivative, including every term.

The finite maps are then computed through the complete target coordinate
identities. For every formal source row, the construction verifies equality
between differentiation followed by finite coordinates and finite source
coordinates followed by the finite differential. In particular every one of
the 922 and 3,095 source relations has zero finite derivative. Exact graded
nilpotence is verified before the colour identities are imposed.

Direct rational ranks of the full finite matrices give

```text
rank(1169 × 1896 incoming matrix) = 825
rank(1896 × 1588 outgoing matrix) = 1037
1896 - 1037 = 859                  complete closed space
859 - 825 = 34                    complete classical cohomology.
```

The normalized representative is the full rational vector printed in
`generator.json.gz`. It has 848 nonzero finite-basis coordinates. Substitution
into the complete outgoing matrix gives exactly zero. Its independently
expanded trace polynomial contains 12,224 terms. The full stable derivative
has 39,829 terms and 3,177 stable singlet coordinates. The explicit identity in
`closure_relation.json` expresses it as a combination of 2,047 target relations
from the complete 5,021-dimensional relation space.

## 4. Full dual and every proper product

The rational functional is defined on all 1,896 native finite cochains through
the complete basis. Its 349 nonzero coordinates are printed. Exact arithmetic
checks its value zero on all 1,169 incoming columns and on all 33 proper-product
representatives. Its value on the normalized operator is one.

The complete B04 proper-factor certificate covers all 56 pairs, 112 factors and
114 adjacent matrices. Tracelessness restricts factor lengths to `2 + 5` and
`3 + 4`. All nonnegative flavour splittings and every pair of equal dotted
spins are included. Its full target Cartan-injectivity certificate removes
the factors with incompatible colour parity. The character/rank bounds for
all 112 factors are saturated by their graviton products. Therefore every
proper product is represented by the 33-dimensional span tested by the dual.
The exact pinned archives and their replay scope are recorded in the README.

The current package reconstructs these 33 products as actual non-Abelian trace
polynomials. It checks their full stable derivative, dotted raising equation
and odd charge conjugation exactly. Thus the functional is zero on every
boundary and every proper product, and equal to one on the printed class.

## 5. Complete stable kernel

Every graviton trace recipe of length at most seven is included. There are 447
recipes at dotted weight `(4,4)` and 361 at `(5,3)`. The exact singlet equations
give 86 closed combinations. These are actual stable trace polynomials whose
classical derivatives vanish before finite-colour relations.

The complete stable outgoing map has a nonzero rank-3,348 minor. The incoming
map together with the 86 closed combinations has a nonzero rank-1,643 minor.
All these vectors are exactly closed. Therefore

```text
3348 + 1643 = 4991
```

proves that these vectors span the entire stable kernel over the rationals.
No assumption about the general stable-cohomology theorem is needed.

Their finite images have a rank-858 minor. The rational functional vanishes on
all those images and is nonzero on a finite closed vector. Since the entire
finite closed space has dimension 859, those images have dimension at most
858. The matching lower bound proves that the full finite image of the stable
kernel has dimension exactly 858.

## 6. Complete relation quotient and fortuity

The space of closed source relations has dimension

```text
1643 - 858 = 785.
```

Consequently the image of the entire 3,095-dimensional source relation space
has dimension

```text
3095 - 785 = 2310.
```

The full inverse image of the finite closed space has dimension

```text
3095 + 859 = 3954.
```

Its differential lands in the complete target relation space. Removing the
entire stable kernel gives

```text
3954 - 1643 = 2311
2311 - 2310 = 1.
```

This is the complete connecting relation quotient. The printed operator is
nonzero in it. Indeed, if its stable derivative were the derivative of any
source relation, subtracting that relation would give a stable closed
polynomial with the same finite image. The full dual vanishes on every such
image, whereas its value on this operator is one. This exact contradiction
establishes the nonzero relation class and therefore fortuity.

## 7. Independent check and replay limits

`independent_verify.py` implements the oriented graded trace algebra and the
one-flavour letter differential separately. It reconstructs the entire
printed polynomial and derivative, verifies both spin actions and charge
conjugation, and evaluates every term of the residual at the complete 1,588
target coordinates through a separate NumPy exterior-matrix implementation.
Large polynomial coefficients are combined with unbounded Python integers.
Every residual is zero. All dual pairings are checked independently as exact
fractions.

`reproduce.py` starts with an empty work directory and reconstructs every new
native polynomial and finite-frame entry. The nine core mathematical files
must reproduce byte for byte. The inherited B04 proper-factor minor entries
are distinguished from this fresh reconstruction throughout.

