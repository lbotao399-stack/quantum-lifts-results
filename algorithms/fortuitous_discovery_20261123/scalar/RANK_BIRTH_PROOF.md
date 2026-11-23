# Strict rank birth: (7,7,3) and (9,5,3) in SO(15)

These are the selected two sectors at j=34, Y=16. The (7,5,5) control is
excluded from novelty: its defect was already nonzero in an earlier
SO(13) product calculation.

| Shape | SO(13), I17 / Q0(I16) | SO(15), I17 / Q0(I16) |
|---|---|---|
| (7,7,3) | 0, exact rational relation preimage | nonzero, exact dual value 1 |
| (9,5,3) | 0, exact rational relation preimage | nonzero, exact dual value 1 |

Vanishing on the left means vanishing of the **fortuitous quotient**.
The corresponding finite cycle may still represent a monotone class.
No vanishing of its entire smaller-rank cohomology class is asserted.

## Complete lower-rank frames

The SO(13) character gives source/terminal dimensions 81/45 for (7,7,3),
and 179/90 for (9,5,3). The self-contained `verify_rank_birth.py` rebuilds
these bounds.

`SO13_773_EXACT_BOUNDARY.json` contains 44 actual Q0-image rows and one
actual scalar-product relation. The 45 rows are independent on the
specified 45 probes and exhaust the character bound. They give an
injective coordinate frame for all terminal relations. The rational
preimage uses only the first 44 rows: the product coordinate is exactly
zero. It matches all 50 entries, including the unused five.

`SO13_953_EXACT_BOUNDARY.json` analogously contains 88 actual Q0-image
rows and two actual product relations. The 90 independent rows exhaust
the character bound. Both product coordinates of the rational preimage
are exactly zero, and all 100 entries match.

The deterministic C++ `mode=1` recipes differentiate sixteen-letter
finite identities of Pfaffian degrees 14,15,16, all greater than thirteen.
The `mode=4` recipes multiply a degree-fifteen Pfaffian by a scalar
bilinear. All rows are full finite-colour relations. The probe banks use
three epsilon triples and four or two AB wedges, respectively; the second
bank has four remaining fixed A slots. They select precisely the highest
weights (7,7,3) and (9,5,3).

## Identification with the explicit SO(15) defects

Odd-degree invariants of three auxiliary SO(3) vectors are the determinant
times polynomials in the six Gram pairings, with determinant-square equal
to Gram determinant. The symmetric algebra of Sym^2(C^3) has one copy of
each all-even three-row diagram. Thus the all-singleton Pfaffian family
has one highest-weight line for each all-odd three-row shape.

The raw lower-rank seed and the explicit operators D^3 ZZ^2, D^3 AA^2 ZZ
lie on the same line. `RAW_TO_EXPLICIT_PFAFFIAN.json` checks their ratios
on two fully expanded multilinear probes:

```
raw_773 = (-4/3969) explicit_773,
raw_953 = (256/1091475) explicit_953.
```

Dividing the stored raw preimage by the stated nonzero ratio gives the
preimage of the exact defect used in the SO(15) representative. This is
a comparison of specified coefficient formulas, not a rescaling of Q0.

## A priori integer bound for CRT

In the modes used here the C++ program uses only integer additions and
multiplications. The declared modular exponentiation function is never
called. The circular marked-tile formula absorbs the trace-log
denominators; the probe sums have integer epsilon/wedge coefficients.
Every printed residue is therefore the reduction of a definite integer.

The following bound holds before reconstruction:

1. Seed parameters lie in [-4,4]. Singleton entries have absolute value
   at most4. A marked path has at most three original slots. Replacing
   its matter slot by a commutator may produce four scalar sites but no
   extra parameter. Its matrix/slot factors give at most4^4, and the
   original/reversed pair gives at most2*4^4=512. Closed-bilinear entries,
   including their cyclic/reversal cuts, also obey512.
2. There are seventeen physical sites. Bound every tile entry by512,
   including all much smaller singletons. Expanding the auxiliary matrix
   traces gives at most2^17 summands. One product of connected coefficients
   is bounded by2^17*512^17=2^170 before marked placements.
3. There are at most four labelled marked objects. Each is absent or has
   a starting site among seventeen. Its length is fixed by its recipe;
   all remaining sites have the unique unmarked singleton. Thus at most
   18^4<2^17 placements occur. The circular algorithm identifies a tiling
   uniquely by the tile containing site zero, so it adds no cyclic
   multiplicity. Counting a closed spectator by a starting site only
   increases this upper bound.
4. The larger bank has6^3*2^4=3456<2^12 flavour assignments. Allow an
   additional factor17!<2^49 for every possible labelled-slot or repeated
   factor multiplicity, even though the implementation already accounts
   for these within its circular and probe sums.

Hence

```
abs(entry) < 2^(170+17+12+49) = 2^248 < 2^256.
```

The product of the sixteen stored CRT primes exceeds2^257. Balanced CRT
therefore determines every integer uniquely in the allowed interval.
The seventeenth prime is an extra replay check. Exactness does not depend
on the small recovered entries or on that extra prime. The maximum
recovered absolute entries are1267913309626368 and828941948003328.

## Rank support and products

For fixed universal defects the ideals obey I_SO(15) subset I_SO(13).
The exact SO(13) preimages also lie in every odd orthogonal relation
ideal at N<=13. Thus the fortuitous components vanish at every such rank.
They are nonzero at N=15 by the complete small rational matrices.
For odd N>=17 there are no length-seventeen relations: the first possible
orthogonal metric relation has length N+1. The displayed nonzero scalar
coefficients imply that the universal defect is nonzero, so the same
primitive is not even Q0-closed there. The fortuitous support of these
fixed representatives among odd orthogonal groups is exactly SO(15).
This statement is not extended to even orthogonal groups.

At SO(15), a nontrivial product of length sixteen has factors of length
at most fourteen. Their Q0 images have length at most fifteen, within
the injective specialization range. Each finite closed factor lifts to
a stable closed factor; their product is monotone and is killed by both
exact duals. Proper ordinary BMN raising ancestors are also monotone:
at length sixteen and level thirty-two all letters are scalars, and at
lower levels no such cochain exists.

Historical novelty remains an inventory statement; the algebraic claims
above are exact. No quantum lifting or all-loop BPS survival is asserted.
