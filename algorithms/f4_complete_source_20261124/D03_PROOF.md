# Complete sixth-power second-map decision

**Result.** The entire incoming second differential to the sixth F4 glueball
has rank zero over the rationals in the relative holomorphic jet complex.
The complete-rank input is the previously published theorem; the direct
closed-factor calculation and literal-tensor checks below are fresh.

## Complete source and matching rank bounds

The full Taylor Ward identity fixes the only possible incoming highest weight
to dotted/flavour degree `(4,4,2,2,2)`, hence Taylor diagram `(4,4,3,3)`.
The ten spin/flavour singlet types are retained in the character inventory.
An equivariant nonzero map to the pure sixth-power target would have to be
nonzero on a compatible highest-weight vector. The other nine types have
incompatible highest weights. This argument uses the simultaneous raising
equations and covariance; it does not apply an idempotent to a cochain.

At lengths twelve, thirteen and fourteen the selected cochain dimensions are

```text
10007, 3720, 704.
```

The independently published rational outgoing rank is 703. Its upper bound is
the nonzero Cartan polynomial with value `1728*5184=8957952`; its lower bound
comes from order-703 minors with residues 1005 and 734. The exact dimension of
the full classical kernel is therefore

```text
3720 - 703 = 3017.
```

The public fixed kernel minor has 3016 actual classical boundaries and the
explicit product discussed below. Its common rational definitions have regular
reductions at both primes. The published numerical certificates are:

| Prime | Old boundary minor, order 3002 | Fourteen-boundary extension | Product extension, order 15 | Complete kernel minor, order 3017 |
|---:|---:|---:|---:|---:|
| 1009 | 833 | 587 | 737 | 449 |
| 1013 | 615 | 953 | 346 | 60 |

```text
833*737 mod1009 = 449
615*346 mod1013 = 60
3016 + 1 = 3017 = 3720 - 703.
```

Either nonzero complete minor proves rational independence. Thus these cycles
span **every** classical source cycle. The proof does not assume that the
product is a nonzero classical cohomology class: it may itself be another
boundary. Nor does it claim a complete basis for all 10007 predecessors.

All pins and hashes are in `proof_dependencies.json`. The complete order-3017
minor is a historical attributed lemma, not a claimed fresh native replay.

## The actual product closes without a subspace projector

The polynomial is the exceptional nine-scalar invariant times a scalar
bilinear and a marked bilinear, with the coefficient constraints imposed
directly on their products. Exact rational elimination gives

```text
172 equations on 96 coefficients
150 raising equations + 22 classical equations
rank 95, kernel dimension 1
87 nonzero primitive integer coefficients.
```

The 22 classical equations are imposed for each fixed pair of spectator
polynomials. The 45 active spectator pairs therefore each multiply a separately
closed marked bilinear. The full coefficients and all zero residuals are in
`closed_product_rational.json`. No Taylor idempotent is used to construct this
highest-weight polynomial.

The complete invariant bilinear calculation is small enough to retain in full.
There are five first Taylor jets and twelve second Taylor jets, hence 60
two-letter invariant cochains of total Taylor degree three. The nonzero native
alternating three-colour tensor supplies 20 classical output monomials. The
exact matrix has

```text
classical rank 20
60 - 20 = 40 closed bilinears.
```

The first quantum operation needs all three scalar flavours. It can act only
on the three complementary matter-fermion/scalar columns. The complete native
adjoint-Casimir contraction equals 18 times the identity on all 52 colours,
checked in all 2704 matrix entries. On the three columns the classical
coefficients are `(1,-1,1)` and the quantum coefficients are `(36,-36,36)`.
More strongly, on **all 60 columns** the quantum row is exactly 18 times the
three-scalar row of twice the classical matrix. Its product with the entire
40-dimensional classical kernel is exactly zero. The matrices, ordered
letters and identity are in `closed_bilinears_rational.json`.

Every higher bare operation needs at least three input letters, so it vanishes
on a bilinear. Consequently all 40 closed bilinears are strictly closed under
the full bare perturbative differential. For the remaining first-jet factors,
an order-r positive bare operation requires 3r flavour slots on r+1 input
letters; first jets supply at most r+1 slots, and `r+1<3r` for every positive r.
These spectator invariants are also strictly closed.

The matched quantum product satisfies the binary chain identity on the full
complex. Iterated matched multiplication of the three closed factors therefore
gives a formally closed completion whose leading term is the stated ordinary
product. This uses the archived all-loop manuscript's `eq-chain` and
`b-00513`, pinned in `proof_dependencies.json`. Cross-factor quantum
contractions are part of that product; ordinary tree multiplication is not
being substituted for them. No broad assertion about projected products is
needed for this explicit product.

## All earlier corrections remain free

The complete first-correction space has dimension 4109 and its classical map
has rank 633. Every permissible first correction retains all

```text
4109 - 633 = 3476
```

homogeneous choices. The target's first Ward identity annihilates the quantum
image of every classically closed difference of first corrections. Thus none
of these 3476 freedoms can change the second-page target pairing.

Each classical boundary has a closed completion as the full quantum variation
of its predecessor. The explicit product has the closed completion above.
Their leading terms span the entire 3017-dimensional kernel, so the complete
incoming second map vanishes. The conclusion covers all admissible sources
and earlier corrections; it uses no restriction to the tested extension bank.

## What was freshly reconstructed

All 49 literal Albert boundary polynomials and the product were evaluated on
all 3017 published colour assignments: 150850 exact residues. The integer
evaluator uses native matrices, invariant forms, Albert cubic tensors, finite
alternating sums and scalar contractions. All 50 entries agree with the
separate public evaluator at the control points; classical split substitutions
vanish. This checks the full extension bank, but not the historical old
3002-column block or the full determinant.

The outgoing prerequisite was rechecked through its saved matrices, rational
characters and literal Cartan Gram determinants. The earlier optional graph
regeneration invoked a historical colour projector and is excluded from
permitted fresh evidence. See `METHOD_AUDIT.md`.

The further published equality `1862-1434-428=0` excludes the third-order
source, and the Taylor bound `r<=6-3` excludes every higher incoming order.
Together with the first-page Ward identity, the strict closure of the pure
target and the complete second-map zero, these results prove formal
perturbative permanence. The third-source minors were not freshly rebuilt in
this packet. Its independent direct-factor calculation establishes the
complete formal second-map zero using the stated published rank prerequisites.
No physical instanton protection or new positive lift is inferred.
