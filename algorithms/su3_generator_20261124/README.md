# Explicit rational SU(3) outgoing generator

This package supplies the representative, complete rational dual, and complete
finite-colour relation test requested by **B05** of the frozen research roadmap.
The input is the complete relative holomorphic jet complex of SU(3), at flavour
degrees `(3,0,0)`, dotted degrees `(4,4)`, dotted spin zero, charge conjugation
odd, and length seven. The level is 30. The doubled native physical charge tuple
`(2J1,2J2,2q1,2q2,2q3)` is `(4,4,14,8,8)`.

The normalized operator has **12,224 explicit trace terms**, using **848** of
the **1,896** finite-colour basis cochains. Its rational dual has **349** nonzero
coordinates. Every coefficient is printed in
[`certificate/generator.json.gz`](certificate/generator.json.gz); its trace-word
dictionary and complete cochain bases are printed in
[`certificate/basis.json.gz`](certificate/basis.json.gz).

The exact results are

```text
1896 - 825 - 1037 = 34            complete classical cohomology
34 - 33 = 1                      quotient by every proper product
4991 - 1896 = 3095               all source finite-colour relations
6609 - 1588 = 5021               all target finite-colour relations
4991 - 3348 = 1643               complete stable classical kernel
1643 - 858 = 785                 stable closed source relations
3095 - 785 = 2310                images of every source relation
(3095 + 859) - 1643 = 2311      target relations hit by the stable differential
2311 - 2310 = 1                 complete fortuitous relation quotient

dual(representative) = 1
dual(every classical boundary) = 0
dual(every proper product) = 0
```

The complete relation test establishes fortuity independently of ring
primitivity. See [PROOF.md](PROOF.md) for the exact completeness and rank argument.

## Complete spaces

| Length | Stable dotted singlets | SU(3) dotted singlets | All finite-colour relations |
|---:|---:|---:|---:|
| 6 | 2,091 | 1,169 | 922 |
| 7 | 4,991 | 1,896 | 3,095 |
| 8 | 6,609 | 1,588 | 5,021 |

At the target charges, three scalar jets of flavour one and four fermionic
superfield jets exhaust the allowed letters. Every allowed derivative
distribution is retained. The highest-weight basis is the exact integer kernel
of the raising equations in the complete oriented trace-word space.
Independent symmetric-group characters establish both columns of dimensions.

The native spinor fields retain their dotted indices. The integer pair attached
to a letter in the machine dictionary records its divided Taylor coefficient;
it is a storage coordinate for the symmetric dotted-spinor jet. Both dotted
raising and dotted lowering annihilate the full printed operator exactly.
Flavour raising annihilates it because only the highest flavour occurs.

## Normalization and dictionary

The letter entry `[a,b,mask]` records the powers of the two dotted coordinates
and the exterior flavour mask. Only `mask=0` and `mask=1` occur in this complete
charge block. Mask zero is fermionic; mask one is the scalar of flavour one.
The undifferentiated mask-zero letter is absent in the relative complex.

The printed coefficients use the published divided-jet convention
`Q(Psi)=Psi^2` in
[`fortuity/full/letters.py`](../su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/full/letters.py).
An undivided jet equals `a! b!` times the corresponding printed letter.
After this invertible factorial change, the differential here is minus the
undivided `Q0` in the pinned B04 native driver. The operator, finite closure,
nonexactness and all ranks are unchanged by this global differential sign.

Each sparse coefficient is stored as `[coordinate,numerator,denominator]`, with
a positive denominator. A word is an ordered product of oriented traces. Graded
cyclic signs and graded signs between traces are part of the dictionary.
`maps.json.gz:7.source_basis` selects and orders the 1,896 basis cochains on
which the dual is defined. The complete finite coordinate relations make this
a functional on the entire finite-colour cochain space.

The common denominator of the printed operator is **2,687,385,600**; the resulting
integer coefficients have greatest common divisor one. Its stable derivative
contains **39,829** trace terms. In the complete target it is exactly zero.
[`closure_relation.json`](certificate/closure_relation.json) expresses this
derivative in **2,047** of the complete target relation basis vectors.

## Proper-product coverage

The complete product classification is the immutable B04 input at
[`a90728a6ca04d74f5475adae4bf37c8d94f34ac0`](https://github.com/lbotao399-stack/quantum-lifts-results/blob/a90728a6ca04d74f5475adae4bf37c8d94f34ac0/algorithms/su3_j30_classical_ranks/instanton_probe/outgoing_channels.md):
**56 factor pairs, 112 factor sectors, 114 classical matrices**, and a
33-dimensional proper-product span. Its raw minor archives are retained in
[`certificate/b04/`](certificate/b04/), with exact input hashes in the manifest.
Those historical minor entries are inherited evidence; this package does not
claim to have rebuilt those 114 matrices.

All 33 proper-product representatives have been reconstructed here as actual
closed trace polynomials, with exact symmetry checks. The full dual annihilates
every one. In addition, all 447 stable graviton recipes in dotted weight `(4,4)`
were enumerated. Their exact dotted-singlet equations produce 86 closed
combinations. These 86 combinations together with the complete stable incoming
map exhaust all 1,643 stable cycles. The dual annihilates their entire finite
image. The number 86 counts an overcomplete collection of closed combinations.

## Verification

Dependencies are Python 3.12, a C++17 compiler, and
[`requirements.txt`](requirements.txt). From the repository root:

```sh
python3 -B algorithms/su3_generator_20261124/verify.py --rational-ranks
python3 -B algorithms/su3_generator_20261124/independent_verify.py
```

The first command validates the archived exact certificate and all decisive
ranks; `--rational-ranks` additionally recomputes the two adjacent ranks directly
over the rationals. It verifies archived frame minors and states that their
native entries are reused. The second command imports none of the construction
algebra or frame modules. It independently differentiates the entire printed
operator, checks both spin actions and charge conjugation, and checks the full
residual at every one of the 1,588 target coordinates with integer arithmetic.
It also checks all 1,169 incoming columns, 33 proper-product columns, 86 stable
graviton columns, and the normalized pairing one.

To reconstruct the complete source data and every new native entry in a fresh
directory:

```sh
python3 -B algorithms/su3_generator_20261124/reproduce.py /tmp/su3-generator-replay
```

The directory must be empty. This command enumerates every trace cochain,
reconstructs all 9,038 finite relations in 54 complete frames, reconstructs both
full native differential matrices, extracts the vector and dual, and compares
all nine core mathematical artifacts byte for byte with the published payloads.
It finishes by running both verifiers. The measured execution receipts are in
[`certificate/`](certificate/).

## Scope

This is an exact classical cohomology, ring-generator and finite-rank fortuity
certificate. Its superconformal parent, radial-adjoint representative, physical
inner product and renormalized instanton matrix elements remain separate tasks.
The calculation makes no new higher-order or nonperturbative lifting claim.

