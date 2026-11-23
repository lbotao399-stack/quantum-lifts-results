# A closed two-probe dual for minimal odd-orthogonal Pfaffian classes

## Scope and unchanged conventions

This note uses the scalar lane's existing notation and normalization: `D`,
`AA`, `ZZ`, `PAIR`, the auxiliary matrices `H,E,F`, the cyclic matter letters,
and `Q0 psi_A=[B,C]`, `Q0 psi_B=[C,A]`, `Q0 psi_C=[A,B]`.
`AUXILIARY_CONVENTIONS.md` specifies them completely. No physical field,
supercharge, Pfaffian, or representative is rescaled here.

Let the existing highest-weight label be an all-odd partition
`lambda1 >= lambda2 >= lambda3 >= 3`. The defining orthogonal dimension is
`N=lambda1+lambda2+lambda3-2`, which is odd. The source has length `N+1`,
with `N` scalars and one matter fermion; its scalar defect has length `N+2`.
All assertions concern classical Q0 and the trace covering. No statement
about survival under Q1 or higher corrections follows.

The complete target marker already implemented by `operators(shape)` is

```
D^lambda3 * AA^((lambda1-lambda2)/2) * ZZ^((lambda2-lambda3)/2).
```

It acts by ordinary differentiation on the degree-N+2 term of the original
Pfaffian generating series. That term is zero in SO(N). Apply the explicit
equivariant odd-trace primitive in the scalar code to obtain the finite
closed cochain. Its universal Q0 image is precisely the specified defect.

## The two probes and closed coefficients

Read the coefficients of the following two canonical multitraces:

```
Tr(ABC)   Tr(AA)^((lambda1-1)/2) Tr(BB)^((lambda2-1)/2) Tr(CC)^((lambda3-1)/2),
Tr(ABC)^3 Tr(AA)^((lambda1-3)/2) Tr(BB)^((lambda2-3)/2) Tr(CC)^((lambda3-3)/2).
```

The first coefficient of the scalar Pfaffian defect is exactly

```
(-1)^((lambda1-lambda3)/2)
* (lambda1+2)!! * (lambda2+1)!! * lambda3!!
/ 4^((lambda3-1)/2).
```

The second is the first multiplied by

```
2*(lambda1+1)*lambda2*(lambda3-1)/3.
```

Both are nonzero. Assign dual coefficient `3/(2*first)` to the first probe
and `-1/(2*second)` to the second probe. The pairing with the defect is
`3/2-1/2=1`.

For each source relation, its two Q0 coefficients are its respective
defect coefficients multiplied by the same constant and by the numbers
`1` and `3`. Its dual pairing is therefore zero. The derivation below proves
this factor of three without solving a cochain kernel or a large matrix.

## Exact auxiliary identities

For this calculation only, use the source code's exponent parameters:
`a=(lambda1-lambda2)/2`, `b=(lambda2-lambda3)/2`, and `c=lambda3`.
The determinant differential operator on the nine scalar marker variables
obeys

```
det(d/dX)[D^c AA^a ZZ^b]
 = c*(c+2*b+1)*(c+2*b+2*a+2) * D^(c-1) AA^a ZZ^b.
```

The three source marker recipes, omitting the second if `a=0` and the third
if `b=0`, are

```
D^(c-1) AA^a ZZ^b * PAIR,
D^(c-1) AA^(a-1) ZZ^b * sum_i dot(X_A,X_i)dot(X_A,Y_i),
D^(c-1) AA^a ZZ^(b-1) * sum_i dot(Z,X_i)dot(Z,Y_i).
```

Contract one Y derivative with the covariant auxiliary metric and one X
derivative with the same physical flavour index. Their divergences are
`D^(c-1) AA^a ZZ^b` multiplied respectively by

```
lambda1+lambda2+lambda3+6,  lambda1+4,  lambda3.
```

`universal_auxiliary_identity.py` is a finite exact certificate of these
parameter-uniform identities. After division by the common marker
polynomial, repeated Leibniz expansion has total degree at most three in
`a,b,c`, with coefficients in the rational-function field of X. The code
checks all twenty full-polynomial identities on the simplex of parameter
nodes with `a+b+(c-1)<=3`. Its 20-by-20 evaluation determinant is exactly
`110592`, so these checks imply the identity for all parameters. This is
polynomial interpolation with a proved degree bound, not an extrapolation
from selected physical ranks. The divergence identities have degree one
and are certified on four affine-independent parameter nodes each.

## Why the actual Q0 coefficients have the required trace-count factor

In the trace-log expansion, the coefficient of `Tr(X_i X_i)` is
`-dot_cov(X_i,X_i)/2`, that of `Tr(ABC)` is `D(X)`, and that of
`Tr(X_i psi_i)` is `-dot_cov(X_i,Y_i)`.

For a probe with an odd number `r` of `Tr(ABC)` factors, a Q0 predecessor
must replace one cubic trace by `Tr(X_i psi_i)`: a quadratic trace with
identical scalar letters has zero commutator predecessor.

```
Q0 Tr(A psi_A) = Tr(ABC)-Tr(ACB) = 2 Tr(ABC),
Q0 Tr(B psi_B) = Tr(BCA)-Tr(BAC) = 2 Tr(ABC),
Q0 Tr(C psi_C) = Tr(CAB)-Tr(CBA) = 2 Tr(ABC).
```

The predecessor has `1/(r-1)!` in the exponential expansion, whereas the
target has `1/r!`. Their ratio is `r`. Moving the scalar parameter
multiplication from `dot_cov(X_i,Y_i)` onto the marker polynomial gives
exactly the divergence above; its sign is minus. Moving one multiplication
by D onto the target marker gives the determinant-derivative identity.
Consequently, the Q0 coefficient divided by the scalar-defect coefficient
at that probe is

```
-2*r * [the corresponding divergence multiplier]
/ [lambda3*(lambda2+1)*(lambda1+2)].
```

Taking `r=1,3` proves both boundary cancellation and the normalization-one
pairing claimed above, in the original physical conventions.

## Derivation of the nonzero coefficient formula

Move the single D factor of the first probe onto the target marker by the
determinant-derivative identity. The remaining quadratic exponential is
the moment-generating function of three independent formal Gaussian
vectors with covariance minus the covariant auxiliary metric.

The original inverse auxiliary metric has determinant `-4`.
Changing Gaussian integration coordinates to a Euclidean metric gives
the following factors, without changing any physical operator:

```
D^(c-1) -> 4^(-(c-1)/2) * (det of the real Gaussian matrix)^(c-1),
AA^a    -> (-1)^a * (first Gaussian row squared norm)^a,
ZZ^b    -> (-1)^b * (first-two-row Gram determinant)^b.
```

Orthogonally resolving the three independent Gaussian rows leaves three
independent squared lengths with 3, 2 and 1 real components. Their required
moments are respectively

```
lambda1!!,  (lambda2-1)!!,  (lambda3-2)!!.
```

Multiplication by the determinant-derivative factor
`lambda3*(lambda2+1)*(lambda1+2)` produces exactly

```
(-1)^((lambda1-lambda3)/2)
* (lambda1+2)!!*(lambda2+1)!!*lambda3!!
/ 4^((lambda3-1)/2).
```

For any odd `r<=lambda3`, moving r factors of D gives a falling product
of the same determinant factors and `1/r!`. Increasing r by two changes
the coefficient by

```
4*(lambda1-r+2)*(lambda2-r+1)*(lambda3-r)/((r+1)*(r+2)).
```

This yields the displayed second-probe coefficient at `r=1`.
`closed_dual.py --verify` independently evaluates the actual trace
coefficients and actual Q0 predecessors, including 19 highest weights and
50 odd trace counts, and agrees exactly with these formulas.

## The completeness gate for a cohomological statement

The dual only proves global nonboundary if the displayed source relations
exhaust the relevant highest-weight part of the first relation module.
The invariant-theory argument is as follows; its source-completeness and
parity conditions were checked separately from the coefficient calculation.

At the first orthogonal relation degree `N+1`, the only omitted metric
tensor diagram is `(2^(N+1))`. Transpose it and cancel the adjoint
antisymmetrizer sign. Schur-Weyl duality identifies this relation functor
with invariants of the three-dimensional symmetric square of the
two-dimensional auxiliary representation, hence auxiliary SO(3) vector
invariants. The concrete Pfaffian coefficient map realizes this
identification, rather than merely having the same character: its image is
in the finite relation space, and its coefficients on products of
quadratic traces pair the auxiliary marker with products of Gram
quadratics. These Gram products span the whole even-degree invariant
space. The differential-coefficient (Fischer) pairing restricted to that
space is nondegenerate: after a complex invertible change of auxiliary
basis it is the positive Euclidean pairing on compact O(3) invariants.
This proves injectivity of the concrete Pfaffian map. The tensor-functor
dimension equality then proves surjectivity. The same argument applies
with one distinguished matter slot; there is only one odd marker, so
there is no additional odd-odd interchange sign.

The marker degree is even. With three scalar vectors and one matter
marker, SO(3) invariants are therefore O(3) invariants. To first order in
the matter marker they form a free module generated by the nine
`dot(X_i,Y_j)` over the six scalar Gram polynomials. There is no linear
syzygy: three generic scalar vectors are linearly independent. Relations
from a four-vector Gram determinant first occur at matter degree two.

As a flavour module, the generators transform as the product of the
fundamental and its exterior square, namely determinant times End of the
fundamental. The scalar Gram algebra is multiplicity-free with all-even
row diagrams. For an all-odd target lambda, removing the determinant
leaves an all-even diagram. Tensoring by End can return this diagram only
from that same all-even diagram: moving a box between distinct rows makes
two row lengths odd. Its multiplicity is the number of addable corners
in three rows, exactly

```
1 + indicator(lambda1>lambda2) + indicator(lambda2>lambda3).
```

The displayed source polynomials have this number of independent
directions. Independence follows on generic scalar vectors from the
identity endomorphism and, when present, the two rank-one directions along
X_A and along the auxiliary cross vector Z. Thus the two-probe dual kills
the whole source relation space. Combined with the exact odd-trace
primitive and finite Pfaffian identity, this gives at least one rational
fortuitous class for every all-odd lambda with lambda3>=3 at
`N=sum(lambda)-2`.

This is an existence statement. It does not assert that every terminal
relation class or the entire cohomology census has been exhausted. The
SO15 numerical certificates additionally give complete quotient dimensions
one for the two selected weights, and their separate SO13 preimages give
the stronger strict rank-support result.

The corresponding SO tower and brane-fusion interpretation were already
anticipated in the pinned upstream work. The contribution here is the
explicit parameter-uniform separating functional and its proof, rather
than a claim that the physical family itself was unexpected.
