# Native letters, complete scalar identities and actual cochains

## Physical letters and Q0

The scalar labels are the native Taylor letters

```
A = Psi_(0,0,(0)), B = Psi_(0,0,(1)), C = Psi_(0,0,(2)).
```

The internal cyclic matter labels used by the coefficient engine are

```
a = Psi_(0,0,(1,2)),
b = -Psi_(0,0,(0,2)),
c = Psi_(0,0,(0,1)).
Q0(a)=BC-CB, Q0(b)=CA-AC, Q0(c)=AB-BA.
```

The `--native-sorted-theta` query option converts the middle sign.
`NATIVE_CONVENTIONS.json` records a direct comparison to the pinned native
`Q_letter`. This is a change of coefficient encoding only; the native Q0
and physical letter normalizations are unchanged.

Words denote ordinary defining traces. A dot separates trace factors.
Cyclic trace identification and orthogonal reversal are imposed exactly.
Sources have one odd matter letter, so no hidden odd-odd permutation
occurs in the source coefficient arithmetic.

## Finite Pfaffian and universal trace expansion

The auxiliary matrices are

```
H=[[1,0],[0,-1]], E=[[0,1],[0,0]], F=[[0,0],[1,0]],
J=[[0,1],[-1,0]].
```

For scalar-marker rows X_A,X_B,X_C and odd matter-marker rows Y_a,Y_b,Y_c,
in component order H,E,F, use

```
Pf(J tensor 1_N + sum_(i,s) (J M_s) tensor
   (x_(i,s) phi_i + eta_(i,s) psi_i)) / Pf(J tensor 1_N).
```

The odd marker is placed before its odd matrix. Only coefficients linear
in the odd markers are needed. Each J M_s is symmetric, so the displayed
matrix is skew. The defining Pfaffian has physical degree at most N.
Its degree-sixteen and degree-seventeen coefficients therefore vanish
at N=15.

The universal coefficient series is

```
exp((1/2) Tr_(aux,color) log(1 + sum M_s tensor (...))).
```

The coefficient engine expands this in ordinary cyclic/reversal trace
words, retaining the exact orbit factors and factorials for repeated
trace factors. Marker coefficient extraction uses ordinary derivatives;
all derivative factorials are retained.

## Explicit marker differential operators

The notation in this section concerns auxiliary polynomial differential
operators only. It defines no new physical normalization.

```
D = det(X_A,X_B,X_C),
dot(U,V) = U_H V_H + 2 U_E V_F + 2 U_F V_E,
Z_H = 2(X_A,E X_B,F - X_A,F X_B,E),
Z_E = X_A,H X_B,E - X_A,E X_B,H,
Z_F = X_A,F X_B,H - X_A,H X_B,F,
AA=dot(X_A,X_A), ZZ=dot(Z,Z),
T=dot(X_A,Y_a)+dot(X_B,Y_b)+dot(X_C,Y_c).
```

Replace X and Y by their marker derivatives and then set markers to zero.
For (7,7,3) the terminal relation and complete predecessor basis are

```
terminal: D^3 ZZ^2,
source1: D^2 T ZZ^2,
source2: D^2 [sum_i dot(Z,X_i)dot(Z,Y_i)] ZZ.
```

For (9,5,3) they are

```
terminal: D^3 AA^2 ZZ,
source1: D^2 T AA^2 ZZ,
source2: D^2 [sum_i dot(X_A,X_i)dot(X_A,Y_i)] AA ZZ,
source3: D^2 [sum_i dot(Z,X_i)dot(Z,Y_i)] AA^2.
```

The terminal relations have seventeen scalar letters. The predecessors
have fifteen scalars and one matter fermion, hence sixteen letters.
The verifier expands the complete marker polynomials and checks their
weights and both flavour raising actions exactly, including the
contragredient matter action.

## Character completeness

At degree sixteen in fifteen colours the only missing defining-vector
invariant diagram is (2^16). At degree seventeen the missing diagrams
are (2^17) and (4,2^15). Their transposes are (16,16), and (17,17),
(16,16,1,1).

For an adjoint permutation with m_k cycles of length k, choose t_k crossed
two-slot cycles. The defining-vector permutation has 2(m_k-t_k) cycles
of length k and t_k cycles of length 2k. Transposition of the diagram
cancels the antisymmetric-adjoint projector sign. The exact relation
character is the sum over these choices with weight

```
product_k binomial(m_k,t_k) / 2^(sum_k m_k).
```

Exact Murnaghan-Nakayama characters and character inner products give
terminal multiplicities 3 and 4. Exterior-two-flavour Pieri removal for
one matter letter gives

```
(7,7,3): scalar shapes (6,6,3),(7,6,2), each multiplicity1;
(9,5,3): scalar shapes (8,4,3),(8,5,2),(9,4,2), each multiplicity1.
```

The actual predecessor images have ranks 2 and 3 on the displayed
coefficient probes. They consequently exhaust the full source spaces.
The complete terminal dimension is saturated by adjoining the specified
scalar Pfaffian. No trace-factor cutoff is used in the identities.

## Explicit equivariant stable primitive

Every terminal word has odd total length seventeen and thus contains an
odd-length trace. Reverse that trace by adjacent swaps in the fixed
sequence `end=L-1,...,1`, `i=0,...,end-1`. Replace each swap difference
by its matter commutator primitive, with a factor one half. The sum
telescopes to `(word-reversed_word)/2`. Orthogonal odd reversal gives
exactly the original trace. Therefore Q0 of this explicit raw primitive
is the original word in the stable ordered-trace algebra.

The delivered `primitive()` averages over odd trace factors, all cyclic
rotations of the chosen trace and its antisymmetric reversal. It is
independent of trace presentation, retains Q0 primitive(word)=word, and
intertwines flavour transformations. The verifier checks both the
telescoping identity and raising/lowering intertwining on the probe words.
Applying it coefficientwise to the terminal Pfaffian gives the actual
highest-weight cycle used in this package.

The inverse coefficient oracle is complete. A source trace containing
one matter letter can only arise by fusing its complementary scalar pair
at one of the fixed reversal swaps. Enumerate every cyclic/reversed
presentation of that source trace, every compatible swap position, and
both ordered scalar pairs. Undo the recorded permutation of positions,
then canonicalize and deduplicate. Every contributing scalar word is
obtained by reversing its specific fusion step. The oracle evaluates
the complete forward primitive on each candidate, so duplicate
presentations do not duplicate coefficients. No arbitrary trace-word
cutoff is involved.

## Connecting functional

The stable derivative of each delivered cochain is its exact scalar
Pfaffian. It becomes zero in fifteen colours, giving a finite Q0 cycle.
Changing a stable lift adds a relation in I16; its derivative changes
the scalar relation by Q0(I16). A finite boundary also leaves this class
unchanged. A stable closed lift has zero derivative. Therefore each
displayed dual defines a functional on finite cohomology, annihilates
the entire stable cohomology image and takes value one on the specified
cycle. This proves both non-exactness and fortuity.

The dotted Taylor degrees are nonnegative and the representatives have
zero dotted weight. An additional boundary from outside this direct
summand of the full relative jet complex cannot remove either class.

The lower-rank contraction engine reuses the finite-Pfaffian tile
construction at commit78d170a26afc2a16cae43f123df86d5880be6196,
`reviews/orthogonal-mixed-lifts/bmn13/bmn_relation_probe.cpp`, with its
flavour probe bank extended by the stated AB wedges. Its outputs are
used to construct new exact lower-rank preimages, not to recount the
old product classes as discoveries.
