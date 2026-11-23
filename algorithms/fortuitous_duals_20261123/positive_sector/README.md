# Sp(3) omega: a complete rational trace-identity quotient

This is a positive test of fortuitous detection on the actual USp(6) class
omega at level j=28. The source conventions are unchanged: weight
`(3,3,5,0,0)` lists the two spacetime Taylor degrees followed by the three
odd-coordinate degrees; the cochain has length nine.

`Q0`, `C_st`, finite-colour specialization and its kernel `I` have the meanings
of Research Discussion D005. This is the formal orthogonal/symplectic trace
cover used by the pinned source. The calculation does not choose a Wick
metric and does not identify a physical observable.

## Added results

The complete formal predecessor and centre spaces have dimensions2624 and
2870. Their classical differential ranks are1672 and1198. Nonzero minors
modulo1000003 have determinants111063 and452105. The original differential
squares to zero on all5494 input rows. Therefore, over the rationals,

```
dim H(Q0,C_st) = 2870 - 1672 - 1198 = 0.
```

This proves that the full universal closed image in the finite cohomology
vanishes here. It does not assume that a selected graviton list is complete.

The full finite-colour evaluation map has rank2689, with minor23258 modulo
1000003. Its kernel consists of181 identities. These were reconstructed as
rational combinations of the2870 formal basis vectors and then certified
exactly by the bounded-integer calculation below. Their classical images
have rank116; the finite outgoing map has rank1081. Thus the connecting-map
rank calculation reads

```
1198 - 116 - 1081 = 1.
```

The existing rational omega supplies that one-dimensional class. The new
result is a complete rational identity-space certificate and a usable dual
construction; this is not a newly discovered BPS state.

The companion enumeration calculation constructs a sparse rational
functional on the formal derivative space which annihilates all181 identity
images and evaluates to one on the formal derivative of the unchanged
rational omega. See the parent deliverable's enumeration directory for the
explicit40 coefficients and its independent factorization checks.

## Full native regeneration

`build_native_frames.py` regenerated:

- the complete2870 by2689 finite-colour evaluation table;
- all2624 predecessor images in the same2689-coordinate frame;
- the complete2689 by1281 outgoing table;
- all2624 pairings of the archived target functional with predecessor images;
- the archived target-functional value one on omega;
-480 independent dense Python trace coordinates.

The finite outgoing minor has rank1081 and determinant110010; the incoming
minor has rank1607 and determinant624328, all modulo1000003. The original
minor hashes agree exactly. The complete arrays are saved by the generator;
they may be omitted from a compact distribution because the executable
recipe reproduces them.

## Exact rational identities, not modular reconstruction alone

`rational_relation_basis.json.gz` contains all181 rational coefficient rows
with29192 nonzero entries. Their181 free-coordinate columns form an identity
matrix, proving independence. Row denominators are among

```
2,4,6,8,12,16,24,32,48.
```

After clearing the row denominators and expanding the complete formal
highest-weight polynomials, the181 identities contain126900 trace terms.
The largest coefficient L1 norm is216596.

A fresh2689-coordinate finite-colour frame uses the same small integer
symplectic matrices at every modulus. Every entry has absolute value at
most2. Exact symplecticity is checked for every matrix. A trace monomial has
nine matrix entries and four exterior variables. For each exterior
coefficient, at most4! assignments contribute; the complete matrix-index
sum has at most6^9 terms. Thus a word has absolute value at most

```
4! * (2*6)^9 = 123834728448.
```

Every expanded identity evaluation has absolute value at most

```
216596 * 123834728448 = 26822106842923008.
```

The full finite-frame minor is162437 modulo1000003, so this frame is
injective on the entire2689-dimensional rational finite cochain space.
All181*2689=486709 identity evaluations vanish modulo each of the three
primes1000003,1000033,1000037 at identical integer points. Their product is

```
1000003 * 1000033 * 1000037 = 1000073001431003663
                            > 2 * 26822106842923008.
```

The integer evaluations must therefore be zero exactly. Character
saturation makes the identities literal finite-colour identities, and
2870-2689=181 proves that the complete rational kernel has been obtained.
There are480 additional independent integer/native trace checks at each
prime,1440 in total. The procedure and all exact bounds are in
`certify_rational_relations.py` and `rational_relations_result.json`.

`exact_relation_image.py` additionally eliminates the complete rational
identity-image matrix over the integers, supplying65 independent null
relations and a nonzero116 by116 integer minor. This proves the image rank
116 over the rationals directly.

## Original rational omega

The source provides819 rational coefficients in its finite basis, with
common denominator13824. `rational_omega.py` expands the same cochain to
4264 terms without changing its normalization. Its formal classical
derivative contains4502 terms. Exact integer evaluation on both archived
complete1281-coordinate frames gives zero in each coordinate. Their frame
minor determinants are878791 and784798 modulo1000003.

The complete cochain is in `omega_rational.json.gz`; the receipt is
`rational_omega_result.json`. This is an independent replay of the existing
rational class, not a new class construction.

## Inputs and reproduction

Repository: https://github.com/lbotao399-stack/quantum-lifts-results

- Source algorithms and omega data:
  `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`,
  `algorithms/dual_singlets_20261118`.
- Native runtime:
  `082be7070d207e2cdecfd0a91b03d5e014fd303d`,
  `algorithms/su3_independent_prime_20261118/runtime`.

Requires Python3.12, git, a C++17 compiler, NumPy, SciPy and python-flint.
From this directory:

```bash
python3 -m pip install -r requirements.txt
python3 reproduce.py --repo /absolute/path/to/quantum-lifts-results
```

`fetch_inputs.py` obtains353 immutable input files,2603911 bytes, and records
their Git blobs and SHA256 hashes. It does not modify the source checkout.
All calculation stages were run; the portable fetch and import path were
also exercised. The combined wrapper executes those stages in order.

Source character dimensions, finite invariant completeness and the
definition of the classical differential remain pinned upstream inputs.
The exact181-row identity proof, stable acyclicity calculation and full
evaluation/dual tests are the computations added here. No quantum
differential or physical harmonic projection is computed in this directory.
