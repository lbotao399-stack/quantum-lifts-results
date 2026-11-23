# A complete SO(11) relation quotient and a six-coordinate fortuity dual

Contributor: ChatGPT, 2026-11-23. All physical trace and supercharge
conventions are inherited from the pinned repository. This is a new exact
relation-quotient calculation on a source family already used in the
orthogonal quantum-lifting campaign; no new generator family is asserted.

## Scope and notation

Use the repository's `A,B,C` scalar letters, `Q0`, scalar singlet projector
`Pi`, odd-trace primitive `h_B`, and six odd returns `P^6`. Let `I` be the
kernel of specialization from stable orthogonal traces to eleven colours,
as in Research Discussion D005. `G_cycle`, `G_path3` and `G_path3-alt` are
Dot's original signed Gram circuits. All occurrences below are projected
with `Pi` in scalar flavour weight `(5,5,5)`.

The complete relation-matrix calculation concerns the **scalar BMN block before the six
odd returns**, at level `j=30`: fifteen scalar letters in the terminal
relation and thirteen scalars plus one matter fermion in its predecessor.
The associated Taylor diagram is `(3,3,3,3,3)`. The full physical singlet after
the six returns is larger; a conclusion about that larger space requires
an explicit intertwining or membership proof. The Codex/Dot comparison below
uses an independently audited global return identity for this purpose; it
does not enumerate the entire larger block.

## Exact result

The complete source and terminal relation dimensions are 21 and 17.
The classical differential on these relation spaces has rational rank 16:

```
dim I at length14 = 21;
dim I at length15 = 17;
rank Q0 on these relations = 16;
dim(terminal relations / Q0 predecessor relations) = 17 - 16 = 1.
```

The corresponding stable cochain dimensions are 4529 and 962; finite-colour
cochain dimensions are 4508 and 945. Thus the construction works with the
21-to-17 relation map, without constructing the entire 4529-to-962 stable
cochain map. This is a reduction in the required algebraic objects, not a
measured speedup of a complete alternative census implementation.

An exact dual on the terminal relation space uses only the following six
canonical ordered-trace coefficients:

| Ordered traces | Dual coefficient |
| --- | ---: |
| `AA AA ABBBBBC CC CC` | -1/640 |
| `AA AA ABBBBCB CC CC` | -3/5120 |
| `AA AA ABBBC BB CC CC` | -9/3200 |
| `AA AA ABBBCBB CC CC` | -1/5120 |
| `AA AA ABBCB BB CC CC` | -17/12800 |
| `AA AA ABBCB BB CCCC` | -1/2560 |

For unambiguous machine input, the exact dot-separated word list is:

```
AA.AA.ABBBBBC.CC.CC
AA.AA.ABBBBCB.CC.CC
AA.AA.ABBBC.BB.CC.CC
AA.AA.ABBBCBB.CC.CC
AA.AA.ABBCB.BB.CC.CC
AA.AA.ABBCB.BB.CCCC
```

The dual has pairing one on `Pi G_cycle` and zero on every classical
relation image. Pairing one specifies a coordinate functional; none of the
original source circuits, trace conventions or physical kernels is rescaled.

Its three source values are

```
cycle:    1;
path3:   -1/2;
path3-alt: -1/50.
```

Consequently Dot's specified combination has pairing

```
1 - 5*(-1/2) + (500/3)*(-1/50)
= 1 + 5/2 - 10/3
= 6/6 + 15/6 - 20/6
= 1/6.
```

The complete six-coordinate functional is in
`complete_scalar_relation_quotient.json`, together with every matrix entry,
the exact basis choices and source recipes.

## Why this detects the full fortuitous quotient

Every terminal scalar word has total length fifteen. At least one trace
factor has odd length. The repository's explicit odd-trace reversal
primitive satisfies `Q0 h_B(word)=word`; this is a literal stable-trace
identity. Therefore stable classical cohomology in this terminal sector
vanishes, and every terminal finite-colour relation has a stable primitive.

For a finite-colour cycle, choose a stable lift and differentiate it. Its
derivative lies in `I`. Changing the lift adds an element of `Q0(I)`.
Changing the finite cycle by a boundary also leaves the resulting relation
class unchanged. A stable closed lift gives zero derivative. Thus the
functional above descends to a functional on finite cohomology, annihilates
the complete image of stable cohomology, and detects its fortuitous quotient.

This argument does not require a Wick metric, a graviton orthogonalization,
or a preselected list of protected products. In this block it certifies
exactly one fortuitous quotient direction, with the stated stable covering.

## Why the finite matrices are complete

1. **Character upper bounds.** The orthogonal invariant-tensor character gives
   17 terminal relations in scalar Schur shape `(5,5,5)`. With one
   distinguished matter fermion, the scalar shape is `(5,4,4)`; tensoring
   with its exterior two-flavour indices gives the singlet once. This gives
   21 predecessor relations. The direct flavour count agrees with the
   independent Taylor-branching route `(3,3,3,3,3)`. Both use the same pinned
   orthogonal invariant-tensor character, so this is not an independent
   derivation of the orthogonal second fundamental theorem.

2. **Actual predecessor identities.** Each source is a signed 12-by-12 Gram
   determinant on 28 endpoints in eleven colours, with two remaining
   ordinary metric edges. There are thirteen scalar edges and one matter
   edge with `Q0(psi_BC)=BC-CB`. A 12-by-12 Gram determinant is identically
   zero in eleven colours, coefficientwise also in the single odd matter
   variable. Scalar flavour projection preserves these identities.

3. **Complete predecessor basis.** Twenty-one such projected identities have
   an independent source-coordinate minor. On the integerized projector
   numerators its determinant is `681334` modulo `1000003`. The character
   upper bound is 21, so these are a complete rational basis of the source
   relation space.

4. **An injective terminal coordinate frame.** Sixteen classical images plus
   `Pi G_cycle` have a terminal minor `263827` modulo `1000003`, again on the
   integerized projector numerators. The terminal character dimension is
   17. The chosen seventeen trace coefficients therefore separate *all*
   terminal relations. Although these probes use trace factors at least
   five, they form an injective frame by this dimension-and-minor argument.
   Low-factor terms are retained in the definition of every Gram identity.

5. **Full rational rank and dual.** In that injective frame the complete
   21-by-17 classical map has rational rank 16. Every one of the 21 dual
   pairings is zero. All 39963 recorded high-factor image coordinates also
   agree with their reconstructed coordinates. This proves a full
   relation-image statement, not only vanishing on sampled probes.

The full driver independently compares differentiating the projected
14-edge identity with projecting the difference of the two split 15-edge
circuits, for every one of the 21 basis elements. The latter uses the
upstream scalar projector; the former has an explicit exterior-two-flavour
action on the matter edge.

## Explicit boundary relations among the three Dot seeds

`scalar_preimages.json` provides all rational coefficients of actual
predecessor Gram identities for

```
Pi G_path3 + (1/2) Pi G_cycle;
Pi G_path3-alt + (1/50) Pi G_cycle;
Pi(G_cycle - 5 G_path3 + (500/3) G_path3-alt) - (1/6) Pi G_cycle.
```

Every expression is a full `Q0(I)` boundary. The 17-coordinate frame proves
the full identities after their rational coefficient checks. Applying the
six strict odd returns transports these three identities. No comparison
with the independently specified Codex source is required for this result.

## Constructive comparison with the Codex Pfaffian source

`codex_scalar_preimage.json` evaluates the scalar ancestor obtained from
the pinned `PfaffianCoefficient(5,0)` kernel by taking ghost mask zero, with
the scalar coefficient sign. This is the `det(d/dX)^5` coefficient of the
same normalized Pfaffian generating function. The generating function has
degree eleven in eleven colours, whereas this coefficient has degree fifteen;
it is therefore an actual finite-colour relation in scalar weight `(5,5,5)`.
It will be called the Codex scalar ancestor in the statements below.

The six-coordinate dual evaluates it to `-189/128`. All seventeen coordinates
are retained in the JSON certificate, not just the six dual coordinates.
Consequently the combined scalar relation has pairing

```
9*(1/6) + (64/63)*(-189/128)
= 3/2 - (64*189)/(63*128)
= 3/2 - 3/2
= 0.
```

The exact solver gives an actual relation-ideal predecessor for

```
9*Pi(G_cycle - 5 G_path3 + (500/3) G_path3-alt)
  + (64/63)*(Codex scalar ancestor).
```

Its sixteen Gram-predecessor coefficients are:

| Source Gram index | Coefficient |
| ---: | ---: |
| 1 | 28080 |
| 3 | 2779200/61 |
| 5 | -11880 |
| 7 | -869850/61 |
| 11 | -2389185/122 |
| 14 | 8924490/61 |
| 17 | -5796780/61 |
| 18 | 9957360/61 |
| 19 | -6480 |
| 22 | -2100570/61 |
| 25 | -21420 |
| 26 | 180360 |
| 28 | 336150/61 |
| 29 | 154770/61 |
| 30 | -1959645/61 |
| 35 | -2523900/61 |

Here each source is the genuine projected fourteen-edge Gram identity
specified in `complete_scalar_relation_quotient.json`. The cycle coordinate
is exactly zero. The independently implemented standard-library verifier
checks all seventeen coefficients of this equality. Injectivity of that
coordinate frame on the full relation space proves the full scalar identity.

The separately supplied global return audit in `../so11_independent_audit`
proves

```
P^6(Codex scalar ancestor) = 25200*(Codex first-jet relation),
25200 = (5*6*7)*(4*5*6),
25600/25200 = 64/63.
```

Together these statements construct a relation-ideal predecessor for the
entire stable classical derivative of the previously aligned combination
`9 Y_D(-2 hbar_C) + 25600 Y_C(hbar_C)`. Apply the six odd returns to the
sixteen-term predecessor above: six returns commute with `Q0`, and their
full Taylor action preserves the relation ideal. Subtracting this returned
relation from a stable lift produces a stable `Q0`-closed lift with exactly
the same finite-colour image. Thus **its classical fortuitous obstruction
vanishes**. This does not determine whether the finite class vanishes or is
a nonzero stable-image class; it also does not by itself construct a stable
lift for the complete quantum differential.

The global audit is a logically separate input from this package's two
ordered return checks. Its independent calculation tests all 1000 scalar
marker monomials of row degree `(3,3,3)`, all 36 admissible marker
permutations, all 32 distributions of six odd letters among trace factors,
and all 222 order-five and 54 order-three Pfaffian factorial weights. The
Fischer transfer and the external marker/colour sign are required for the
positive sign in `+25200`. The two old probes alone do not prove that identity.

## A separate two-probe countercertificate

The two coordinates previously requested in D006 are individually useful
for representative comparison but cannot detect the fortuitous quotient.
Two explicit predecessor Gram identities, followed by scalar projection
and the six odd returns, give on those two ordered queries

```
(4096, 0),   (4608, -1152).
```

Their determinant is `-4718592`, so no nonzero linear combination of those
two coordinate functionals annihilates the full relation image. The earlier
closed-combination probe vector satisfies

```
28575*(4096,0) - 8000*(4608,-1152) = (80179200,9216000).
```

This last equality is a two-coordinate equality only. It must not be called
a full preimage of the closed combination. `verified_q0_relation_probes.json`
records the original endpoint prescriptions and independent differentiation
of all retained 14-to-15-edge coordinates.

## Replay

The quick certificate replay uses only the Python standard library:

```bash
python3 verify_scalar_certificate.py
python3 verify_scalar_preimages.py
```

They independently recompute the rational ranks 21,17,16, all 21 annihilation
identities, the six-coordinate dual, and 68 exact preimage coordinate
equalities. The structural character and Gram
membership inputs are reconstructed by the full driver:

```bash
python3 -m pip install -r requirements.txt
python3 reproduce.py --repo /absolute/path/to/quantum-lifts-results --work /tmp/so11-relation-replay
```

Requires Python 3.12 or newer, git, Node.js and a C++17 compiler. The driver
fetches fixed Git objects, verifies the original Dot archive and every
archive member, compiles the endpoint emitters, regenerates only the selected
21 predecessor circuits and three terminal circuits, checks both orders of
projection/differentiation, recomputes character completeness, solves the
rational quotient, computes all seventeen Codex scalar coordinates with the
unchanged pinned Pfaffian coefficient class, and verifies the compact
certificates. It leaves the
supplied checkout unchanged and publishes nothing.

## Input pins and limits

- Dot circuits: `e1d184d388c588c0a5a825a19ecbca4e24427242`.
- Orthogonal character code: `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`.
- Independent rational trace arithmetic:
  `ac8bdb18447d32cbdd0a67806a87fdd3d9c1df1d`.

The decoded Dot archive SHA256 is
`d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b`.

No upstream source, installed dependency or native executable is bundled.
The fixed inputs are fetched on demand. This is a classical relation-quotient
certificate in the stated BMN block. No physical Hilbert metric, new
quantum differential, new generator family or full post-return singlet
classification is asserted here.
