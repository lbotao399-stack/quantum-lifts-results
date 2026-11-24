# Complete physical radial form for the first Sp(2) fortuitous block

## Result and coupling scope

The prospectively selected undeformed four-dimensional N=4 block has gauge group
USp(4), classical dimension 9, left spin 1 at its highest weight, dotted spin 0,
and flavour SU(3) irreducible representation with Dynkin labels (1,0).
Its highest flavour weight is `(R1,R2,R3)=(3,2,2)`.

At the first interacting supercharge coefficient, the complete multiplicity
space has **23 physical harmonic states**: **22 in the stable image and one
fortuitous state**, for the universal symplectic trace covering. The entire
flavour triplet contains **69 states**. Every representative and every entry
of the positive 69 by 69 radial matrix is provided in the files below.

This is the free radial form together with the leading interacting physical
supercharge, hence the physical one-loop zero-energy condition. It does not
assert an exact finite-coupling state, a renormalized higher-order metric, or
the composite-operator comparison required by F03.

| Native length | Native weight `(z+,z-,theta1,theta2,theta3)` | Universal highest-weight dimension | Finite dimension | All finite relations | Outgoing rank | Incoming rank | Harmonic dimension | Stable image | Fortuitous quotient |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | (3,3,1,0,0) | 35 | 35 | 0 | 17 | 18 | 0 | 0 | 0 |
| 6 | (2,2,2,1,1) | 427 | 379 | 48 | 156 | 202 | 21 | 21 | 0 |
| 7 | (1,1,3,2,2) | 816 | 605 | 211 | 230 | 374 | 1 | 0 | 1 |
| 8 | (0,0,4,3,3) | 261 | 149 | 112 | 45 | 103 | 1 | 1 | 0 |
| Total | Complete declared representation | 1539 | 1168 | 371 | 448 | 697 | 23 | 22 | 1 |

The finite dimensions by fermion number are independently reproduced by an
exact root-character calculation:

| Native length | Fermion number and multiplicity |
|---:|---|
| 5 | 4:35 |
| 6 | 2:188, 4:181, 6:10 |
| 7 | 0:59, 2:380, 4:166 |
| 8 | 0:19, 2:130 |

## Exact meaning of completeness

The number 1168 is the **entire multiplicity space of the specified compact
irreducible representations**, not the unrestricted Cartan weight space.
Other flavour representations that contain the weight (3,2,2) have different
flavour Casimirs. They are outside the representation fixed in
`declaration.json` before the metric was calculated. No distinction between
stable-image and fortuitous states is used to define a physical symmetry.

The Lorentz labels are `(left spin,dotted spin)=(1,0)`. The selected supercharge
requires left highest weight 1. All native letters keep their symmetric
undotted and dotted spinor indices; the two integer jet degrees record the
weights of these symmetric spinors. Every allowed derivative, gaugino,
matter fermion, scalar and field-strength letter is included. The states
need not be superconformal primaries: descendants with these charges remain
in the block.

The full flavour representation has the three weights
`(3,2,2)`, `(2,3,2)`, `(2,2,3)`. For every one of the 1539 universal basis
rows, `export_block.py` checks exact lowering and raising identities for the
entire fundamental triplet. It also checks that lowering agrees with an
explicit flavour permutation, with its Grassmann sign. The single-letter
radial form is invariant under this permutation. Consequently every one of
the three diagonal weight matrices is the exported 23 by 23 multiplicity
matrix; different weights have zero pairing by their distinct Cartan
charges. All 69 trace-polynomial states and all 4761 matrix entries are
exported explicitly. The complete finite cochain representation has 3504
states and 1113 finite-rank null relations.

All native lengths with these physical charges are present. To see the
exhaustion directly, equality of the two angular charges gives equal jet
totals. If that common total is `t`, the native weight and length must be
`(t,t,4-t,3-t,3-t)` and `8-t`. Nonnegative letter weights give exactly
`t=0,1,2,3`. These four possibilities are the four rows of the table.
The complete predecessors have the same native weights and one fewer
letter; every predecessor is enumerated. Every outgoing image is evaluated
as a full trace polynomial, with no truncated target basis.

The highest-weight equations are solved simultaneously over the rational
numbers. This constructs kernels of raising maps. No averaging or extremal
operator is used to construct states. `character_check.py` separately counts
the finite dimensions using the eight roots of type C2, the two zero weights
of its adjoint, the Weyl character denominator, and symmetric/exterior
cycle-index coefficients. These character coefficients count dimensions;
they are not operators acting on states.

## Physical letters, coupling and radial conjugation

The native letters are precisely the Taylor coefficients used by the pinned
`fortuity.full.letters` module. Their physical identification is the Taylor
expansion of the Chang--Yin generating field

\[
\Psi(z,\theta)=-i\left[
 \lambda(z)+2\theta_i\phi^i(z)
 +\epsilon^{ijk}\theta_i\theta_j\psi_k(z)
 +4\theta_1\theta_2\theta_3 f(z)\right],
\]

with symmetric covariant derivatives and the generating gaugino specified
in equations (2.7)--(2.11) of [Chang--Yin, arXiv:1305.6314v2](https://arxiv.org/pdf/1305.6314v2).
The dotted spinor contraction in `lambda(z)` is retained. In particular its
two derivative terms and their mixed radial pairing are both included.

The canonical native letter with jet degrees `a,b` is the corresponding
Taylor coefficient multiplied by `(4*pi/g_YM)*2^(a+b)`. Equivalently, take
the coefficient of the same monomial in `(4*pi/g_YM)*Psi(2z,theta)`.
The exact classical identity `Q Psi = Psi^2` therefore gives the first
interacting physical coefficient

\[
 [Q_{\rm physical}]_{\text{order }g_{\rm YM}}
 =\frac{g_{\rm YM}}{4\pi}\,Q_{\rm native}.
\]

This equality includes the action on covariant derivatives. The native
Grassmann signs are those obtained by expanding this generating-field
identity, not by a separate guessed field rule.

For a native coefficient containing `s` distinct flavour variables, the
single-letter radial norm is

\[
 \frac{(a+b+s-1)!}{a!\,b!},
 \qquad a+b+s\geq1.
\]

This is the unified norm (3.13)/(E.30) of
[Chang--Feng--Lin--Tao, arXiv:2306.04673v2](https://arxiv.org/pdf/2306.04673v2)
after the stated canonical rescaling and division by Taylor factorials.
The derivation below independently resolves an inconsistent elementary
field-strength factor in that preprint; the unified formula is not accepted
solely by citation.

The radial bra uses the Hermitian generators with
`tr(T_A T_B)=delta_AB`. In the fundamental four-dimensional representation,
the complete colour covariance is

\[
 \sum_A (T_A)_{ij}\overline{(T_A)_{kl}}
 =\frac12\left(\delta_{ik}\delta_{jl}+J_{il}J_{kj}\right),
 \qquad
 J=\begin{pmatrix}0&1_{2}\\-1_{2}&0\end{pmatrix}.
\]

`independent_checks.py` constructs all ten compact Hermitian generators and
checks every one of the 256 entries of this identity. The matrix is a
physical free-field contraction. It is never replaced by a coefficient
functional on cohomology.

### Independent field-strength normalization

The gauge action (D.9) of arXiv:2306.04673v2, in the stated trace convention,
has coefficient `1/(2*g_YM^2)` multiplying `tr(F_mu_nu F_mu_nu)`. Its
gauge-fixed free propagator is

\[
 \langle A_\mu^A(x)A_\nu^B(0)\rangle
 =\frac{g_{\rm YM}^2}{8\pi^2|x|^2}\delta^{AB}\delta_{\mu\nu}.
\]

Differentiating both fields and applying the two tensor inversions in the
radial bra gives

\[
 \langle F_{\mu\nu}^A\mid F_{\rho\sigma}^B\rangle
 =\frac{g_{\rm YM}^2}{2\pi^2}\delta^{AB}
  (\delta_{\mu\rho}\delta_{\nu\sigma}-\delta_{\mu\sigma}\delta_{\nu\rho}).
\]

Applying precisely (A.5) and (D.4), including normalized
antisymmetrization and `F_{alpha beta}=-(sigma^{mu nu})_{alpha beta}F_{mu nu}/16`,
then gives the covariant spinor result

\[
 \langle F^{\alpha\beta,A}\mid F_{\gamma\delta}^B\rangle
 =\frac{g_{\rm YM}^2}{128\pi^2}\delta^{AB}
   \delta^\alpha_{(\gamma}\delta^\beta_{\delta)}.
\]

All 256 curvature entries and all 16 spinor-pair entries are checked in exact
arithmetic, after radial conjugation and before using the unified native
norm. Thus the primary norm of `f` is `g_YM^2/(128*pi^2)` and the canonical
native field-strength norm is exactly `16*16/128=2`.
The field-strength line in (E.12)/(E.13) of the cited preprint is smaller by
a factor two; combined with (E.29), it would disagree with that paper's own
(3.13)/(E.30). It is not used here.

The scalar and matter-fermion primary norms are
`g_YM^2/(64*pi^2)`. The dotted-gaugino primary norm is the same, with its
Kronecker spinor pairing. The exact radial commutator recurrences (E.21)--(E.23)
give every descendant norm. The checker uses these recurrences independently
for all 48 letters encountered in the complete block and its predecessors,
including the mixed gaugino terms. It recovers the unified factorial norm
above for every letter.

## Full matrix, all nulls and the two physical conditions

`wick.cpp` sums every permitted Wick matching, with its fermion sign, and
both terms in every symplectic colour contraction. Its signed index
identifications evaluate every colour sum exactly. The largest contracted
length is nine. The absolute raw integer bound
`9! * 2^9 * 4^9 = 48704929136640` is below the signed 64-bit limit.
All trace-polynomial coefficients and matrix linear algebra use exact
arbitrary-precision rational or integer arithmetic.

`full_gram_Y*.json.gz` contains the complete universal centre plus all
incoming images, including all their mixed pairings. `finite_metric_Y*.json.gz`
contains every finite-rank relation, the nonsingular finite Gram matrix and
the exact expression of every universal row in that finite basis.
The equality reconstructing the full Gram matrix from these expressions
is checked entry by entry. A nonsingular modular minor proves independence
over the rationals; the exact reconstruction proves the upper rank bound.
The independent character count proves completeness.

Positivity of the complete cochain form follows from its explicit compact
Hermitian colour basis and the positive radial oscillator form, together
with the certified quotient by every null vector. On the harmonic space,
an exact fraction-free factorization additionally supplies a positive
diagonal and verifies the entire factorization identity.

For each finite basis vector, the full outgoing image is evaluated.
`outgoing_gram_Y*.json.gz` is the physical Gram matrix of all those images.
The equations solved simultaneously are:

1. Zero pairing with every outgoing image under that positive form. This is
   equivalent to the entire physical outgoing vector vanishing.
2. Zero pairing with every incoming exact vector under the centre's radial
   form. This is precisely the radial-adjoint condition, since every allowed
   predecessor is present.

No state is obtained by applying a projection. The saved exact solution
columns obey both entire matrices of equations. The independently computed
incoming ranks, outgoing ranks and harmonic dimensions add to the complete
finite dimensions in every row of the table.

At this order, radial conjugation is the free radial conjugation and the
one-loop Hamiltonian coefficient is

\[
 [D-D_{\rm classical}]_{\text{order }g_{\rm YM}^2}
 =\frac{g_{\rm YM}^2}{8\pi^2}
   \{Q_{\rm native},Q_{\rm native}^{\dagger_0}\}.
\]

Every exported state is annihilated by both operators in the anticommutator.
The radial-adjoint condition is established from actual Wick pairings,
rather than from a coefficient dual.

## Fortuity of the entire block

The covering is the universal graded symplectic trace algebra, with cyclic
and symplectic reversal identities and no finite-rank trace relations.
For every length, the entire finite relation space is differentiated in
that covering. The exact obstruction quotient dimensions are

\[
\begin{array}{c|rrrr}
\text{length}&5&6&7&8\\\hline
\operatorname{rank}Q_{\rm universal}&17&199&365&103\\
\operatorname{rank}Q_{\rm universal}(\text{all finite relations})&0&43&134&58\\
\operatorname{rank}Q_{\rm finite}&17&156&230&45\\\hline
\text{fortuitous quotient dimension}&0&0&1&0
\end{array}
\]

Subtracting the second and third lines from the first gives the last line. This
counts finite cycles which admit no closed universal representative.
All incoming boundaries descend from universal boundaries and therefore
have zero obstruction. In length seven the entire finite cohomology is
one-dimensional, so the exported harmonic state is necessarily fortuitous.
In lengths six and eight every cohomology class lies in the stable image.
No Cartan restriction or index subtraction replaces this calculation.

In the primitive integral trace normalization of
`physical_representatives.json.gz`, the length-seven state has 8489 trace
terms and positive radial norm

\[
1178509826288953100895665738079738765312000.
\]

The single length-eight state has 2040 trace terms and norm

\[
1210668214956478925280385267177499383024406541829920363133440000.
\]

All 21 length-six states and their 441 pairings are retained in the full
matrix. The mixed native-length pairings vanish in this free form because
Wick pairings require the same number of letters. This zero is derived at
the stated order; it is not declared to be an exact physical selection rule.
An exact invertible basis change mixes all 23 multiplicity directions,
including the fortuitous direction and both neighbouring length blocks.
The transformed entire Gram matrix, its unchanged determinant, and an
unchanged quadratic norm are checked and exported.

## Files and replay

| File | Complete content |
|---|---|
| `declaration.json` | Theory, charges, full symmetry representation and coupling scope, selected before the metric |
| `frames_Y*_offset*.json.gz` | All universal centre/predecessor bases and exact supercharge images |
| `full_gram_Y*.json.gz` | All centre and incoming-image radial pairings |
| `finite_metric_Y*.json.gz` | All finite nulls, independent basis, complete finite Gram and exact reconstruction |
| `outgoing_gram_Y*.json.gz` | Gram of every outgoing image |
| `harmonic_Y*.json.gz` | Exact simultaneous physical supercharge/radial-adjoint solutions |
| `physical_representatives.json.gz` | All 23 multiplicity representatives and all 69 flavour weight states, as explicit trace polynomials |
| `physical_metric.json.gz` | Entire 23 by 23 and 69 by 69 forms, with an explicit mixed basis change |
| `character_receipt.json` | Independent exact finite-character completeness |
| `fortuity_receipt.json` | All four complete universal-to-finite obstruction calculations |
| `independent_checks.json` | Independent colour expansion, radial spinor normalization and descendant recurrences |
| `block_receipt.json` | Complete flavour identities, positivity and mixed-basis covariance |
| `input_manifest.json` | Public base pin and hashes of every imported repository source |
| `execution_receipt.json` | Executed replay steps and result |

The physical trace-polynomial encoding is a list of pairs
`[integer coefficient, list of traces]`; each trace contains native letters
`[a,b,list of flavour indices]`. Flavour indices 0,1,2 encode the source
paper's indices 1,2,3. All factorials and physical normalizations are fixed
above. There are no omitted polynomial terms.

From the repository root, with Python 3.12, NumPy, SymPy, python-flint 0.9.0,
and a C++17 compiler:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algorithms/physical_block_20261124/reproduce.py --full
```

This regenerates every native frame and every Wick matrix before checking
the exact certificate. `--verify` retains the already calculated Wick
matrices and repeats the complete rational reductions, independent checks,
character counts, fortuity comparison and explicit state export. It clearly
records this narrower replay mode. Both modes use one worker.

The execution receipt is marked `running` before any regeneration and after
each completed step. A stopped process therefore cannot leave an earlier
success receipt as the status of the current packet. After a successful
reproduction, freeze or check all code, inputs, outputs and the log with

```bash
python3 algorithms/physical_block_20261124/freeze_manifest.py --write
python3 algorithms/physical_block_20261124/freeze_manifest.py
```

`MANIFEST.json` records the exact executed command and every packet file's
SHA-256. Reproduction changes measured timings, so the manifest must be
explicitly refreshed only after a successful run. The compiled Wick library
is rebuilt from the pinned C++ source and is not a portable input artifact.

## Remaining physical matching

F03 is not supplied by this packet. The leading physical supercharge and
the free radial form fix the one-loop zero-energy states above. They do not
fix the quantum correction to the action on renormalized composite local
operators, the matched current contacts or local products, or the next
metric correction. A formal higher-order survivor statement cannot replace
those physical data. The exported complete bases, all finite relations and
both incoming/outgoing interfaces provide the required fixed starting data
for that separate calculation.
