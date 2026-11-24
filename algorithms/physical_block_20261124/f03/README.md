# Complete quadratic nonscalar channel: physical comparison at first quantum order

## Scope and review status

This packet is a **candidate channel-specific F03 comparison**, supplied for
independent review. No roadmap status or score is changed. It concerns
undeformed four-dimensional N=4 super-Yang--Mills with gauge group USp(4),
the complete source charges

\[
E=\frac52,\qquad (R_1,R_2,R_3)=\left(\frac12,\frac12,\frac12\right),
\qquad (J_L,J_R)=\left(\frac12,0\right),
\]

and the complete target charges

\[
E=3,\qquad (R_1,R_2,R_3)=(1,1,1),\qquad (J_L,J_R)=(0,0).
\]

The calculated coefficients are the coefficient of `g_YM` and the coefficient
of `g_YM^2` in the action of the physical supercharge on canonically rescaled
normal products. The second is its first nonzero quantum coefficient in this
channel. A spectral-sequence order is not identified with a dilatation loop
order. This is not the orthogonal third-order benchmark, a fortuitous lift,
or a comparison on the separate 69-state F02 block.

The exact algebraic calculations in `channel.py` are complete. The physical
argument below uses the standard stress-tensor and Konishi multiplets and a
specified supersymmetric minimal-subtraction scheme; the review must assess
the local-current/product argument as well as the numerical receipt.

## 1. Complete native spaces, including every derivative alternative

Use exactly the physical letters, spinor conventions and canonical Taylor
rescaling of the parent F02 packet. A word of length `L`, total jet degrees
`a,b` and total flavour degrees `n1,n2,n3` has

\[
2E=2a+2b+n_1+n_2+n_3+L,
\qquad 2J_L=a+b+n_1+n_2+n_3-L,
\qquad 2J_R=a-b,
\]

\[
2R_i=L-n_1-n_2-n_3+2n_i.
\]

These are additive letter charges. The exact nonnegative solutions, before
discarding zero traces, are:

| Space | Native weight `(a,b,n1,n2,n3)` | Length | Number of independent traces |
|---|---|---:|---:|
| Source | (1,1,0,0,0) | 1 | 0 |
| Source | (0,0,1,1,1) | 2 | 3 |
| Target | (1,1,0,0,0) | 2 | 1 |
| Target | (0,0,1,1,1) | 3 | 1 |

The one-letter source is a covariant derivative of the gaugino. Its gauge
trace vanishes identically. Thus it is explicitly accounted for, not removed
by restricting to derivative-free expressions. No physical spacetime
derivative is omitted. The target gaugino pair is retained in addition to
the cubic scalar word. Its native `a,b` degrees encode the two primary
dotted-gaugino components.

The three physical sources are the three individual operators
`Tr(phi^i psi_i)`, with no summation in an individual source. The complete
source contains both independent differences and their sum. The cubic
target obeys symplectic reversal

\[
\operatorname{Tr}(\phi^1\phi^2\phi^3)
 +\operatorname{Tr}(\phi^1\phi^3\phi^2)=0.
\]

All other possible traces contain a one-letter trace and vanish. The second
target is the covariant spinor contraction
`Tr(lambda^{dot alpha} lambda_{dot alpha})`. There are no additional finite
rank relations: the complete native source Gram matrix is `10 I_3`, and the
complete target Gram matrix is `diag(15,10)` in the polynomial basis exported
by the script. The entire 5 by 5 pairing is independently recomputed using
all ten compact Hermitian colour generators. Its nonzero diagonal is

\[
10,\quad10,\quad10,\quad15,\quad10.
\]

The physical free bound is additive and nonnegative on the elementary
letters. An operator with these saturated charges can therefore only use
the retained letters. This excludes additional non-BPS letters as possible
local counterterms or descendants at the same classical dimension.

## 2. Renormalization and the operator dictionary

Use supersymmetric dimensional reduction in `4-2 epsilon` dimensions,
minimal subtraction of poles, and no additional finite change between the
bare cubic scalar and gaugino bilinear. This is the normal-product scheme
used for the first mixing orders in [Eden, 0911.4884v2, sections 1--2 and
appendix](https://arxiv.org/pdf/0911.4884v2). That paper's symbols and field
phases are not silently substituted for the CFLT convention. In particular,
the relative phase is derived directly in section 3 below.

At the required order the complete freedom is as follows.

| Object | Complete possible mixing | Scheme or physical condition used |
|---|---|---|
| Three source operators | A singlet factor and one common factor on the two independent traceless combinations | Minimal subtraction; traceless flavour-current normalization is fixed by its charge action |
| Cubic and gaugino targets | A general 2 by 2 matrix with even powers on its diagonal and odd powers between the two lengths | Pole subtraction in the bare normal-product basis; the finite order-`g_YM` coefficients of the physical multiplet eigenoperators are calculated below |
| Physical supersymmetry current | An overall charge normalization and identically conserved improvements | Charge normalization fixed by CFLT (D.3); improvements have zero charge and cannot change the channel's residue |
| Current--source contact counterterms | Local operators of dimension at most two with target R-charges | There are none, by the complete charge count |
| Derivatives of delta-function contact terms in its divergence | Operators of dimension at most two with target R-charges | There are none, by the same count |

The coupling/field sign symmetry `g_YM -> -g_YM`, every canonical field
`-> -field`, is preserved by this scheme. It forces the stated parity in
the target mixing matrix. In particular there is no order-`g_YM` source
renormalization, and no coefficient of `g_YM^2` from a quadratic source to
the cubic target.

There is no diagram for the cubic--gaugino separated two-point function at
order `g_YM`: a single Yukawa vertex that could join two identical
R-symmetry-index gauginos would require `Phi_44`, which is zero. A gauge
vertex or scalar vertex has the wrong field content. Hence both the pole
and finite mixed two-point coefficients at this order vanish. This is the
explicit field-content reason for the `BF1=0` entry in [Eden, equation
(17)](https://arxiv.org/pdf/0911.4884v2). It is group independent; no large
rank assumption is used. Thus there is no order-`g_YM` divergent conversion
of either target normal product. A finite conversion is excluded by the
stated minimal-subtraction definition, not presumed to vanish in every
scheme.

The exact leading native differential and its free radial adjoint are

\[
\begin{pmatrix}2&2&2\\0&0&0\end{pmatrix},
\qquad
\begin{pmatrix}3&0\\3&0\\3&0\end{pmatrix}.
\]

Their product on the source has eigenvalues `0,0,18`. With the parent
packet's physical normalization, the nonzero one-loop anomalous dimension
is therefore

\[
\frac{g_{\rm YM}^2}{8\pi^2}\,18
 =\frac{9g_{\rm YM}^2}{4\pi^2}.
\]

The two zero directions are descendants of the conserved traceless flavour
currents; the sum is a descendant of the Konishi scalar. Since supersymmetry
is manifest in the chosen regulator, a primary and its descendant have the
same pole renormalization. The complete coefficient of `g_YM^2` in the
source renormalization matrix is consequently

\[
\frac{3}{8\pi^2\epsilon}
\begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix}.
\]

Its singlet eigenvalue is `9/(8*pi^2*epsilon)` and its two other eigenvalues
are zero. The pole coefficient follows from the displayed anomalous
dimension and the minimal-subtraction relation in `4-2 epsilon`
dimensions; it agrees with the explicit calculation of [Eden, equations
(20)--(23)](https://arxiv.org/pdf/0911.4884v2), with the colour factor fixed
here by exact USp(4) contractions. The complete target pole coefficient at
this order is `diag(9/(8*pi^2*epsilon),0)`. The protected target has no pole;
the cubic target has the same Konishi pole. No other operator exists in
either charge sector to carry an omitted mixing.

These order-`g_YM^2` pole factors first enter the nonzero supercharge action
at order `g_YM^3`, because its free action on the source is zero. They are
recorded to specify the renormalized operator dictionary, but are not
mistaken for the order-`g_YM^2` anomaly.

## 3. Independent physical determination of the anomaly, including its phase

The following uses the actual N=4 field transformations [CFLT,
(3.1),(D.3), and spinor convention (A.4)](https://arxiv.org/pdf/2306.04673v2).
Keep the dotted spinor indices contracted covariantly. They give

\[
[Q_{\dot\alpha 1},\phi^1]=\lambda_{\dot\alpha},
\qquad
\{Q_{\dot\alpha 1},\lambda_{\dot\beta}\}
 =-2\epsilon_{\dot\alpha\dot\beta}[\phi^2,\phi^3].
\]

The second equation follows by summing both `p=2,3` terms in
`[Phi_1p,Phi^{4p}]`, using the scalar reality convention. Consequently

\[
Q_{\dot\alpha 1}Q_{\dot\beta 1}\operatorname{Tr}((\phi^1)^2)
 =-\epsilon_{\dot\alpha\dot\beta}
  \operatorname{Tr}\left(
   \lambda^{\dot\gamma}\lambda_{\dot\gamma}
    +4\phi^1[\phi^2,\phi^3]\right).
\]

This is the descendant of a protected half-BPS stress-tensor-multiplet
primary. Its displayed relative coefficient is fixed by the elementary
supersymmetry transformations. Loop corrections cannot modify this first
relative coefficient: it is the order-`g_YM` cubic term after canonical
field normalization. This is the CFLT-convention counterpart of the
protected descendant of [Eden, equations (4)--(7)](https://arxiv.org/pdf/0911.4884v2).

The source sum is the actual first supersymmetry descendant of the Konishi
scalar, as is verified directly in [CFLT, (4.13)](https://arxiv.org/pdf/2306.04673v2).
It has the nonzero anomalous dimension just calculated. Its image under
the chosen physical supercharge therefore belongs to the long Konishi
multiplet. Superconformal covariance of separated two-point functions
requires zero pairing between this image and the protected descendant.
This physical restriction is also the starting point of [Eden--Jarczak--
Sokatchev--Stanev, 0501077v2, sections 1--2](https://arxiv.org/pdf/hep-th/0501077v2),
but all colour norms and the complete channel used here are calculated
independently.

In the exported native target order, cubic first and gaugino second, the
protected direction through its first two coefficients is

\[
\begin{pmatrix} i g_{\rm YM}/(2\pi)\\ 1\end{pmatrix}.
\]

The complete metric to the needed order is `diag(15,10)`; its order-`g_YM`
coefficient is zero by the diagram argument in section 2. With the classical
coefficient fixed, physical orthogonality now gives the entire long
direction:

\[
\begin{pmatrix}-i g_{\rm YM}/(2\pi)&1\end{pmatrix}
\begin{pmatrix}15&0\\0&10\end{pmatrix}
\begin{pmatrix}1\\3i g_{\rm YM}/(4\pi)\end{pmatrix}
 =-\frac{15i g_{\rm YM}}{2\pi}
   +\frac{30i g_{\rm YM}}{4\pi}=0.
\]

There is no omitted target direction in this equation. Flavour covariance
makes all three source columns identical; both independent source
differences are included and map to zero. In the canonical native source
and target bases, the full coefficient matrices are

\[
[g_{\rm YM}]\,Q_{\rm physical}
 =\begin{pmatrix}
 1/(2\pi)&1/(2\pi)&1/(2\pi)\\0&0&0
 \end{pmatrix},
\]

\[
[g_{\rm YM}^{2}]\,Q_{\rm physical}
 =\begin{pmatrix}
 0&0&0\\3i/(8\pi^2)&3i/(8\pi^2)&3i/(8\pi^2)
 \end{pmatrix}.
\]

Here `[g_YM^n]` denotes the coefficient at exactly that power, not a claim
that all higher coefficients vanish. The full radial-adjoint coefficients
are obtained with the same complete physical metrics and complex
conjugation. In particular its action on the displayed protected target
vanishes at the two required orders. `channel.py` checks the full matrices,
both independent kernel columns and this adjoint identity.

In CFLT's unrescaled physical fields this yields, for an individual source
and with no summation over `i`, the first two field-degree-compatible terms

\[
\{Q^4_\alpha,\operatorname{Tr}(\phi^i\psi_{\beta i})\}
 =\epsilon_{\alpha\beta}
 \left[
  2i\operatorname{Tr}(\phi^1[\phi^2,\phi^3])
  -\frac{3i g_{\rm YM}^2}{16\pi^2}
       \operatorname{Tr}(\lambda^{\dot\gamma}\lambda_{\dot\gamma})
 \right]
\]

at precisely the canonical coefficient orders specified above. The index
choice `alpha=-, beta=+`, with `epsilon_{-+}=-1`, gives the positive
`+3i*g_YM^2/(16*pi^2)` anomaly in the fixed source convention. This phase
is derived from CFLT's transformations and radial conjugation. A direct
substitution of differently phased fields from another paper can reverse
it incorrectly.

## 4. Independent point-split residue and the native formal differential

For the source normal product, use both Wilson factors in the adjoint
trace: transport the fermion from one end of the split to the scalar's end
and transport back around the trace. Varying the Wilson factors produces
the commutator insertion. The relevant short-distance graph contains the
Yukawa vertex with the opposite scalar and matter-fermion fields; all
uncontracted fields are the two dotted gauginos.

The exact finite convolution, before the coincident limit, is

\[
\int_{\mathbb R^4}
 \frac{(y+\varepsilon/2)_\mu\,d^4y}
 {|y-\varepsilon/2|^2|y+\varepsilon/2|^4}
 =\frac{\pi^2\varepsilon_\mu}{|\varepsilon|^2}.
\]

Feynman parametrization gives the coefficient directly:

\[
2\int_0^1 t(1-t)\,
 \frac{\pi^2}{2t(1-t)|\varepsilon|^2}\,dt
 =\frac{\pi^2}{|\varepsilon|^2}.
\]

The two uncontracted identical gauginos are antisymmetric in their dotted
indices. Their spinor contraction eliminates the direction dependence of
`epsilon_mu epsilon_nu/epsilon_squared` and supplies the factor `1/2`.
No angular averaging is needed. The script verifies every coefficient of
the corresponding quadratic separation polynomial in the exact CFLT sigma
matrix convention.

The Wilson variation contributes magnitude 2, the two terms of the CFLT
Yukawa vertex contribute 128/`g_YM^2`, each contracted propagator contributes
`g_YM^2/(64*pi^2)`, and the compact colour contraction is 6. Therefore the
entire residue magnitude is

\[
2\,\frac{128}{g_{\rm YM}^{2}}
 \left(\frac{g_{\rm YM}^{2}}{64\pi^{2}}\right)^2
 \pi^2\,\frac12\,6
 =\frac{3g_{\rm YM}^{2}}{16\pi^{2}}.
\]

Its relative phase is fixed independently by the physical multiplet
calculation in section 3. This check is not a replacement of a physical
bra by a coefficient dual. It independently verifies the magnitude of the
same anomalous local product, using the action and an exact finite
short-distance integral.

The complete compact colour calculation checks all 100 entries of

\[
\sum_A\operatorname{Tr}([T_A,T_B][T_A,T_C])=-6\delta_{BC},
\qquad \operatorname{Tr}(T_A T_B)=\delta_{AB}.
\]

Applying the pinned native pair rule and this exact contraction to every
source gives

\[
\begin{pmatrix}0&0&0\\-12&-12&-12\end{pmatrix}.
\]

Thus the native rule used by the repository, whose overall constant was
deliberately set to one, matches the physical first quantum coefficient
with multiplier `-i*g_YM^2/(32*pi^2)`. The classical rule has multiplier
`g_YM/(4*pi)`. These are explicit channel normalizations; a formal loop
parameter is not silently identified with a physical coupling power.

The version matters: [Budzik--Kulp, 2512.07771v3, equations (10),(70)--(71)
and section 3.4](https://arxiv.org/html/2512.07771v3) has corrected numerical
factors relative to v2. The repository's pair implementation removes the
common coefficient intentionally. This calculation fixes that coefficient
from physical operators rather than importing either version's convention
without a translation.

## 5. Local products, current contacts and the remaining review question

The local product dictionary is the minimal-subtracted gauge-covariant
normal product specified in section 2. Every source and target product is
exported as an entire trace polynomial. The finite order-`g_YM` basis change
to the protected and long physical target operators is the two-column
matrix

\[
\begin{pmatrix}
1&i g_{\rm YM}/(2\pi)\\
3i g_{\rm YM}/(4\pi)&1
\end{pmatrix}.
\]

The source's full order-`g_YM^2` pole matrix and the target's corresponding
pole matrix were given explicitly above. Finite diagonal renormalizations
are fixed to zero in the normal-product scheme. An additional finite
off-diagonal normal-product redefinition would change the displayed
quantum differential by the corresponding change of basis, and is not
being silently allowed. The first relative mixing of the physical target
multiplets is not a removable omitted operator: both target directions are
retained.

The current under comparison is the physical supersymmetry Noether current
of the same N=4 action, with its charge normalized to CFLT (D.3), regulated
with manifest N=1 supersymmetry. This is an ordinary supersymmetry Ward
identity for composite operator insertions, not the restriction of a
Konishi rescaling Ward identity to chiral insertions. The Konishi anomaly
modifies the composite operator's transformation; it is not an anomaly in
that fundamental supersymmetry Ward identity. This distinction and the
renormalization of its composite descendant are explicit in [Eden,
section 1, especially (11)--(14)](https://arxiv.org/pdf/0911.4884v2).

For the current--source product, the current has dimension `7/2` and the
source dimension `5/2`. A local contact counterterm in this product would
require a dimension-two operator with target R-charges `(1,1,1)`. No such
operator exists. A derivative of a delta function needs still lower
dimension. After taking the divergence, the only possible contact is a
delta function times a dimension-three target, and its complete target
space is exactly the two-dimensional space already calculated. The
current's identically conserved improvements cannot change its sphere
charge; an improvement-induced local derivative contact is excluded by
the same dimension and charge argument.

The resulting delta-contact residue is the physical 2 by 3 supercharge
matrix in section 3. It is determined by an independently fixed physical
source (the Konishi/flavour descendants), the elementary protected-target
identity, and the full positive physical target pairing. It is not fixed
merely by naming a Ward identity or by counting the target dimension.
Nonchiral physical bras are included: the full target space and its radial
adjoint, not just chiral vacuum expectation values, determine the residue.
There is no room for an additional local residue with zero pairing against
that complete positive target form.

**Question reserved for independent review.** The packet has not expanded
the entire renormalized Noether-current density in component fields or
recomputed all current-insertion Feynman graphs in supersymmetric
dimensional reduction. The argument above instead uses the manifest
supersymmetry Ward identity, the explicit normal-product scheme, absence
of admissible contact counterterms, and the independent physical
multiplet calculation. If the F03 test requires a direct graph-by-graph
current-density comparison beyond this channel-level determination, that
calculation remains unprovided and F03 must remain partial. No all-order
untwisted comparison theorem is assumed.

## Replay

From the repository root:

```bash
PYTHONPATH=/workspace/scratch/24bd17e8b09b/research_deps OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algorithms/physical_block_20261124/f03/channel.py
```

`channel_receipt.json` contains every source/target trace polynomial, both
complete matrices, the two-dimensional source kernel, the physical
adjoint, all 25 independent Wick pairings, all 100 colour contractions,
and the exact point-split integral and spinor checks. The separate F02
manifest is unchanged.
