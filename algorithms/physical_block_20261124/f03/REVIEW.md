# Independent acceptance review: complete first quantum-order channel

**Decision: F03 is met for the complete USp(4) channel specified in README.md,
through the coefficient of the second power of the physical coupling.**
This is a channel-specific comparison. It does not extend the comparison to
the 69-state block in the parent directory, to higher perturbative orders,
to instantons, or to the orthogonal third-differential benchmark.

## Exact replay and independent physical checks

The root reviewer ran `channel.py` again. All three source polynomials, both
target polynomials, 25 independently evaluated Wick entries, 100 adjoint
colour contractions, both source-kernel columns, both complete charge
matrices and their radial adjoints passed. The source and target forms are
respectively `10 I_3` and `diag(15,10)`. The script exhausts the physical
charges before discarding the identically zero one-letter trace.

The physical comparison was checked against the primary sources, rather
than accepted from the numerical receipt:

- Chang, Feng, Lin and Tao, [2306.04673v2](https://arxiv.org/html/2306.04673v2),
  equations (3.1), (3.2), (4.13), (D.3), and (D.9).
- Eden, [0911.4884v2](https://arxiv.org/pdf/0911.4884v2), equations (7),
  (11), (17), and (20)--(23). In particular the separated cubic--gaugino
  pairing has zero coefficient at the first power of the coupling.
- Eden, Jarczak, Sokatchev and Stanev,
  [hep-th/0501077v2](https://arxiv.org/abs/hep-th/0501077v2): the physical
  Konishi and protected descendants are different renormalized mixtures.

The root reviewer independently translated the three polynomial types from
the superfield, keeping the spinor contraction covariant. The exported
native source is

\[
-\frac{64\pi^2}{g_{\rm YM}^2}\operatorname{Tr}(\phi^i\psi_i),
\]

with no sum in a source column. The exported native cubic and gaugino
targets, in that order, are

\[
\frac{256i\pi^3}{g_{\rm YM}^3}
\operatorname{Tr}(\phi^1[\phi^2,\phi^3]),
\qquad
-\frac{32\pi^2}{g_{\rm YM}^2}
\operatorname{Tr}(\lambda^{\dot\alpha}\lambda_{\dot\alpha}).
\]

The protected descendant therefore has native coordinates
`(i*g_YM/(2*pi),1)` through these two coefficients. Orthogonality to it,
using the complete physical form and the vanishing first-order mixed
two-point coefficient, fixes the Konishi direction to
`(1,3*i*g_YM/(4*pi))`. With the independent classical normalization, every
source column of the physical supercharge is

\[
\begin{pmatrix}g_{\rm YM}/(2\pi)\\
3i g_{\rm YM}^2/(8\pi^2)\end{pmatrix}
\]

at the displayed powers only. The independent Wilson-line convolution
reproduces the magnitude. The phase follows from the protected descendant
and the explicit superfield dictionary, and agrees with (4.13).

## Why the current comparison is determined in this channel

The F03 gate permits a complete channel-specific derivation; it does not
require a graph-by-graph current expansion when an independent Ward
argument fixes the same operation. The current is the supersymmetry current
of the specified physical action, normalized by its elementary charge
action in a regulator preserving the chosen supersymmetry.

The current--source product has total classical dimension six. A local
delta-function counterterm would require a physical gauge-invariant
operator of dimension two and target R-charges `(1,1,1)`. The complete charge
count excludes it. Derivatives of a delta function require still lower
dimensions and are excluded as well. Identically conserved improvements
do not change the sphere charge.

After taking the divergence, the delta residue has dimension three. Its
entire physical target is the two-dimensional space already enumerated.
The independently fixed source descendants, protected-target identity and
positive target form determine all three charge-action columns. An
additional residue cannot be invisible to all physical bras in this
complete nondegenerate space. This uses the supersymmetry Ward identity
for general physical insertions, rather than a chiral-only Konishi identity.

Finite target changes are scheme dependent. The minimal-subtracted
normal-product basis and the first physical eigenoperator mixing are
specified separately in README.md. The source and target pole matrices
are given there; their effects on the nonzero supercharge action begin at
the third coupling power. They cannot alter the two coefficients accepted
here. No higher-order exactness of the anomaly is inferred.

## Scope retained after acceptance

The full component Noether-current density has not been expanded, and a
separate graph-by-graph dimensional-reduction current calculation has not
been performed. The accepted result is the above channel-level Ward and
physical-multiplet determination, with an independent point-split magnitude
check. A different sector requires its own complete mixing/contact analysis.
