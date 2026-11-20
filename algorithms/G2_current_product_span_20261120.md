# Two $G_2$ singlet source sectors

The exact stress-tensor multiplet supplies the currents $\mathbf L_{Ii}$, $\mathbf G^I_j$ and the gaugino bilinear $S$ appearing in (4.83)–(4.85). The primary cohomological product of two such currents is closed by the supersymmetry Ward identity. At separated points the supercharge acts on both factors. In two complex dimensions the degree-zero holomorphic OPE extends across coincidence. Its classical term is the ordinary product. Holomorphic derivatives commute with the supercharge. These are the physical inputs in applying the following rank test to an instanton differential. [Primary product, §2; physical currents, §4.3.1](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

At doubled charges $(1,1,4,4,4)$ and $(1,1,6,6,6)$, the two derivative degrees are equal. Their common value determines the length and the three flavour degrees. All possibilities are included below. The incoming and outgoing ranks are lower bounds over the rationals from integer minors at the primes $1{,}000{,}003$ and $1{,}000{,}033$.

| Doubled flavour charges | Length | Cochains | Incoming rank | Outgoing rank | Cohomology upper bound |
|---|---:|---:|---:|---:|---:|
| $(4,4,4)$ | 5 | 3 | 0 | 3 | 0 |
| | 4 | 5 | 1 | 3 | 1 |
| | 3 | 1 | 0 | 1 | 0 |
| $(6,6,6)$ | 8 | 15 | 6 | 9 | 0 |
| | 7 | 48 | 28 | 20 | 0 |
| | 6 | 30 | 17 | 11 | 2 |
| | 5 | 3 | 0 | 3 | 0 |

At length four the singlet $\epsilon^{ij}\mathbf L_{Ii}\mathbf G^I_j/6$ gives $4\eta$ on the matter-free Cartan jet, where $\eta$ is its one exterior generator. All classical boundaries vanish there. The archived source gives the same nonzero value. A one-dimensional cohomology leaves just one line for the leading term of this closed product. Finite operator mixing within the line preserves this statement. [Four-letter calculation](verify_g2_current_product.py)

At length six use the two physical products in the order

$$
\frac{S}{6}\epsilon^{ij}\mathbf L_{Ii}\mathbf G^I_j,
\qquad
\epsilon^{IJK}\epsilon^{ik}\epsilon^{jh}
\mathbf L_{Ii}\mathbf L_{Jj}\partial_k\mathbf L_{Kh}.
$$

Here $i,j,k,h$ label the two spin components and $I,J,K$ the three flavours, with the alternating tensors of the current formulas. In the physical component convention the two Cartan functionals, each zero on classical boundaries, give

$$
\det\begin{pmatrix}-64&0\\16&576\end{pmatrix}=-36{,}864.
$$

The first row contains one matter fermion and its only nonzero entry comes from the matter term of (4.84). Reversing the matter-component orientation reverses this row and preserves the rank. The second row contains a scalar derivative and survives with the matter field set to zero. Thus both physical products are independent in classical cohomology. They exhaust its upper bound of two. The previous component comparison left a relative sign open for a literal equality of representatives. This span calculation needs only the rank. [Component comparison](G2_protected_incoming_currents_20261118.md) · [Replay with both primes](verify_g2_current_span.py)

For a closed exact product, taking its leading term and then its successive corrections gives a solution of each outgoing descent equation. The dimension of a later source page cannot exceed that of its classical source. The exact products span both classical source spaces, so the outgoing instanton differential vanishes on their perturbative survivors under the primary-product input above. This excludes the proposed incoming pairing into $S^2$ or the scalar Gram determinant from the first singlet sector, and into $S^3$ from the second. The argument determines these outgoing ranks from symmetry and the product. Direct instanton matrix elements remain unevaluated. [Exact dimensions and determinants](../data/N4_G2_protected_source_span_20261120.json)
