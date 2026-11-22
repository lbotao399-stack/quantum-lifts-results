<a id="np-section-14"></a>

## 16. Four-dimensional determinant lift and microscopic incoming class

<a id="np-69504578f9778a53-1"></a>

### Completely Higgsed magnetic $SU(2)$

A Higgs patch is a region where nonzero scalar expectation values give masses to gauge bosons.

Complete Higgsing means that the gauge transformations fixing those scalar values form the identity subgroup.

Here the magnetic group is $SU(2)$, and the massive flavor supplies a nonzero quark–antiquark contraction.

The light fields below are coordinates on the region where that contraction is nonzero.

Take electric $N_f=N_c+2$ and add a mass $m\ne0$ for the last flavor.

The heavy F-equation gives $q_{N_f}\widetilde q^{\widetilde{N_f}}=-\mu m$, completely Higgsing magnetic $SU(2)$.

Hats mark the retained $N_c+1$ light flavors.

The source one-instanton calculation gives

$$
W_{\mathrm{pert}}=\frac1\mu
\widehat M^i{}_{\tilde\jmath}\widehat q_i
\widehat{\widetilde q}^{\tilde\jmath},
$$

$$
W_{\mathrm{inst}}
=\frac{\widetilde\Lambda^{4-N_c}\det(\widehat M/\mu)}
{q_{N_f}\widetilde q^{\widetilde{N_f}}}
=-\frac{\widetilde\Lambda^{4-N_c}}{\mu^{N_c+2}m}\det\widehat M
=-\frac{\det\widehat M}{\Lambda_L^{2N_c-1}},
\qquad
\Lambda_L^{2N_c-1}=m\Lambda^{2N_c-2}.
$$

This microscopic instanton superpotential is known.

The following calculation gives an explicit consequence in its effective F-term ideal.

[Seiberg, §4](https://arxiv.org/abs/hep-th/9411149); [Intriligator–Seiberg, (5.11)–(5.17)](https://arxiv.org/abs/hep-th/9509066)

**Proposition 16.1. **

**Vanishing of the determinant class.**

At fixed magnetic light coordinates,

$$
\widehat M^i{}_{\tilde\jmath}
\frac{\partial W}{\partial\widehat M^i{}_{\tilde\jmath}}
-\widehat q_i\frac{\partial W}{\partial\widehat q_i}
=-(N_c+1)\frac{\widetilde\Lambda^{4-N_c}}{\mu^{N_c+2}m}
\det\widehat M.
$$

> **Proof.**
>
> The cubic coupling has Euler weight $1$ in each term, so its two contributions cancel.
>
> The determinant has degree $N_c+1$, giving the right-hand side.
>
> In the F-term quotient of $W_{\mathrm{pert}}$, setting $\widehat q=\widehat{\widetilde q}=0$ defines an evaluation into $\mathbb C[\widehat M]$, hence $[\det\widehat M]\ne0$.
>
> In the parent theory, the heavy Higgs fields satisfy F- and D-flatness, the light quarks vanish, and $\widehat M$ is an arbitrary invertible matrix.
>
> Perturbative nonrenormalization preserves this family.
>
> Its supersymmetric vacuum expectation $\det\langle\widehat M\rangle\ne0$ excludes perturbative exactness.
>
> At a fixed nonzero instanton coefficient, the identity places the determinant in the exact effective F-term ideal.
>
> $\square$

> **Example 16.1. **
>
> **The smallest determinant instanton.**
>
> The smallest example is $SU(2)$ with $N_f=4\to3$:
>
> $$
> \Lambda^2\widetilde\Lambda^2=\mu^4,
> \qquad
> W_{\mathrm{inst}}=-\frac{\widetilde\Lambda^2}{\mu^4m}\det\widehat M,
> \qquad
> \Lambda_L^3=m\Lambda^2.
> $$
>
> The one-instanton zero-mode saturation is
>
> $$
> (4\text{ gaugino}+8\text{ quark})
> -(2\text{ gaugino}+2\text{ quark})
> -6\text{ quark}=2\text{ supersymmetry modes}.
> $$
>
> Three light-quark mass insertions give $\det(\widehat M/\mu)$.
>
> Electric-normalized baryons satisfy $B_i=\sqrt{\Lambda_L^3/\mu}\,\widehat q_i$.
>
> This rescaling is singular at vanishing magnetic instanton weight.
>
> The perturbative comparison therefore holds $\widehat q,\widehat{\widetilde q},\mu,m$ fixed.

In the retained antichiral convention, $\gamma$ are even letters, $\beta$ are their odd dual letters, and $Q\beta=-\partial\bar W/\partial\gamma$.

Here an even letter is a commuting scalar coordinate and its odd dual letter anticommutes.

The relation $Q\beta=-\partial\bar W/\partial\gamma$ records an F-equation as a $Q$-variation.

Therefore a polynomial combination of F-equations gives an explicit odd operator whose $Q$-variation is the corresponding scalar relation.

A jet is a field together with a specified number of spacetime derivatives; the zero-jet part contains the undifferentiated fields.

Barred flavor slots belong to the conjugate flavor spaces.

The identity gives the effective primitive

<a id="np-a3c60c764266802e-1"></a>

$$
Q\left[
\frac{\bar\mu^4\bar m}{3\overline{\widetilde\Lambda}^{\,2}}
\left(
\overline{\widehat M}^{\bar i}{}_{\bar{\tilde\jmath}}
\beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}
-\overline{\widehat q}_{\bar i}\beta^{q\,\bar i}
\right)\right]
=\det\overline{\widehat M}.
$$

The primitive is polynomial in dynamical fields; its coefficient contains the inverse instanton weight.

The coefficient ring specifies which parameter-dependent coefficients an operator may contain.

A formal series in $\overline{\widetilde\Lambda}^{\,2}$ has nonnegative powers; its inverse is absent.

At fixed nonzero $\mu,m$, the same effective identity yields parameter torsion.

> **Proof: the formal determinant relation.**
>
> Multiplying the explicit primitive identity by the instanton parameter gives
>
> $$
> Q\left[\frac{\bar\mu^4\bar m}{3}\left(\overline{\widehat M}^{\bar i}{}_{\bar{\tilde\jmath}}\beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}-\overline{\widehat q}_{\bar i}\beta^{q\,\bar i}\right)\right]=\overline{\widetilde\Lambda}^{\,2}\det\overline{\widehat M}.
> $$
>
> If the determinant itself had a primitive with nonnegative formal powers, reducing that identity at zero instanton parameter would give a perturbative primitive.
>
> The nonzero perturbative determinant class excludes it.
>
> Hence the effective class is nonzero over the formal parameter ring and is annihilated by $\overline{\widetilde\Lambda}^{\,2}$.
>
> In a model with a fixed nonzero instanton coefficient, division by that coefficient gives the displayed primitive.
>
> $\square$

The formal family with nonnegative powers therefore has instanton-parameter torsion.

In the fiber at a fixed nonzero parameter, the determinant is exact.

<a id="determinant-full-torsion"></a>

### Every scalar class removed by the determinant interaction

**Proposition 16.2.** Fix nonzero $\mu,m$ and keep the magnetic light coordinates of Proposition 16.1. Take the one-instanton weight $\widetilde\Lambda^{4-N_c}$ as an independent formal parameter, including at a dimensionless threshold. In the scalar F-term ring the kernel of inverting this parameter is exactly

$$
(\det\widehat M)=
(\det\widehat M)\,\mathbb C[\widehat M].
$$

Multiplication by $\det\widehat M$ identifies this ideal, as a module over the meson polynomials, with one free copy of $\mathbb C[\widehat M]$. Its annihilator in the scalar F-term ring is the ideal generated by $\widetilde\Lambda^{4-N_c}$ and by all entries of $\widehat q,\widehat{\widetilde q}$. Every scalar loss occurs on the first instanton page, meaning the first induced differential on the zero-weight cohomology.

> **Proof.** Write the magnetic quark as a row, the antiquark as a column, and let $n=N_c+1$ be the size of the light meson matrix. The scalar F-equations are
>
> $$
> \widehat q\,\widehat M=0,\qquad
> \widehat M\,\widehat{\widetilde q}=0,\qquad
> \widehat q_i\widehat{\widetilde q}^{\tilde\jmath}
> =\frac{\widetilde\Lambda^{4-N_c}}{\mu^{N_c+1}m}
> \frac{\partial\det\widehat M}{\partial\widehat M^i{}_{\tilde\jmath}} .
> $$
>
> The derivative is the signed complementary minor. Multiplying the first two equations by the adjugate gives $\widehat q\det\widehat M=0$ and $(\det\widehat M)\widehat{\widetilde q}=0$. Proposition 16.1 gives $\widetilde\Lambda^{4-N_c}\det\widehat M=0$. At zero weight and zero quarks, evaluation sends every meson polynomial to itself. Thus no additional coefficient annihilates the determinant.
>
> Work first over polynomials in the weight. The difference between the number of quarks and the number of antiquarks grades the quotient. Replace a quark–antiquark pair by its signed minor using the last F-equation. Each replacement removes two quark letters. A homogeneous polynomial therefore has a representative with quarks alone, antiquarks alone, or mesons alone, in each case allowing meson and parameter coefficients.
>
> For positive quark degree consider the ideal generated by the $n$ entries of $\widehat q\widehat M$ in the polynomial ring on $\widehat q,\widehat M$. This ideal is radical. On $\widehat q\ne0$, each of the $n$ columns of $\widehat M$ satisfies one independent linear equation. This is an irreducible smooth locus of dimension $n+n(n-1)=n^2$. The locus $\widehat q=0$ also has dimension $n^2$ and is smooth in the equation scheme at invertible $\widehat M$, where the derivative with respect to $\widehat q$ has rank $n$. These are the two possible irreducible components. The $n$ equations have height $n$ in a polynomial ring, hence form a complete intersection. The complete-intersection unmixedness theorem and generic reducedness give radicality.
>
> At each point of $\widehat q\ne0$, $\widehat q\widehat M=0$ implies $\det\widehat M=0$. At rank $n-1$ every column of the cofactor matrix belongs to the one-dimensional left nullspace, written as a column. The last F-equation determines $\widehat{\widetilde q}$ from any nonzero entry of $\widehat q$, and the adjugate identity gives $\widehat M\widehat{\widetilde q}=0$. At rank at most $n-2$ take $\widehat{\widetilde q}=0$. Thus every such point extends to a solution at every nonzero weight.
>
> A positive-quark-degree representative that vanishes after inverting the weight vanishes at all these points. Varying the nonzero weight makes each polynomial coefficient vanish there. Its positive degree makes it vanish also on $\widehat q=0$. Radicality and the Nullstellensatz put each coefficient in $(\widehat q\widehat M)$, so the representative was zero already. Transposition proves the same result for positive antiquark degree.
>
> In degree zero the representative depends on the meson and parameter alone. Every singular meson admits a nonzero left null vector and hence occurs at nonzero weight. Its defining prime ideal is generated by the irreducible determinant. Varying the weight again shows each coefficient is a multiple of $\det\widehat M$. This proves the asserted kernel at every polynomial degree. Flatness of parameter-adic completion for a Noetherian ring preserves the kernels of multiplication by each fixed power. Completion acts trivially on the ideal killed by the first power.
>
> For the page statement, a degree-one Koszul cochain is a polynomial combination of the odd partners of the scalar fields. If it is closed at zero weight, its full variation is the first power of the weight times a scalar, since the superpotential is linear in that weight. The resulting scalar class is therefore in the ideal just calculated. Conversely, multiplying the Euler cochain of Proposition 16.1 by any meson polynomial realizes its determinant multiple on the first page. At any later page a scalar target has a representative whose multiple by a positive power is a full boundary. Reducing the computed torsion ideal at zero weight gives a target already removed on the first page. $\square$

Eliminating an inverse to the weight independently gives the saturated ideal obtained by adjoining $\det\widehat M$ to the three displayed sets of F-equations. The exact Gröbner checks at matrix sizes two and three give reduced lexicographic bases of sizes $9$ and $22$. The proof gives all sizes in the stated effective presentation. [Calculation](../checks/scalar_saturation.py)

<a id="np-f028df5dda216e10-1"></a>

### UV current, triangle coefficient, and polynomial completion

The microscopic calculation keeps the original magnetic quarks, mesons, gauge ghosts, and their derivatives.

Its operator $\mathbf Q_{\mathrm{pert}}$ is the full perturbative supercharge of that field algebra.

The effective calculation has already eliminated the massive gauge and matter fields.

Matching the two requires a rule for replacing microscopic local operators by effective ones, including terms supported when insertion points coincide.

A Konishi current generates a rescaling of matter fields.

The one-loop anomaly is the quantum failure of its classical conservation equation; the calculation below supplies the local correction needed to make the incoming operator closed.

In the four-flavor magnetic parent, $h=q_4\widetilde q^{\tilde4}$ denotes the heavy Higgs contraction.

Here $h$ is a scalar; the earlier homotopy operator belongs to a separate construction.

Take the light-row variation

$$
\delta q_i=q_i,
\qquad \delta M^i{}_{\tilde\jmath}=-M^i{}_{\tilde\jmath},
\qquad i=1,2,3,\quad \tilde\jmath=\tilde1,\ldots,\tilde4.
$$

Including the heavy column in this variation gives $\delta W_{\mathrm{UV}}=0$ off shell.

In the CDSW superspace convention, $S=-\operatorname{Tr}_{\mathrm{fund}}W_\alpha W^\alpha/(32\pi^2)$, and the corresponding renormalized Konishi current $J$ obeys

$$
\bar D^2J=-3S.
$$

This anomaly arises at one loop, so the naive Euler cochain requires perturbative completion.

[CDSW, §3](https://arxiv.org/abs/hep-th/0211170); [Seiberg, (2.3)](https://arxiv.org/abs/hep-th/0212225)

In the holomorphic Ward complex used here, $\gamma_i,\beta^i$ are dual letters for the matter representation $R$; $t_D$ are its generators, and $D,E$ are Lie-algebra component labels.

Direct contraction of the fundamental triangle gives

$$
Q_1(\gamma_i\beta^i)
=\frac{g^2}{64\pi^2}\operatorname{Tr}_R(t_Dt_E)
(\partial_{+\dot\alpha}\mathsf c^D)
(\partial_+{}^{\dot\alpha}\mathsf c^E).
$$

$g$, $\mathsf c$, and the dotted-spinor contraction retain the manuscript normalization; $\hbar$ is the loop-counting parameter.

> **Proof: coefficient.**
>
> The fiber derivative signs of both oriented matter cycles are $+1$; the two interaction coefficients contribute $(-ig)^2$.
>
> Combining the source symmetry factor, triangle kernel, simplex integral, and momentum conversion gives
>
> $$
> \frac1{2!}\left(-\frac1{8\pi^2}\right)
> \frac12\frac14(-2g^2)=\frac{g^2}{64\pi^2}.
> $$
>
> Keeping the color factor $\operatorname{Tr}_R(t_Dt_E)$ preserves the source representation normalization.
>
> $\square$

The source Euler cochain $A$ contracts the meson letters in all four columns of the three light rows, then subtracts the three light-quark contractions:

$$
A=\overline M^{\bar i}{}_{\bar{\tilde\jmath}}
\beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}
-\overline q_{\bar i}\beta^{q\,\bar i},
\qquad Q_0A=0.
$$

The ghost action on the antichiral scalar column uses the transpose representation:

$$
\mathsf c_\gamma=\mathsf c^At_A^{\mathsf T},
\qquad Q_0\mathsf c_\gamma=-ig\mathsf c_\gamma^2,
$$

$$
Q_1A=-\frac{3g^2}{64\pi^2}
\operatorname{Tr}_{\mathrm{fund}}
[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)].
$$

Let $v$ and $t$ be the heavy antichiral Higgs column and dual row, with $\bar h=tv$.

In a color component frame,

$$
G=\begin{pmatrix}v^1&-t_2\\v^2&t_1\end{pmatrix},
\qquad \det G=\bar h,
\qquad Q_0G=-ig\mathsf c_\gamma G.
$$

For $\bar h\ne0$, the manuscript Higgs-frame cochain is

$$
K=\operatorname{Tr}_{\mathrm{fund}}
[(\partial_{+\dot\alpha}G)G^{-1}
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)].
$$

> **Proof: polynomial cycle at every perturbative order.**
>
> Differentiation and graded cyclicity give
>
> $$
> Q_0K=-ig\operatorname{Tr}_{\mathrm{fund}}
> [(\partial_{+\dot\alpha}\mathsf c_\gamma)
> (\partial_+{}^{\dot\alpha}\mathsf c_\gamma)].
> $$
>
> With $a,b$ denoting abstract Higgs-color slots, its polynomial numerator is
>
> $$
> \bar hK=
> t_a(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b
> (\partial_{+\dot\alpha}v)^b
> -(\partial_{+\dot\alpha}t)_a
> (\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b v^b.
> $$
>
> The light-flavor betas and heavy Higgs scalars have different flavors.
>
> The triangle contractions give $Q_1(\bar hA)=\bar hQ_1A$.
>
> Each Wick edge consumes a $\beta$ or $\mathsf b$ terminal.
>
> An $n$-loop graph has $n+1$ interaction vertices and $2n+1$ edges.
>
> The root $\bar hA$ supplies one terminal, and $\bar hK$ supplies zero, so
>
> $$
> Q_n(\bar hA)=0\quad(n\ge2),
> \qquad Q_n(\bar hK)=0\quad(n\ge1).
> $$
>
> The first correction cancels the anomaly exactly:
>
> $$
> -\frac{3g^2}{64\pi^2}
> +\frac{3ig}{64\pi^2}(-ig)=0,
> $$
>
> $$
> \mathbf Q_{\mathrm{pert}}
> \left(\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right)=0.
> $$
>
> For $m\ne0$, the zero-jet detector sets positive jets, ghosts, and light quark scalars to zero.
>
> The heavy F-equations give $\bar h=-\bar\mu\bar m$; light mesons and their odd partners remain.
>
> Every differential image evaluates to zero, whereas the cycle evaluates to
>
> $$
> -\bar\mu\bar m\,
> \overline{\widehat M}^{\bar i}{}_{\bar{\tilde\jmath}}
> \beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}\ne0.
> $$
>
> Thus the incoming polynomial class is nonexact in the specified full perturbative Ward complex.
>
> $\square$

The cochain is invariant under constant gauge transformations, and its contraction with a constant ghost shift vanishes.

It therefore belongs to the corresponding relative complex.

<a id="np-4f8eecadb9b2d576-1"></a>

### Local source matching and the zero-mass fiber

Taking the source derivative at fixed magnetic light fields gives

$$
\left.\frac{\partial W_{\mathrm{eff}}}
{\partial\log\widetilde\Lambda^2}\right|_{\widehat q,\widehat{\widetilde q},\widehat M}
=-\frac{\widetilde\Lambda^2}{\mu^4m}\det\widehat M.
$$

This fixes the coefficient of the Wilsonian chiral insertion.

Conditional on microscopic matching in the same local source/contact prescription, the candidate instanton image of the polynomial incoming cycle is

$$
-\bar\mu\bar m\,
\frac{3\overline{\widetilde\Lambda}^{\,2}}{\bar\mu^4\bar m}
\det\overline{\widehat M}
=-\frac{3\overline{\widetilde\Lambda}^{\,2}}{\bar\mu^3}
\det\overline{\widehat M}.
$$

Calculating this microscopic identity remains open.

The matching requires source-dependent D-term descendants, coincident contact terms, quantum products, and a primitive in the finite polynomial jet algebra.

At fixed electric baryons, the source derivative becomes $W_{\mathrm{eff}}$.

Its difference from the fixed-magnetic-variable result is

$$
W_{\mathrm{pert}}
=\widehat q_i\frac{\partial W_{\mathrm{eff}}}{\partial\widehat q_i},
$$

This effective descendant requires matching its local representative together with the current.

At $m=0$, the heavy meson odd letter satisfies

$$
Q_0\beta^{M^4{}_{\tilde4}}=-\frac{\bar h}{\bar\mu},
\qquad
\bar hA=Q_0(-\bar\mu\beta^{M^4{}_{\tilde4}}A).
$$

Secondary ghost-jet data therefore determine the incoming class in the massless limit.

The source cochain $C$ is

$$
C=\bar\mu\beta^{M^4{}_{\tilde4}}
\operatorname{Tr}_{\mathrm{fund}}
[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)]
+\frac{i}{g}\bar hK.
$$

> **Proof: formal massless class.**
>
> Graded differentiation and terminal counting give
>
> $$
> \mathbf Q_{\mathrm{pert}}C
> =-\bar\mu\bar m\operatorname{Tr}_{\mathrm{fund}}
> [(\partial_{+\dot\alpha}\mathsf c_\gamma)
> (\partial_+{}^{\dot\alpha}\mathsf c_\gamma)].
> $$
>
> Thus the cochain is closed at $m=0$.
>
> Set all scalar jets to zero and retain the heavy singlet beta and two first Cartan ghost jets.
>
> Here $t_3$ is the source fundamental Cartan generator, and $\xi_{\dot+},\xi_{\dot-}$ are independent odd detector variables.
>
> Set
>
> $$
> \partial_{+\dot+}\mathsf c_\gamma=t_3^{\mathsf T}\xi_{\dot+},
> \qquad
> \partial_{+\dot-}\mathsf c_\gamma=t_3^{\mathsf T}\xi_{\dot-}.
> $$
>
> All $Q_0$-images, including those of higher ghost jets, vanish through scalar factors or vanishing Cartan structure constants.
>
> The source spinor contraction gives
>
> $$
> C\longmapsto
> 2\bar\mu\operatorname{Tr}_{\mathrm{fund}}(t_3^2)
> \beta^{M^4{}_{\tilde4}}\xi_{\dot+}\xi_{\dot-}\ne0.
> $$
>
> Thus $C$ represents a nonzero $Q_0$-class.
>
> If $C=\mathbf Q_{\mathrm{pert}}U$ with $U$ a formal $\hbar$ series with nonnegative powers, its constant coefficient would imply $C=Q_0U_0$, a contradiction.
>
> Hence $[C]\ne0$ over $\mathbb C[[\hbar]]$.
>
> $\square$

The full light/heavy cross-contraction calculation gives

$$
\left[\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right]_{m=0}
=\frac{3\hbar g^2}{64\pi^2}[C].
$$

The first positive power of $\hbar$ allows a one-loop image of a tree-level closed operator. The complete source calculation below decides this possibility for every power.

The heavy-row alternative also has cross terms.

Retain the calculation's $b=\beta^{M^4{}_{\tilde4}}$, the heavy-row Euler cochain $A_H$, and the preceding quadratic ghost-jet trace $Z$.

Writing the antichiral matching parameter as $\bar\mu$ gives

$$
Q_1(bA_H)
=\frac{g^2}{64\pi^2}bZ
+\frac{ig}{64\pi^2\bar\mu}
(\partial_{+\dot\alpha}t)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)v,
$$

$$
Q_1(\bar hA_H)=-\frac{3g^2}{128\pi^2}\bar hZ.
$$

Both equations satisfy $Q_0Q_1+Q_1Q_0=0$.

The heavy Higgs cross contractions supply the $3/2$ factor.

The independent test below detects the resulting closed class.

The zero-jet/F-term detector establishes perturbative survival of the massless determinant itself.

Assuming the undeformed full electric scalar chiral-ring presentation, operator-algebra Seiberg duality, and the corresponding RG/twist comparison, the magnetic three-by-three minor maps to the electric rank-two polynomial identity and vanishes in the exact SCFT scalar ring.

This conclusion is conditional on the stated operator-ring inputs.

Local matching determines its microscopic instanton order and coefficient.

<a id="massless-all-powers"></a>

### The massless microscopic current survives every loop power

**Theorem 16.3.** In the relative polynomial jet complex of the massless four-flavor magnetic theory, with fixed $g,\mu\ne0$,

$$
\hbar^k[C]\ne0\quad(k=0,1,2,\ldots),\qquad
\left[\bar h A+\frac{3i\hbar g}{64\pi^2}\bar hK\right]_{m=0}\ne0.
$$

Independently permuting the two flavor indices of the heavy meson gives sixteen classes spanning a free rank-sixteen submodule over $\mathbb C[[\hbar]]$.

> **Proof.** A monomial has a letter count and two total derivative orders, one for each holomorphic coordinate. The massless interaction is cubic. An $r$-loop graph changes the letter count by $3(r+1)-2(2r+1)=1-r$ and adds $r$ derivatives in each coordinate by momentum homogeneity. At $r=0$ the tree transformation preserves both derivative orders. [Source, equations B.15 and derivative homogeneity](../sources/all-loop/main.tex)
>
> Retain the coefficient of the heavy singlet odd letter times the two Cartan ghost jets used in the preceding tree detector. The target has three letters and derivative orders $(1,1)$. If $\mathbf Q_{\rm pert}U=\hbar^k C$ with $k\ge1$, its coefficient at order $k$ can receive a positive-loop contribution only from $Q_1$ on the three-letter, zero-derivative part of $U_{k-1}$. Here $U_j$ is the coefficient of $\hbar^j$ in $U$. The detector kills $Q_0 U_k$. At order $k-1$, projecting onto four letters and zero derivatives gives that the indicated part of $U_{k-1}$ is tree closed. This projection applies to arbitrary polynomial jets.
>
> For each flavor pair, assign independent scaling weights to every quark column and antiquark row. Each meson entry carries their inverse product, and each odd dual letter carries the inverse weight of its scalar. The cubic vertices and Wick contractions preserve these algebraic weights. The source has the weight of one flavor-four quark and one flavor-four antiquark. Constant-gauge relativity excludes the undifferentiated ghost and requires color invariance. The complete list of even sources at this weight is as follows. The index $i$ in the first two lines runs over four flavors, and a charged odd partner contracts with the quark color slot.
>
> | Source type | Flavor choices | Number |
> |---|---|---:|
> | Meson odd partner, quark, quark odd partner | heavy singlet times a diagonal quark current, or $\beta^{M^i{}_{\tilde4}}$ times the heavy quark and the dual of flavor $i$ | $4+4-1$ |
> | Meson odd partner, antiquark, antiquark odd partner | heavy singlet times a diagonal antiquark current, or $\beta^{M^4{}_{\tilde\imath}}$ times the heavy antiquark and the dual of flavor $\tilde\imath$ | $4+4-1$ |
> | One meson and two distinct meson odd partners | for meson row $i$, column $\tilde\jmath$, the odd rows have multiset $\{4,i\}$ and the odd columns have multiset $\{\tilde4,\tilde\jmath\}$ | $2\cdot9+3+3$ |
> | Gauge antifield and the heavy pair | the traceless color contraction $t\,\mathsf b\,v$ using the source invariant color form | $1$ |
>
> With two odd letters the weight fixes the three matter types in the first three rows. With zero odd letters the weight requires the heavy pair and one neutral gauge antifield. Two compatible color slots have one invariant pairing, and the traceless part of a column times a row occurs once. These counts include the heavy/light cross contractions.
>
> To read the tree-closure equations, set all mesons, singlet odd letters, ghosts and gauge antifields to zero after applying $Q_0$. Choose the heavy column and row to be the first color basis vector and its dual. For a light flavor $i$, set its quark to the second basis vector and retain the dual component of its odd partner. In a second evaluation exchange quark with antiquark, taking the second dual basis vector and its paired odd component. All other light fields vanish. For the heavy flavor use the first color in both evaluations. In each case take the coefficient of the retained odd variable.
>
> The tree variation of the heavy meson odd letter is $-1/\bar\mu$ at these values. The mixed meson contraction for a light flavor vanishes by the vanishing first–second color pairing. Terms with two singlet odd letters leave a singlet odd letter. The gauge antifield produces the traceless matter moment from (B.15), with opposite signs for a fundamental column and a dual row. Their two evaluations cancel, including its trace subtraction. If $a_i,\widetilde a_i$ are the coefficients of the heavy singlet odd letter placed to the left of the diagonal quark and antiquark currents in the source under test, tree closure therefore requires
>
> $$
> 0=-\frac{a_i+\widetilde a_i}{\bar\mu}
> \qquad(i=1,2,3,4).
> $$
>
> The one-loop target contains the heavy singlet odd letter. No interaction vertex supplies this letter. It must remain on the root, so the surviving graph on the two remaining root letters is the fundamental triangle. A gauge vertex cannot absorb a meson, and a nondiagonal current cannot close its flavor line. The dual-representation triangle has the same quadratic trace as the fundamental. Moving the odd variation past the heavy singlet gives the common minus sign. With $b=\beta^{M^4{}_{\tilde4}}$ and $Z$ the ordered ghost trace above, the projected one-loop coefficient multiplying the Cartan evaluation of $bZ$ is
>
> $$
> -\frac{g^2}{64\pi^2}
> \sum_{\substack{i=1,\ldots,4\\\mathrm{quark\ flavors}}}
> (a_i+\widetilde a_i)=0.
> $$
>
> The same projection of $C$ is $2\bar\mu\operatorname{Tr}_{\rm fund}(t_3^2)b\,\xi_{\dot+}\xi_{\dot-}\ne0$. This contradicts the coefficient equation for every $k\ge1$. The preceding tree argument covers $k=0$. Each nonzero scalar series is a unit times a power of $\hbar$. Independent row and column permutations give the other fifteen components, with distinct preserved flavor weights, proving the module statement. The earlier polynomial identity gives the stated nonzero incoming class. $\square$

Expanding the tree differential in both colors gives $142$ distinct output monomials from $39$ invariant sources. The tree matrix and the matrix with the triangle row appended both have rational rank $38$. The calculation also verifies the four separate evaluation identities on each of the $39$ columns. [Exact source matrix](../checks/massless_source.py)

<a id="mass-pole-ward"></a>

### A necessary mass pole in microscopic matching

**Proposition 16.4.** Suppose the local supercharge acts linearly over the parameter ring and its first magnetic-instanton coefficient and perturbative boundary comparisons admit nonnegative loop powers. Let $D$ denote its coefficient at $\hbar^0$, with the instanton weight stripped off. Matching the proposed image at nonzero $m$ forces the scalar projection of the light-row current $A$ to obey

$$
\operatorname{ev}_{\bar m}(D A)
=\frac{3}{\bar\mu^4\bar m}\det\overline{\widehat M}.
$$

In particular, a mass-regular polynomial action on $A$ through the massless fiber is incompatible with this proposed matching in those loop conventions.

> **Proof.** Nilpotency gives $Q_0 D+D Q_0=0$. Off shell at tree level,
>
> $$
> Q_0(-\bar\mu b A)=(\bar h+\bar\mu\bar m)A,\qquad
> \bar h A=Q_0(-\bar\mu b A)-\bar\mu\bar m A .
> $$
>
> Evaluate zero-derivative scalars on the algebraic family with heavy column the first color basis vector, heavy row $-\bar\mu\bar m$ times its dual, and arbitrary light meson block $\overline{\widehat M}$. Set every other scalar, odd letter, gauge field and positive jet to zero. Denote this substitution by $\operatorname{ev}_{\bar m}$. All tree F-equations hold, including the heavy one. Equation (B.15) then makes the substitution zero on every tree boundary. It is polynomial in $\bar m$ all the way to $\bar m=0$.
>
> The $\hbar^0$ term of the polynomial source is $\bar h A$. Applying $D$ to its tree identity and evaluating gives
>
> $$
> \operatorname{ev}_{\bar m}\bigl(D(\bar h A)\bigr)
> =-\bar\mu\bar m\,\operatorname{ev}_{\bar m}(D A).
> $$
>
> The proposed determinant coefficient on the left is $-3\det\overline{\widehat M}/\bar\mu^3$. This proves the residue. A tree-boundary change of either local representative evaluates to zero at each nonzero mass, even if its primitive uses inverse powers of the mass. A mass-regular $D A$ would give a mass-divisible left-hand side, in contradiction with evaluation of the independent meson block at the identity. $\square$

The physical insertion on the open Higgs region may contain an inverse heavy contraction. The test fixes its required scalar residue under the stated loop expansion. Determining its coefficient from local instanton insertions and coincident subtractions is the remaining matching calculation.

<a id="mass-loop-base-change"></a>

### Mass and loop order together

Keep the antichiral mass and loop parameter independent, with coefficients $R=\mathbb C[[\bar m,\hbar]]$, holding $g,\mu\ne0$. Write $H_R$ for the cohomology of the relative polynomial jet complex over $R$. The ghost trace $Z$ and the two odd cochains $A,C$ obey

$$
\begin{gathered}
\mathbf Q_{\rm pert}Z=0,\qquad
\mathbf Q_{\rm pert}A=-\frac{3\hbar g^2}{64\pi^2}Z,\qquad
\mathbf Q_{\rm pert}C=-\bar\mu\bar m Z,\\
\bar h A+\frac{3i\hbar g}{64\pi^2}\bar hK
=-\mathbf Q_{\rm pert}(\bar\mu b A)
+\frac{3\hbar g^2}{64\pi^2}C-\bar\mu\bar m A .
\end{gathered}
$$

**Theorem 16.5.** The scalar annihilator of $[Z]$ in $R$ is exactly $(\bar m,\hbar)$. The scalar annihilator of the polynomial incoming class is zero. This class gives a nonzero element in the kernel of the natural specialization map

$$
H_R/(\bar m,\hbar)H_R\ \longrightarrow\
H\!\left(\left.\mathbf Q_{\rm pert}\right|_{\bar m=\hbar=0}\right).
$$

The left-hand side takes cohomology first and sets both parameters to zero afterwards. The right-hand side first sets both coefficients to zero in the differential.

> **Proof.** At zero mass and loop order, the pure Cartan ghost evaluation gives $Z\mapsto2\operatorname{Tr}_{\rm fund}(t_3^2)\xi_{\dot+}\xi_{\dot-}\ne0$ and kills all tree boundaries. Thus $[Z]\ne0$ over $R$. Its two displayed primitives kill it after multiplication by the two generators of the maximal ideal. Every series outside that ideal is a unit. For the incoming current, extend the mass-polynomial evaluation of Proposition 16.4 by keeping the odd partners of the light meson block. This again kills all tree images because the retained odd letters are gauge singlets. Every positive-loop image has positive derivative order. The current evaluates to $-\bar\mu\bar m\,\overline{\widehat M}^{\bar i}{}_{\bar{\tilde\jmath}}\beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}$. Its product with any nonzero scalar series is nonzero. This proves both annihilator statements.
>
> On a massless formal class take its constant loop coefficient and extract the coefficient of $b\xi_{\dot+}\xi_{\dot-}$ in the Cartan evaluation of Theorem 16.3. This is well defined modulo boundaries. Two vanishing properties follow.
>
> 1. For any class $[T]$ killed by $\hbar$, write $\mathbf Q_{\rm pert}U=\hbar T$ at zero mass. The tree part of $U$ is closed. The four flavor evaluations of Theorem 16.3 set its one-loop triangle projection to zero, and the tree image of the next coefficient also evaluates to zero. This proves the assertion for arbitrary $T$.
> 2. For any class extending to a closed cochain over $R$, write $U_{00},U_{10}$ for its coefficients of $1,\bar m$ at loop order zero. The mass term in the tree variation acts as $-\bar m\,\partial/\partial b$ because the spacetime derivatives of the mass vanish. Closure in the next mass coefficient gives
>
> $$
> \left.Q_0\right|_{\bar m=0}U_{10}=\frac{\partial U_{00}}{\partial b}.
> $$
>
> Evaluate on the same two Cartan jets, now discarding $b$. The left side is a tree boundary. The coefficient of $\xi_{\dot+}\xi_{\dot-}$ on the right extracts precisely the preceding $b\xi_{\dot+}\xi_{\dot-}$ coefficient. Thus the projection also vanishes on restrictions of extendable classes.
>
> If the incoming class were $\bar m[V]+\hbar[U]$ in $H_R$, its zero-mass restriction would give
>
> $$
> \hbar\left(\frac{3g^2}{64\pi^2}[C]-[U|_{\bar m=0}]\right)=0.
> $$
>
> Test 1 kills the projection of the parentheses, and test 2 kills the projection of $[U|_{\bar m=0}]$. The projection of the remaining summand is $3g^2(2\bar\mu\operatorname{Tr}_{\rm fund}t_3^2)/(64\pi^2)\ne0$. Thus the incoming class lies outside $(\bar m,\hbar)H_R$. Setting both parameters to zero directly gives $\bar h A=Q_0(-\bar\mu bA)$, proving the kernel statement. $\square$

> **Explicit second differential for joint specialization.** The Koszul resolution implements setting two parameters to zero using exterior generators $e_m,e_\hbar$ with variations $\bar m,\hbar$, respectively. The two parameters form a regular sequence on the jet coefficient modules, where each acts on the formal coefficients of the monomials. Place the exterior generators on the left in the displayed order. The total differential acts on a product by the Koszul variation of the exterior factors plus $(-1)^{\text{number of exterior factors}}$ times the physical variation. Directly,
>
> $$
> \begin{aligned}
> D_{\rm tot}\!\left(e_m e_\hbar Z-\frac{e_\hbar C}{\bar\mu}
> +\frac{64\pi^2}{3g^2}e_m A\right)
> &=-\frac{\hbar C}{\bar\mu}+\frac{64\pi^2\bar m}{3g^2}A .
> \end{aligned}
> $$
>
> The terms with one exterior generator cancel by the two primitive equations. Taking physical cohomology first, $[Z]$ at exterior degree two survives the first Koszul differential because both parameters kill it. There is no exterior degree three. The displayed residual, in exterior degree zero, represents $-64\pi^2/(3g^2\bar\mu)$ times the polynomial incoming class modulo $(\bar m,\hbar)H_R$, nonzero by Theorem 16.5. It is therefore a nonzero second differential in the specialization calculation. Its order counts two Koszul steps in the perturbative theory. $\square$
