<a id="np-review"></a>

## 3. Review: operator formalism and required background

The calculations use three descriptions of the same proposed operator relation: polynomials in microscopic fields, classes under a square-zero transformation, and an effective theory obtained after eliminating massive fields.

This chapter gives the maps and equations connecting those descriptions.

<a id="np-review-003"></a>

### 3.1. Field polynomials and the square-zero transformation

In the holomorphic description of [Budzik–Gaiotto–Kulp–Williams–Wu–Yu](https://arxiv.org/abs/2306.01039), $b$ has even coefficients, $c$ has odd coefficients, and the ghost $c$ records gauge transformations.

For an expression $A$, its parity $|A|$ is zero for even degree and one for odd degree; the tree transformation obeys the graded derivation rule

$$
Q_0(AB)=(Q_0A)B+(-1)^{|A|}A(Q_0B),
\qquad
Q_0c=c^2,
\qquad Q_0b=cb-bc.
$$

Matrix multiplication retains its displayed order; the sign concerns the even/odd coefficients.

> **Proof of square-zero action on the generators.**
>
> $$
> Q_0^2c=(c^2)c-c(c^2)=0,
> $$
>
> $$
> \begin{aligned}
> Q_0^2b
> &=(c^2)b-c(cb-bc)-(cb-bc)c-b(c^2)\\
> &=c^2b-c^2b+cbc-cbc+bc^2-bc^2=0.
> \end{aligned}
> $$
>
> Since $[Q_0,\partial_1]=[Q_0,\partial_2]=0$, the same square-zero action holds on every derivative and hence on their polynomial algebra.
>
> $\square$

Gauge-invariant traces select operators unchanged by the gauge transformations encoded by $c$.

The cohomology identifies $Q_0$-closed expressions differing by a $Q_0$-variation; an explicit primitive demonstrates that a class is zero.

<a id="np-review-015"></a>

### 3.2. From tree classes to successive loop equations

The loop parameter $\hbar$ is formal: equality of two series means equality of every coefficient.

Expanding the condition $Q_{\mathrm{pert}}^2=0$ gives

$$
Q_0^2=0,
\qquad Q_0Q_1+Q_1Q_0=0,
\qquad Q_1^2+Q_0Q_2+Q_2Q_0=0.
$$

The first correction induces $d_1[\mathcal O]=[Q_1\mathcal O]$ on tree cohomology.

> **Proof that the first induced differential is well defined.**
>
> For $Q_0\mathcal O=0$, the second coefficient equation gives
>
> $$
> Q_0(Q_1\mathcal O)=-Q_1Q_0\mathcal O=0.
> $$
>
> Changing the tree representative by $Q_0x$ changes its image by a tree boundary:
>
> $$
> Q_1(\mathcal O+Q_0x)
> =Q_1\mathcal O-Q_0(Q_1x).
> $$
>
> The next coefficient equation gives
>
> $$
> d_1^2[\mathcal O]=[Q_1^2\mathcal O]
> =[-Q_0Q_2\mathcal O-Q_2Q_0\mathcal O]=0.
> $$
>
> Thus both the induced map and its cohomology are defined on classes.
>
> $\square$

Solving the first correction equation removes its outgoing obstruction:

$$
Q_0(\text{first correction})=-Q_1\mathcal O.
$$

The quotient by incoming $d_1$ images also removes tree classes which become exact at this stage.

The resulting $E_2$ is therefore $\ker d_1/\operatorname{im}d_1$.

**Textbook input: filtered-complex spectral sequence** (Weibel, *An Introduction to Homological Algebra*, chapter “Spectral Sequences”).

For a differential preserving a finite filtration, successive equations induce square-zero maps $d_r$ on $E_r$, with $E_{r+1}=H(E_r,d_r)$; the stable page is the associated graded space of the full cohomology.

“Associated graded” records successive quotients of the filtration, rather than a preferred representative or multiplication rule.

Here fixed weight makes the filtration finite; the source specifies its field counts and loop shifts.

The formal coefficient ring and the computed operator degree are retained separately, since division by $\hbar$ changes which primitives are allowed.

<a id="np-review-038"></a>

### 3.3. Quantum products and secondary classes

Higher-loop operators act on several letters at once, so the product compatible with the full differential is the source quantum product $\star$.

Its defining compatibility is

$$
Q_{\mathrm{pert}}(A\star B)
=(Q_{\mathrm{pert}}A)\star B
+(-1)^{|A|}A\star(Q_{\mathrm{pert}}B).
$$

**External input:** the [source product theorem](../sources/all-loop/theorems.md#b-00513) supplies an associative, graded-commutative product with identity for the specified anomaly-free perturbative pure-gauge complex.

Anomaly-free means that gauge transformations preserve the quantum integration measure; the [product appendix](24-product-kernels.md#np-section-22) states the source interaction and applies the theorem in that setting.

Associativity means independence of bracket placement; graded commutativity supplies a minus sign when two odd classes are interchanged; the identity leaves every class unchanged.

> **Proof that the product respects class equivalence.**
>
> For closed $A,B$, the compatibility equation gives
>
> $$
> Q_{\mathrm{pert}}(A\star B)=0.
> $$
>
> If $A$ is changed by $Q_{\mathrm{pert}}x$, then
>
> $$
> (A+Q_{\mathrm{pert}}x)\star B-A\star B
> =(Q_{\mathrm{pert}}x)\star B
> =Q_{\mathrm{pert}}(x\star B).
> $$
>
> Hence the product class remains unchanged.
>
> $\square$

An identity between classes becomes usable in a further construction once its primitive is known.

Write $L\mathcal O=C_0\star\mathcal O$, and let the odd operator $h$ satisfy $Q_{\mathrm{pert}}h+hQ_{\mathrm{pert}}=L^N$.

For a closed representative $a$ and a primitive $w$ of $La$, the secondary construction is

$$
Q_{\mathrm{pert}}a=0,
\qquad Q_{\mathrm{pert}}w=La,
\qquad [ha-L^{N-1}w].
$$

> **Proof of closure and the stated primitive ambiguity.**
>
> Since $Q_{\mathrm{pert}}L=LQ_{\mathrm{pert}}$,
>
> $$
> \begin{aligned}
> Q_{\mathrm{pert}}(ha-L^{N-1}w)
> &=L^Na-hQ_{\mathrm{pert}}a-L^{N-1}Q_{\mathrm{pert}}w\\
> &=L^Na-L^{N-1}La=0.
> \end{aligned}
> $$
>
> For another primitive $w'$, one has $Q_{\mathrm{pert}}(w'-w)=0$ and
>
> $$
> (ha-L^{N-1}w')-(ha-L^{N-1}w)
> =-L^{N-1}(w'-w).
> $$
>
> At fixed $h$, the secondary class is therefore defined modulo $C_0^{N-1}H(Q_{\mathrm{pert}})$.
>
> $\square$

The exact continuation replaces the homotopy equation by $Qh+hQ=L^N-\nu\Lambda^{3N}\mathbf1$, where $Q$ includes nonperturbative terms and $\mathbf1$ is the identity operator.

The [conditional pairing theorem](09-secondary-pairing.md#np-93111f11b210145b-1) computes the action induced on perturbative classes; the physical comparison supplies this exact identity and its scalar coefficient $\nu$.

<a id="np-review-065"></a>

### 3.4. Trace invariants, finite-dimensional counting, and the Cartan test

The invariant traces supply polynomial coefficients for the leading classes.

Their tree closure follows directly from the gauge action:

> **Proof.**
>
> The derivation rule gives $Q_0(b^k)=[c,b^k]$.
>
> Cyclicity of the matrix trace, with $b$ even, gives
>
> $$
> Q_0\operatorname{tr}(b^k)
> =\operatorname{tr}(cb^k-b^kc)=0.
> $$
>
> Thus multiplication by these traces preserves tree cohomology.
>
> $\square$

**Textbook input: Cayley–Hamilton theorem** (Hoffman–Kunze, *Linear Algebra*, “Annihilating Polynomials”).

A square matrix satisfies its characteristic polynomial.

For the traceless $2\times2$ matrix $b$, with $\mathbf1$ the identity matrix, its usable form is

$$
b^2=\frac{\operatorname{tr}(b^2)}2\mathbf1,
\qquad
\operatorname{tr}(b)=0.
$$

The [leading-module theorem](05-leading-modules.md#np-section-03) establishes

$$
H(Q_0)_{\mathrm{leading\ scalar}}
=\mathbb Q[\operatorname{tr}(b^2)]\,[C_0]
\qquad(SU(2)).
$$

A module allows polynomial coefficients in the displayed traces; a vector space over $\mathbb Q$ allows rational-number coefficients.

Consequently, rank one over this polynomial ring includes infinitely many independent vectors over $\mathbb Q$:

$$
[C_0],\quad[\operatorname{tr}(b^2)C_0],\quad
[\operatorname{tr}(b^2)^2C_0],\quad
[\operatorname{tr}(b^2)^3C_0],\ldots
$$

At fixed $b$-degree $k$, exactly one displayed power has that degree when $k$ is even; none has odd degree.

Fixed weight bounds derivative orders and field counts.

The resulting chain spaces are finite-dimensional, so closure and exactness are kernel and image calculations with rational matrices.

Rank–nullity gives

$$
\dim H(Q_0)
=\dim(\text{specified chain space})
-\operatorname{rank}(\text{outgoing }Q_0)
-\operatorname{rank}(\text{incoming }Q_0).
$$

The [first-correction test](04-complex-and-products.md#np-f1ccd4f158d90ace-1) asks whether a tree cocycle $\mathcal O$ admits a polynomial solution of

$$
Q_0(\text{first correction})=-Q_1\mathcal O.
$$

Cartan restriction keeps traceless diagonal field matrices and their jets, setting all off-diagonal components to zero; these diagonal matrices form $\mathfrak t$.

> **Proof of the necessary Cartan condition.**
>
> Each diagonal ghost entry is odd, so its square vanishes.
>
> Diagonal $b$ entries commute with it.
>
> Hence
>
> $$
> \left.Q_0c\right|_{\mathfrak t}
> =\left.c^2\right|_{\mathfrak t}=0,
> \qquad
> \left.Q_0b\right|_{\mathfrak t}
> =\left.[c,b]\right|_{\mathfrak t}=0.
> $$
>
> Derivatives commute with $Q_0$ and restriction.
>
> The derivation rule therefore annihilates the restriction of every polynomial boundary:
>
> $$
> \left.Q_0A\right|_{\mathfrak t}=0,
> \qquad
> Q_0(\text{first correction})=-Q_1\mathcal O
> \ \Longrightarrow\ 
> \left.Q_1\mathcal O\right|_{\mathfrak t}=0.
> $$
>
> $\square$

A nonzero restricted image excludes a first correction.

A zero restricted image leaves the full polynomial equation to solve.

The $SU(2),SU(3)$ all-degree results establish $E_2$ in their stated leading sectors; the $SU(4)$ classification establishes the Cartan-test kernel.

Full survival requires every subsequent obstruction to vanish and the class to survive incoming differentials.

The [complete representatives](06-su3-lifts.md#np-93544cd67af7d9c7-1) and [higher-partner existence proof](04-complex-and-products.md#np-d40c1e0d3b4949c1-1) supply those additional calculations in their specified sectors.

<a id="np-review-089"></a>

### 3.5. Effective interactions and scalar operator relations

The standard supersymmetric field-theory input is that a holomorphic superpotential contributes a chiral superspace integral and gives, in the conjugate scalar coordinates used below, the F-equations

$$
\int d^4x\,d^2\theta\,W,
\qquad
\frac{\partial\bar W}{\partial\gamma_i}=0.
$$

Here $d^2\theta$ integrates over two anticommuting supersymmetry coordinates.

In the Wilsonian description, which integrates out field modes above a chosen energy scale, the superpotential written in holomorphic fields and couplings receives no perturbative corrections; wave-function renormalization changes the kinetic terms.

These inputs are stated in Weinberg, *The Quantum Theory of Fields*, Volume III, *Supersymmetry*.

The relation to the chosen twisted operator complex additionally uses its specified field dictionary.

For the F-equation part of the effective antichiral complex used here, the commuting coordinates $\gamma_i$ and odd dual partners $\beta^i$ obey

$$
Q\gamma_i=0,
\qquad
Q\beta^i=-\frac{\partial\bar W}{\partial\gamma_i}.
$$

The bar denotes the conjugate superpotential expressed in these conjugate scalar coordinates.

For polynomial $\bar W$, assigning degree zero to $\gamma_i$ and degree minus one to $\beta^i$ gives the scalar cohomology before restricting to gauge-invariant classes

$$
H^0(Q)=\frac{\mathbb C[\gamma]}{(\partial\bar W)}.
$$

The denominator denotes polynomial combinations of all F-equations.

This statement concerns the specified F-equation complex with zero spacetime derivatives; a remaining gauge symmetry also requires restricting to invariant classes.

> **Proof: the quotient and its primitive.**
>
> Every scalar polynomial $\mathcal O(\gamma)$ is closed.
>
> For each coordinate label $i$,
>
> $$
> Q\bigl(\mathcal O(\gamma)\beta^i\bigr)
> =-\mathcal O(\gamma)\frac{\partial\bar W}{\partial\gamma_i}.
> $$
>
> Every degree-minus-one expression is a polynomial combination of these odd partners.
>
> Its image is therefore exactly the displayed F-equation ideal.
>
> In the ADS example, $\gamma^{(0)}=\bar M$ and
>
> $$
> Q\beta^{(0)}=\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2},
> \qquad
> Q\bigl((\gamma^{(0)})^2\beta^{(0)}\bigr)=\bar\Lambda_L^5.
> $$
>
> At fixed $\bar\Lambda_L^5\ne0$, division gives an odd primitive of $1$.
>
> Multiplying that primitive by any closed operator makes it exact.
>
> $\square$

Over formal coefficients $\mathbb C[[\bar\Lambda_L^5]]$, the inverse scale is absent.

On the region $\gamma^{(0)}\ne0$, the same calculation gives

$$
H^{-1}(Q)=0,
\qquad
H^0(Q)=\mathbb C[(\gamma^{(0)})^{\pm1}],
\qquad
\bar\Lambda_L^5H^0(Q)=0.
$$

> **Proof of the formal cohomology.**
>
> Multiplication by $\bar\Lambda_L^5$ is injective before taking cohomology, and $\gamma^{(0)}$ is invertible on this domain.
>
> $$
> Q(\mathcal O\beta^{(0)})=\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2}\mathcal O=0
> \quad\Longrightarrow\quad\mathcal O=0.
> $$
>
> Hence degree-minus-one cohomology vanishes.
>
> The degree-zero quotient sets $\bar\Lambda_L^5$ equal to zero and leaves the Laurent polynomials in $\gamma^{(0)}$.
>
> $\square$

Thus the formal family has parameter torsion, while the separately specified model at fixed nonzero scale has vanishing cohomology.

Allowing inverse scalar coordinates and allowing inverse parameters are distinct choices.

For a parameter $\eta$ in a scalar quotient by an ideal $I$, a polynomial $f$ gives zero after adjoining $\eta^{-1}$ precisely when $\eta^r f\in I$ for some nonnegative integer $r$. The set of all such polynomials is called the saturation $I:\eta^\infty$. Taking its quotient by $I$ retains the actual classes killed by localization. The least such exponent for a fixed class measures how many parameter factors its primitive requires. The successive differential on page $r$ is the obstruction obtained by solving the first $r-1$ coefficient equations and taking the coefficient at order $r$, modulo the images from earlier orders. The determinant and multiplier calculations apply these operations before choosing a nonzero numerical coupling. [Saturation and finite-sector proof](23-operator-data.md#complete-tests)

An unbroken translation acts also on these successive pages. At first order write $Q_0,Q_1$ for the differential coefficients and $P_0,P_1$ for the corresponding translation coefficients. The tree relation is $[Q_0,P_0]=0$. The first coefficient of the exact Ward identity gives $[Q_1,P_0]=-[Q_0,P_1]$. Applied to a tree cycle, the right-hand side is a tree boundary. Thus the induced differential commutes with the induced translation. At higher orders a filtered translation carries each solution of the preceding coefficient equations to another solution, and carries its obstruction to the translated obstruction. [Bulk application](21-seiberg-witten.md#bulk-translation)

<a id="np-review-111"></a>

### 3.6. Duality and elimination of massive fields

[Seiberg duality](https://arxiv.org/abs/hep-th/9411149) relates electric and magnetic operator descriptions.

For $SU(N_c)$ SQCD in the stated flavor range, the electric composite meson is an elementary magnetic field:

$$
M^i{}_{\tilde\jmath}=Q^{ia}\widetilde Q_{a\tilde\jmath},
\qquad
W_{\mathrm{mag}}=\frac1\mu M^i{}_{\tilde\jmath}q_i{}^A\widetilde q_A{}^{\tilde\jmath}.
$$

Here $a,A$ are electric and magnetic color slots, and $i,\tilde\jmath$ label flavors; $Q^{ia}$ is a quark scalar.

Magnetic strong dynamics restores electric composite relations.

Transporting an effective F-equation into the original local supercharge complex requires an operator map that also matches products and terms supported at coincident insertion points.

This remains the missing step for the proposed microscopic determinant pairing.

Eliminating a massive field means solving its stationary equation and substituting its solution into $W$.

For the two-flavor $SU(2)$ constrained theory, integrating out the flavor of mass $m_2$ gives

$$
W_{\mathrm{eff}}=\frac{\Lambda_L^5}{M},
\qquad
\Lambda_L^5=m_2\Lambda^4,
\qquad M\ne0.
$$

> **Proof: the denominator domain.**
>
> After adding $m_1M$, with $m_1\ne0$ and $\Lambda_L^5\ne0$, the stationary equation is
>
> $$
> \frac{\partial W_{\mathrm{eff}}}{\partial M}
> =\frac{m_1M^2-\Lambda_L^5}{M^2}=0,
> \qquad M\ne0.
> $$
>
> Its numerator has two nonzero roots.
>
> Multiplication by an additional $M$ gives the polynomial $M(m_1M^2-\Lambda_L^5)$, whose zero at $M=0$ lies outside the original domain.
>
> The correct algebra retains the inverse meson:
>
> $$
> \frac{\mathbb C[M,M^{-1}]}{(m_1M^2-\Lambda_L^5)}.
> $$
>
> The antichiral scalar cohomology uses its conjugate presentation.
>
> $\square$

Fractional powers arise when the retained condensate $S$ is eliminated from

$$
S^{N_c-N_f}\det M=\Lambda^{3N_c-N_f},
\qquad W_{\mathrm{eff}}=(N_c-N_f)S.
$$

The polynomial equation keeps every condensate branch.

Finally, [Beasley–Witten](https://arxiv.org/abs/hep-th/0409149) exhibit higher F-terms containing extra fermions.

On the smooth Higgs branch of $SU(2)$ SQCD with four flavors, the chiral integrand contains $p=N_f-N_c+1=3$ fermion pairs and supplies an eight-fermion component interaction.

Such terms carry information beyond the scalar equations derived from an ordinary $W$.

After a mass deformation to three flavors, the six-fermion effective interaction coexists with the meson superpotential $W=-\operatorname{Pf}V/\Lambda^3$; both describe the same theory at the specified field scales.

<a id="np-review-127"></a>

### 3.7. Instanton integration and local operator matching

An instanton is a finite-action solution of the Euclidean field equations in a sector of nonzero integer gauge winding.

The Euclidean action is the spacetime integral determining the weight of a field configuration; integer winding distinguishes configurations that cannot be continuously deformed into one another with the prescribed boundary behaviour.

In the coupling convention used for the instanton phases, $g$ is the gauge coupling and $\theta$ is the vacuum angle:

$$
\tau=\frac{\theta}{2\pi}+\frac{4\pi i}{g^2},
\qquad
e^{2\pi i\tau}=e^{-8\pi^2/g^2+i\theta}.
$$

The real exponential is invisible to a power series about zero coupling:

$$
\lim_{g\to0^+}g^{-r}e^{-8\pi^2/g^2}=0,
\qquad r\ge0.
$$

The source instanton scale factors, such as $\Lambda^{3N}$, also include their prescribed dimensions and measure normalization.

Collective coordinates parametrize the instanton’s position, size, gauge orientation, and fermion zero modes; a zero mode solves the homogeneous linearized field equation.

Integration over an odd coordinate extracts its coefficient, so a nonzero fermion integral requires the indicated ordered product of all odd coordinates.

**Textbook inputs: change of variables, integration by parts, and the supersymmetry Ward identity** (Weinberg, *The Quantum Theory of Fields*, Volumes II–III).

The Ward identity is the equality obtained by making an infinitesimal supersymmetry change of variables in the path integral.

The path integral sums over field configurations with their action-dependent weights.

At finite cutoff, integration by parts keeps the boundary flux and the variation of the integration density and regulator.

The source calculation implements those terms in the instanton coordinates; vanishing divergence and vanishing boundary flux are separate conditions.

A Wilson line transports a color-carrying factor to the color frame of another insertion before the indices are contracted.

Its supersymmetry variation contributes to the point-split operator, whose factors are kept at distinct points before a prescribed coincident limit is taken.

The physical derivative is gauge covariant, meaning that it transforms in the same representation as the differentiated field.

With the source connection variation $R_\alpha$, its normalization is fixed by

$$
\delta_R(D_\alpha b)=[R_\alpha,b].
$$

Here $D_\alpha$ is the source covariant derivative, $R_\alpha$ the residual gauge-connection variation after collective motion is removed, and the bracket is the matrix commutator.

The completed connection insertion is one term in the local Ward identity; the other terms come from the remaining field variations, action, integration measure, and coincident-point contributions.

The distinction between a separated correlator and a local product is visible in the source $SU(2)$ equations.

Let $[C_0^2]_{\rm ren}$ be the local product after its short-distance singularities have been subtracted; $x_2$ is the specified odd primitive, and $\nu_{\rm prod},\nu_Q$ are the two scalar coefficients in the same prescription:

$$
C_0\star C_0=[C_0^2]_{\rm ren}+\nu_{\rm prod}\Lambda^6,
\qquad
Qx_2=[C_0^2]_{\rm ren}+\nu_Q\Lambda^6.
$$

> **Proof of the invariant local coefficient.**
>
> Subtracting the two equations gives
>
> $$
> C_0\star C_0-Qx_2
> =(\nu_{\rm prod}-\nu_Q)\Lambda^6.
> $$
>
> A finite change of subtraction by $s\Lambda^6$, with $s$ a scalar coefficient, replaces the common local product and both coefficients by
>
> $$
> [C_0^2]_{\rm ren}\longmapsto[C_0^2]_{\rm ren}+s\Lambda^6,
> \qquad
> \nu_{\rm prod}\longmapsto\nu_{\rm prod}-s,
> \qquad
> \nu_Q\longmapsto\nu_Q-s.
> $$
>
> Therefore
>
> $$
> (\nu_{\rm prod}-s)-(\nu_Q-s)=\nu_{\rm prod}-\nu_Q.
> $$
>
> $\square$

The physical matching calculation fixes this difference in the stated local operator dictionary, together with the remaining contributions to the secondary insertion.

The subsequent chapters separately calculate the perturbative classes, the conditional pairing identity, the evaluated instanton integrals, and the exact effective interactions used in that comparison.
