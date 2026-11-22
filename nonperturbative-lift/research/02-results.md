<a id="np-summary-results"></a>

## 2. Summary of Results

<a id="np-summary-001"></a>

### 2.1. Operators, equivalence classes, and perturbative order

The first calculations use field polynomials at one spacetime point, including finitely many derivatives; these are the polynomial local operators considered here.

In the field-coefficient algebra, even variables commute and odd variables anticommute with other odd variables.

The chosen supersymmetry transformation is denoted by $Q$; it reverses this even/odd grading, and applying it twice gives zero.

$$
Q^2=0.
$$

An operator $\mathcal O$ is **closed** when $Q\mathcal O=0$, and **exact** when it equals $Qx$ for an allowed local operator $x$.

The operator $x$ is called a **primitive** of $Qx$.

The class $[\mathcal O]$ identifies closed operators differing by an exact operator; **cohomology**, denoted by $H(Q)$, is the vector space of these classes.

Here $\ker Q$ is the space of operators sent to zero, and $\operatorname{im}Q$ is the space of values $Qx$.

$$
Q\mathcal O=0,
\qquad
[\mathcal O+Qx]=[\mathcal O],
\qquad
[Qx]=0,
\qquad
H(Q)=\frac{\ker Q}{\operatorname{im}Q}.
$$

The formal parameter $\hbar$ counts loops in perturbation theory; $Q_r$ is the coefficient at loop order $r$, and $Q_0$ is the tree-level transformation.

$$
Q_{\mathrm{pert}}
=\sum_{\substack{r\ge0\\\mathrm{loop\ orders}}}\hbar^rQ_r,
\qquad
Q_{\mathrm{pert}}^2=0.
$$

The spaces $E_1,E_2,E_3$ record classes after successively imposing the tree, first-loop, and second-loop equations; $d_r$ is the correction induced on $E_r$, and $E_\infty$ records survival through all orders in the finite sector under consideration.

$$
E_1=H(Q_0),
\qquad
d_1[\mathcal O]=[Q_1\mathcal O],
\qquad
E_{r+1}=\frac{\ker d_r}{\operatorname{im}d_r}.
$$

The quantum product $A\star B$ multiplies local operators with their prescribed loop corrections; $\mu_r(A,B)$ is its coefficient at loop order $r$.

$$
A\star B=AB+
\sum_{\substack{r\ge1\\\mathrm{loop\ orders}}}\hbar^r\mu_r(A,B),
\qquad
[A][B]=[A\star B].
$$

The classes in each product are closed, and the product respects the equivalence relation above.

The gauge group $SU(N)$ consists of unitary $N\times N$ matrices with determinant one; gauge-invariant operators retain their value under a change of gauge description.

The holomorphic twist is the formulation retaining the chosen $Q$ and derivatives in two complex spacetime directions.

A jet means a field derivative evaluated at the insertion point.

Its field $b$ has commuting coefficients, while the ghost $c$ has anticommuting coefficients encoding gauge transformations; $\partial_1,\partial_2$ differentiate in those two directions, and $\operatorname{tr}$ is the matrix trace, the sum of diagonal entries.

The action on these field and ghost variables includes the gauge-symmetry term:

$$
Q_0c=c^2,
\qquad
Q_0b=[c,b]=cb-bc,
\qquad
C_0=\operatorname{tr}(\partial_1c\,\partial_2c),
\qquad
Q_{\mathrm{pert}}C_0=0.
$$

The operator $C_0$ represents the bilinear in the gaugino, the fermionic partner of the gauge field, in the source normalization.

In the pure-gauge formulas, $N_b,N_c$ count the numbers of $b,c$ factors; the weight $(n_1,n_2)$ counts derivatives in the two directions and assigns an additional $(1,1)$ to each $b$ factor.

A scalar sector consists of combinations invariant under the stated spacetime rotations; a leading two-ghost scalar has two $c$ factors in its first nonzero perturbative coefficient.

The result statements below give the computed relations and their domains; the linked chapters contain their proofs.

<a id="np-summary-028"></a>

### 2.2. The $SU(2)$ product and its nilpotent class

Here $[\hbar^r]$ extracts a loop coefficient, and $\operatorname{tr}$ is the matrix trace.

The subscripts in $C_0,C_2,C_4$ label the displayed families; the gauge group remains fixed within each statement.

$$
C_0=\operatorname{tr}(\partial_1c\,\partial_2c),
\qquad
[\hbar^0]C_2=-\operatorname{tr}(b^2)C_0,
\qquad
[\hbar^0]C_4=\operatorname{tr}(b^2)^2C_0.
$$

The complete representatives satisfy $Q_{\mathrm{pert}}C_0=Q_{\mathrm{pert}}C_2=Q_{\mathrm{pert}}C_4=0$, and their products obey

$$
[C_0][C_2]=0,
\qquad
[C_2]^2=-3[C_0][C_4]\ne0,
\qquad
[C_2]^3=0.
$$

The [external product theorem](../sources/all-loop/theorems.md#b-00513) applies when quantum corrections preserve gauge invariance, as in this pure-SYM complex.

It gives associativity: changing the placement of parentheses leaves the product of classes unchanged.

Here $[C_2]$ is nilpotent of index three: its cube is the first positive power that vanishes.

The equalities hold over rational polynomials in $\hbar$, with the article sign $C_4=-C_{4,\mathrm{source}}$.

<a id="np-summary-034"></a>

### 2.3. Counting all leading two-ghost operators

A jet is a field derivative evaluated at the operator insertion point.

The integers $N_b,N_c$ count occurrences of $b,c$; the weight $(n_1,n_2)$ counts derivatives, with each $b$ contributing $(1,1)$.

The leading sector considered here has

$$
(N_b,N_c)=(k,2),
\qquad(n_1,n_2)=(k+1,k+1),
\qquad k\ge0.
$$

Its operators contain $k$ undifferentiated $b$ letters and two first derivatives of $c$.

Throughout, $\dim$ counts independent classes and $\mathbb Q$ denotes rational coefficients.

The subscript “scalar” selects expressions invariant under rotations of the two derivative directions.

Tree cohomology means solving $Q_0\mathcal O=0$ and identifying additions of $Q_0A$.

A finite-rank free module here means that every tree class is a unique polynomial combination of finitely many displayed trace expressions.

The polynomial coefficients are built from the gauge-invariant traces $\operatorname{tr}(b^2),\ldots,\operatorname{tr}(b^N)$.

This finite number of generators controls arbitrarily large $k$.

| Gauge group | Number of polynomial generators | Their $b$-degrees | $\dim H(Q_0)^{(k,2)}_{\mathrm{scalar}}$ |
|---|---:|---|---|
| $SU(2)$ | $1$ | $0$ | $1$ for even $k$; $0$ for odd $k$ |
| $SU(3)$ | $3$ | $0,1,2$ | $\lfloor k/2\rfloor+1$ |
| $SU(4)$ | $6$ | $0,1,2,2,3,4$ | $(\lfloor k/2\rfloor+1)(\lfloor k/2\rfloor+2)/2$ |

Here $\lfloor k/2\rfloor$ is the greatest integer less than or equal to $k/2$.

For example, the full $SU(2)$ module is

$$
H(Q_0)_{\mathrm{leading\ scalar}}
=\mathbb Q[\operatorname{tr}(b^2)]\,[C_0].
$$

The pages $E_1,E_2,E_3$ record classes before loops, after the first-loop test, and after the second-loop test; $E_\infty$ records survival through every loop.

The first test sends $[\mathcal O]$ to $[Q_1\mathcal O]$.

$$
E_1=H(Q_0),\qquad E_{r+1}=H(E_r,d_r),
\qquad d_1[\mathcal O]=[Q_1\mathcal O].
$$

The all-degree results in the specified leading scalar sector are

$$
\begin{array}{c|c}
SU(2)&\displaystyle
\dim E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}1,&k\text{ even},\\0,&k\text{ odd},\end{cases}\\[6pt]
SU(3)&\displaystyle
\dim E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}0,&k\equiv2\pmod3,\\1,&k\not\equiv2\pmod3.\end{cases}
\end{array}
$$

Cartan restriction sets all off-diagonal matrix components to zero.

A nonzero restricted $Q_1\mathcal O$ obstructs solving the first-correction equation $Q_0(\text{first correction})=-Q_1\mathcal O$.

The composition $\mathrm{Cartan}\circ Q_1$ means first apply $Q_1$, then set off-diagonal components to zero.

Its kernel consists of the classes whose restricted image vanishes.

For $SU(4)$, the complete classification of this test is

$$
\dim\ker(\mathrm{Cartan}\circ Q_1)\big|_{H(Q_0)^{(k,2)}_{\mathrm{scalar}}}
=\begin{cases}0,&k\equiv3\pmod4,\\1,&k\not\equiv3\pmod4.\end{cases}
$$

The vanishing sectors have $E_2=0$.

In the other sectors, $\dim E_2\le1$.

Full first corrections are established at $k=2,4,5$; the $k=6$ primitive is the remaining finite step in the all-degree $SU(4)$ upgrade.

<a id="np-summary-049"></a>

### 2.4. Complete representatives and higher-loop survival

For $SU(3)$, the source representative $a_1$ of $C_1$ is

$$
a_1=-\frac12\operatorname{tr}(b\partial_1c\partial_2c-b\partial_2c\partial_1c),
\qquad Q_0a_1=Q_1a_1=0,
\qquad [C_1]^2\ne0.
$$

The constructed $C_3,C_4$ representatives begin with the following explicit polynomials and contain exactly three nonzero loop coefficients:

$$
\begin{aligned}
[\hbar^0]C_3={}&-\frac{111}{112}\operatorname{tr}(b^2)a_1
+\frac{57}{56}\operatorname{tr}(b^3)C_0\\
&+\frac{15}{14}\bigl[
\operatorname{tr}(b^2\partial_2c)\operatorname{tr}(b\partial_1c)
-\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b\partial_2c)\bigr],\\[4pt]
[\hbar^0]C_4={}&\operatorname{tr}(b^3)a_1
-\frac32\operatorname{tr}(b^2\partial_1c\,b^2\partial_2c)
-\frac32\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c),\\[4pt]
C_3={}&[\hbar^0]C_3+\hbar[\hbar]C_3+\hbar^2[\hbar^2]C_3,\\
C_4={}&[\hbar^0]C_4+\hbar[\hbar]C_4+\hbar^2[\hbar^2]C_4.
\end{aligned}
$$

$$
Q_{\mathrm{pert}}C_3=Q_{\mathrm{pert}}C_4=0.
$$

“Complete” means that every loop coefficient of these last equations vanishes, including the higher-order equations acting on the three displayed coefficients.

The trace appendix supplies every coefficient.

Within the respective scalar spaces with $N_b+\text{loop order}=3,4$, all representatives of a fixed leading class differ by a full perturbative boundary.

A primitive is an operator whose $Q_{\mathrm{pert}}$ image equals a specified expression.

A partner below is an additional nonzero class with more ghost factors, constructed by combining primitives so that their $Q_{\mathrm{pert}}$ images cancel.

Its label records its leading letter counts and weight.

| Gauge group and class | Established result |
|---|---|
| $SU(2)$ first partner, $(N_b,N_c)=(1,5)$, weight $(5,5)$ | Complete representative; one-dimensional specified class space |
| $SU(2)$ partner class, $(3,5)$, weight $(7,7)$ | Complete representative; one-dimensional specified scalar space |
| $SU(2)$ leading $(5,5)$, weight $(9,9)$ | One-dimensional $E_\infty$; explicit coefficients through the third correction |
| $SU(3)$ first partner, $(0,7)$, weight $(5,5)$ | Complete representative; one-dimensional $E_\infty$; multiplication by $C_0$ gives zero; changes of primitive contribute zero through $C_0^2$ |
| $SU(3)$ leading $k=3,4$ | Complete $C_3,C_4$; one-dimensional respective homogeneous scalar spaces |
| $SU(3)$ leading $k=6$ | $\dim E_3^{(6,2)}(7,7)_{\mathrm{scalar}}=1$ |
| $SU(4)$ leading $k=2$ | Complete representative with one nonzero correction |
| $SU(4)$ leading $k=4,5$ | Exact first corrections |

For the $SU(2)$ weight-$(9,9)$ class, the fourth and fifth explicit coefficients remain to be solved.

For $SU(3)$ at $k=6$, later closure equations determine survival beyond $E_3$.

<a id="np-summary-059"></a>

### 2.5. A single formula for the diagonal representatives

Let $z$ be the auxiliary expansion variable, $\mathbf1$ the identity matrix, and $[z^{-1}]$ the coefficient of $z^{-1}$ in the expansion at infinity.

The polynomial

$$
[z^{-1}]\det(z\mathbf1-b)^{(k+1)/N}
$$

has degree $k+2$ in $b$.

Its Hessian is its matrix of second derivatives in the diagonal coordinates.

Inserting $\partial_1c$ and then $\partial_2c$ into its two slots gives the leading Cartan classes for $SU(3),SU(4)$.

Relative to the representatives displayed in the body, the exact factors are

$$
SU(3):\quad\frac{(k+1)^2}{3},
\qquad
SU(4):\quad-\frac{(k+1)^2}{4}.
$$

For every $N\ge2$, the polynomial and its Hessian vanish precisely at

$$
k\equiv N-1\pmod N.
$$

In the source terminology these are twisted periods: contour integrals of a fractional power of the characteristic polynomial, evaluated here by the displayed Laurent-coefficient extraction.

Feigin–Silantyev classify this polynomial family.

The operator calculation establishes its correspondence with the full Cartan test at $N=3,4$.

<a id="np-summary-067"></a>

### 2.6. Computation and remaining conjectures

The results use complete rational polynomial identities and exact rank bounds.

The product appendix supplies the one-loop displacement formula and all three two-loop kernels with zero input-jet order; the trace appendix and electronic certificates retain the representatives and primitives.

Stored loop operators use

$$
Q_r(\mathrm{physical})=(64\pi^2)^{-r}Q_r(\mathrm{stored}).
$$

For the proposed $SU(2)$ tower, the existing labels are $D_0=-C_0$ and $D_j=(2j-1)C_{2j}$, with $i,j$ indexing tower elements.

The source conjectures

$$
D_iD_j=D_0D_{i+j},\qquad i,j\ge1,
$$

with the $i=j=1$ case established by the exact product above.

The additional partner-pattern conjecture and computations modulo primes remain in the [exact-arithmetic section](08-exact-arithmetic.md#np-section-06), with their stated evidence.

<a id="np-summary-073"></a>

### 2.7. A conditional mechanism for removing perturbative classes

The operator $C_0$ is the gaugino bilinear, the product of two fermionic partners of the gauge field.

Write $L=C_0\star{}$ for multiplication by it: $L\mathcal O=C_0\star\mathcal O$.

The power $L^N$ applies this operation $N$ times.

A **multiplication homotopy** is an odd operator $h$ with $Q_{\rm pert}h+hQ_{\rm pert}=L^N$.

For a closed operator $\mathcal O$, this supplies the explicit primitive $h\mathcal O$ of $L^N\mathcal O$.

An **instanton** is a finite-action solution of the Euclidean field equations with nonzero integer gauge winding; this integer distinguishes configurations that cannot be continuously deformed into one another with the prescribed boundary behaviour.

In pure $SU(N)$ SYM, $\Lambda^{3N}$ denotes the one-instanton scale factor.

The proposed exact continuation uses the full supercharge $Q$ and a scalar coefficient $\nu$ in the source normalization:

$$
[Q,L]=0,
\qquad
\{Q,h\}=L^N-\nu\Lambda^{3N}\mathbf1,
\qquad \nu\ne0.
$$

Here $[Q,L]=QL-LQ$, $\{Q,h\}=Qh+hQ$, and $\mathbf1$ is the identity operator.

The displayed identity is the physical hypothesis of the following result.

For a perturbative class $a$ annihilated by $C_0$, choose representatives with

$$
Q_{\rm pert}a=0,
\qquad Q_{\rm pert}w=La.
$$

The **secondary class** combines the primitive $w$ of $La$ with the homotopy $h$.

In the following formula $h,L$ take their perturbative values, and $E_\infty^{\rm pert}$ is the space of classes surviving every perturbative correction:

$$
\langle C_0^{N-1},C_0,a\rangle
=[ha-L^{N-1}w]
\quad\text{in}\quad
E_\infty^{\rm pert}/C_0^{N-1}E_\infty^{\rm pert}
$$

The quotient identifies classes differing by a multiple of $C_0^{N-1}$, precisely the ambiguity caused by changing $w$.

The **instanton differential** $Q_{1\text{-inst}}$ is the first action on these classes induced by the one-instanton correction to $Q$.

Under the exact-product hypothesis,

$$
Q_{1\text{-inst}}\langle C_0^{N-1},C_0,a\rangle
+\langle C_0^{N-1},C_0,Q_{1\text{-inst}}a\rangle
=-\nu a
\quad\bmod C_0^{N-1}E_\infty^{\rm pert}.
$$

When the charges give $Q_{1\text{-inst}}a=0$, the right-hand side is the instanton image of the secondary class.

If $\nu a$ is nonzero in the quotient, this pairs two perturbative classes.

The verified perturbative inputs include $a=C_2$ for $SU(2)$ and $a=C_1$ for $SU(3)$.

[Construction and hypotheses](09-secondary-pairing.md#np-section-07).

<a id="np-summary-085"></a>

### 2.8. Finite homotopies and the formal loop parameter

The $SU(2)$ homotopy has been constructed through fourth order in the loop-counting parameter $\hbar$, on polynomials in undifferentiated $b$ with arbitrary ghost factors.

More generally, the solved coefficient ranges are

| Coefficient | Solved total input jet degree | Number of nonzero component coefficients |
|---|---:|---:|
| $[\hbar]h$ | $0\le m\le6$ | $84$ |
| $[\hbar^2]h$ | $0\le m\le4$ | $324$ |
| $[\hbar^3]h$ | $0\le m\le2$ | $163$ |
| $[\hbar^4]h$ | $m=0$ | $11$ |

Here the input jet degree counts derivatives on the input $b$ fields; a jet is a field derivative evaluated at the insertion point.

The coefficient of $\hbar^r$ is denoted $[\hbar^r]h$.

Extending this finite construction to every derivative and loop order is an open problem.

The first $SU(2)$ secondary class, with $a=C_2$, has the full representative

$$
\hbar M_1+\hbar^2M_0
$$

The polynomials $M_1,M_0$ contain respectively one and zero $b$ fields, and five ghosts each.

The coefficient ring $\mathbb Q[[\hbar]]$ consists of formal series with nonnegative powers of $\hbar$.

Thus $\hbar^{-1}$ lies outside this ring.

A nonzero class is **parameter torsion** when a positive power of $\hbar$ times that class is zero.

Adjoining $\hbar^{-1}$ makes such a class vanish.

The displayed secondary representative retains its explicit factors of $\hbar$.

For correction operators commuting with $Q_0$ when even or anticommuting with $Q_0$ when odd, the source proves a finite-support restriction.

On the full jet algebra of traceless matrices, assume the operator differentiates ghost jets, involves finitely many jet coordinates, and has coefficients polynomial in finitely many ghost jets.

Its positive differential-order part vanishes: its remaining action is multiplication by a $Q_0$-closed polynomial.

Differential order counts differentiations with respect to ghost coordinates.

The homotopy $h$ instead solves the displayed inhomogeneous equation.

<a id="np-summary-093"></a>

### 2.9. Exact one-instanton integrals and color transport

For $SU(2)$, let $x_0\in\mathbb R^4$ and $\rho>0$ be the instanton position and size.

Let $x_1,x_2$ be distinct insertion points and $r=x_1-x_2$.

Euclidean squares measure squared lengths.

The **separated-point integral** sums the two operator profiles over the position and size of the instanton:

$$
\int d^4x_0\int_0^\infty d\rho\,
\frac{\rho^9r^2}{((x_1-x_0)^2+\rho^2)^4((x_2-x_0)^2+\rho^2)^4}
=\frac{\pi^2}{90}.
$$

In the fixed physical-profile convention, the fermion integration supplies a factor $9$.

With $\mathcal N$ the remaining collective-measure normalization,

$$
\langle C_0(x_1)C_0(x_2)\rangle_{1\text{-inst}}
=\frac{\mathcal N\pi^2}{10}\Lambda^6.
$$

A **Wilson line** transports a color vector between two points using the gauge connection.

Here $U(0,r)$ transports the shifted ghost jet back to the origin along the straight segment of length $r>0$ in the fourth coordinate direction $e_4$.

The trace is taken in the adjoint color representation.

The positive normalized average $\mathbb E$ is

$$
\mathbb E[F]=\frac{90}{\pi^2}
\int_0^\infty d\rho\int_{\mathbb R^4}d^4y\,
\frac{\rho^9r^2F}
{(y^2+\rho^2)^4((y+r e_4)^2+\rho^2)^4}.
$$

Here $y$ is the instanton-center displacement and $F$ is the function being averaged.

Varying the transporter with the supercharge produces the **Wilson insertion**, whose integrated contribution satisfies

$$
\frac{\text{Wilson insertion}}{\langle C_0(0)C_0(r)\rangle}
=-\frac14+\frac4{27}\mathbb E[\operatorname{tr}U].
$$

$$
\frac{\text{Wilson insertion}}{\langle C_0(0)C_0(r)\rangle}
\ge\frac{301}{972}-\frac{1600\pi^2}{81081}
>\frac{5470949}{47675628}>0.
$$

The exact closed expression for $\mathbb E[\operatorname{tr}U]$ is given in the [straight-line calculation](10-instanton-measure.md#np-f863edf6248e0361-1).

<a id="np-summary-104"></a>

### 2.10. Regular response and collective-coordinate Ward balance

A **regular Green response** is the change of the instanton field obtained by solving its linearized field equation with a solution finite at the instanton center and the stated condition at infinity.

The response to the external curvature has three angular components.

Their degrees $l$ and color-angular eigenvalues $c$ are $(1,2),(3,5),(2,4)$.

Here $c$ labels the eigenvalue in the radial equation.

With $s=y^2$, the unknown $k(s)$ is the radial coefficient, $f(s)$ is its forcing, and $u(s)$ is the solution regular at the center with zero forcing; primes denote $d/ds$.

The [three-channel table](11-green-response.md#np-a6cb89ed04a182b4-1) supplies these inputs to

$$
4sk''+4(l+2)k'
+\left[\frac{4c}{s+\rho^2}-\frac{8s}{(s+\rho^2)^2}\right]k=f.
$$

$$
k(s)=-u(s)\int_s^\infty\frac{dt}{t^{l+2}u(t)^2}
\int_0^t\frac{dz}{4}z^{l+1}u(z)f(z).
$$

The resulting ghost jets acquire six independent Grassmann directions: four instanton fermion parameters and two external fermion parameters.

A Grassmann parameter anticommutes with the others and squares to zero.

The explicit nonzero six-by-six determinant is given in the [rank calculation](11-green-response.md#np-d147d35f585789c5-1).

**Collective coordinates** parametrize the instanton's position, size, color orientation, and fermion zero modes.

A fermion zero mode solves the homogeneous linearized fermion equation.

The supercharge induces an odd derivative on these coordinates.

**Ward balance** means that its integrated action equals the flux through the regulator boundaries, including the variation of the regulator.

Here $x_2$ is the source pure-ghost primitive with $Q_0x_2=C_0^2$ and $Q_rx_2=0$ for $r\ge1$; the expectation value uses its specified physical profile.

For the regulated size interval $\varepsilon\le\rho\le R$,

$$
\langle Qx_2\rangle_{\varepsilon\le\rho\le R}
=-\frac{\mathcal N\pi^2}{40}
-\left(-\frac{\mathcal N\pi^2}{40}\right)=0.
$$

At spatial cutoff $L$, the completed collective vector $Q_L$ obeys $Q_L^2=0$ on every retained coordinate and has zero divergence with respect to the collective measure.

The completion includes the even boundary values of the gauge potential and the orthogonal projection onto orientation zero modes.

[Finite-cutoff construction](12-finite-cutoff.md#np-section-10).

<a id="np-summary-112"></a>

### 2.11. The evaluated connection contribution and the remaining matching

The **connection insertion** varies the gauge connection inside the covariant derivatives of the partner polynomial.

The variation $\delta_R$ uses the part $R$ left after subtracting motion along the collective coordinates.

Its action includes

$$
\delta_Rc_{20}=[R_1,c_{10}],
\qquad
\delta_Rc_{30}=2[R_1,c_{20}]+[D_1R_1,c_{10}],
\qquad
\delta_R(D_\alpha b)=[R_\alpha,b].
$$

Here $c_{rt}=\partial_1^r\partial_2^tc$ labels the holomorphic jet, while $\delta_R$ acts on its covariant physical representative.

The source dictionary symmetrizes the derivative slots together with the gaugino slot with total weight one.

The derivative $D_\alpha$ transforms in the same gauge representation as its input, and brackets are color-matrix commutators.

Let $b_H$ be the external curvature parameter, $\xi_1,\xi_2,\eta^{\dot1},\eta^{\dot2}$ the instanton fermion parameters, and $e_1,e_2$ the external fermion parameters.

A subscript on square brackets extracts the coefficient of the indicated ordered monomial.

The projection coefficient $P_L$ transfers the orientation component of the response to collective motion at spatial radius $L$: the remaining response has zero finite-ball inner product with each orientation zero mode.

For the Green solution with upper endpoint at infinity and fixed finite $P_L$,

$$
\int_{\mathbb R^4}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}
=-\pi^2\left[
\frac{3364P_L}{10395\rho^4}
+\frac{724901g^2}{27941760\rho^2}\right].
$$

Here $g$ is the physical source coupling.

Substituting the finite-ball projection gives the logarithmic and finite coefficients

$$
\lim_{L\to\infty}
\frac{\rho^2}{\pi^2g^2\log(L^2/\rho^2)}
\int_{|y|\le L}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}
=-\frac{841}{41580}
$$

$$
\begin{aligned}
\lim_{L\to\infty}\bigg[&
\int_{|y|\le L}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}\\
&+\frac{841\pi^2g^2}{41580\rho^2}\log(L^2/\rho^2)\bigg]
=\frac{197929\pi^2g^2}{9313920\rho^2}.
\end{aligned}
$$

These evaluate the connection contribution to $M_1$.

The full secondary insertion is $\hbar M_1+\hbar^2M_0$.

**Local Ward matching** compares the supercharge action on renormalized local operator products with its action on the instanton collective coordinates, in a common subtraction prescription.

The prescription removes coincident-point singularities by local counterterms.

Contact terms are contributions supported where insertion points coincide.

For the renormalized composite $[C_0^2]_{\rm ren}$,

$$
C_0\star C_0=[C_0^2]_{\rm ren}+\nu_{\rm prod}\Lambda^6,
\qquad
Qx_2=[C_0^2]_{\rm ren}+\nu_Q\Lambda^6.
$$

The two scalar coefficients shift together when the local composite is redefined; their difference $\nu_{\rm prod}-\nu_Q$ is invariant.

The full physical coefficient of the secondary pairing remains open.

Its matching includes the normal gaugino variation, the second-order response of $M_0$, the action and measure, antifield terms encoding equations of motion, and the coincident composite/contact contributions.

[Remaining equation](13-connection-matching.md#np-b0509e2b7b2fc17a-1).

<a id="np-summary-124"></a>

### 2.12. Massive matter and the vacuum algebra

For $SU(2)$, a massive adjoint has scalar invariant $u=\gamma^i\gamma_i/2$ and superpotential $W=mu$.

One massive fundamental flavor consists of two color-doublet scalar fields $q,t$, with meson $M=\epsilon_{ab}q^at^b$ and $W=mM$.

Here $i$ labels the three adjoint components, with lowering by the source invariant adjoint metric; $m\ne0$ is the mass, $\gamma$ is the even matter field, $\beta$ its odd dual, and $\epsilon_{ab}$ the invariant alternating color pairing.

The superpotential supplies the scalar interaction and its derivatives give the scalar field equations.

The **Konishi identity** is the supercharge variation of the matter bilinear, including its one-loop gauge contribution:

$$
Q_{\rm pert}(\beta_i\gamma^i)=2mu+\hbar C_0
$$

$$
Q_{\rm pert}(\beta_q q+\beta_t t)=2mM+\frac{\hbar}{2}C_0.
$$

$$
[u]=-\frac{\hbar}{2m}[C_0],
\qquad [M]=-\frac{\hbar}{4m}[C_0].
$$

The low-energy scale $\Lambda_L$ after integrating out massive matter satisfies

$$
\Lambda_L^6=m^2\Lambda^4
\quad\text{(massive adjoint)},
\qquad
\Lambda_L^6=m\Lambda^5
\quad\text{(one massive flavor)}.
$$

Assuming the exact low-energy relation $[C_0]^2=\nu\Lambda_L^6$ gives

$$
[u]^2=\frac{\hbar^2\nu}{4m^2}\Lambda_L^6
=\frac{\hbar^2\nu}{4}\Lambda^4
$$

$$
[M]^2=\frac{\hbar^2\nu}{16m^2}\Lambda_L^6
=\frac{\hbar^2\nu}{16}\frac{\Lambda^5}{m}.
$$

At $m=0$, the adjoint bilinear instead satisfies $Q_{\rm pert}(\beta_i\gamma^i)=\hbar C_0$.

Thus multiplication of $[C_0]$ by $\hbar$ is zero in formal cohomology; after adjoining $\hbar^{-1}$, $C_0$ has primitive $\hbar^{-1}\beta_i\gamma^i$.

At fixed complex $\nu\Lambda^{3N}\ne0$, the **gaugino vacuum algebra** is the polynomial algebra in $C_0$ modulo its exact relation.

It separates into $N$ vacua.

Their distinct gaugino values are $s_k=\nu^{1/N}\Lambda^3e^{2\pi ik/N}$, with $k=0,\ldots,N-1$, and $e_k$ is the polynomial that evaluates to one in vacuum $k$ and zero in every other vacuum:

$$
e_k=\frac1N\sum_{\substack{j=0,\ldots,N-1\\\mathrm{finite\ geometric\ series}}}
\left(\frac{C_0}{s_k}\right)^j.
$$

$$
\mathbb C[C_0]/(C_0^N-\nu\Lambda^{3N})
=\bigoplus_{\substack{k=0,\ldots,N-1\\\mathrm{vacua}}}\mathbb C e_k,
\qquad
C_0e_k=s_ke_k,
\qquad e_ke_l=\delta_{kl}e_k.
$$

Here $\delta_{kl}$ equals one for $k=l$ and zero otherwise.

This algebraic decomposition concerns the gaugino subalgebra at fixed nonzero scale.

The relation to the full local-operator complex uses the preceding conditional pairing and physical matching.

[Matter calculation](14-massive-matter.md#np-section-12).

<a id="np-summary-139"></a>

### 2.13. SQCD: a determinant removed by an instanton

In four-dimensional supersymmetric QCD (SQCD), the gauge fields couple to matter multiplets, each containing a scalar and a fermion.

Here $N_c$ counts colors in the gauge group $SU(N_c)$, and $N_f$ counts pairs of fundamental and antifundamental multiplets, whose color scalars are a column and a dual row.

The meson $M=Q\widetilde Q$ is their gauge-invariant scalar product; here $Q,\widetilde Q$ denote quark scalars, while the supercharge acts on operators.

Baryons $B,\widetilde B$ are gauge-invariant products of $N_c$ quarks or antiquarks.

Seiberg duality describes the same low-energy physics using a magnetic gauge group $SU(N_f-N_c)$, quarks $q,\widetilde q$, and an elementary meson $M$.

The superpotential $W$ specifies the holomorphic interactions, depending on complex scalar fields rather than their conjugates.

Its scalar equations $\partial W=0$, called F-equations, impose relations among the effective scalar operator classes.

The matching parameter $\mu$ fixes the meson normalization; $\Lambda,\widetilde\Lambda$ are the electric and magnetic strong-interaction scales.

For electric $SU(2)$ with four flavors, give the fourth flavor mass $m\ne0$.

Its magnetic scalar expectation values satisfy $q_4\widetilde q^{\tilde4}=-\mu m$, making every magnetic gauge boson massive.

This region of field space is the **magnetic Higgs patch**.

Hats denote the three retained light flavors.

The known one-instanton superpotential of [Seiberg](https://arxiv.org/abs/hep-th/9411149) is

$$
\Lambda^2\widetilde\Lambda^2=\mu^4,
\qquad
W_{\mathrm{eff}}
=\frac1\mu\widehat M^i{}_{\tilde\jmath}\widehat q_i
\widehat{\widetilde q}^{\tilde\jmath}
-\frac{\widetilde\Lambda^2}{\mu^4m}\det\widehat M.
$$

Here $i,\tilde\jmath$ label the three light flavors.

At fixed magnetic coordinates, our explicit relation is

$$
\widehat M^i{}_{\tilde\jmath}
\frac{\partial W_{\mathrm{eff}}}{\partial\widehat M^i{}_{\tilde\jmath}}
-\widehat q_i\frac{\partial W_{\mathrm{eff}}}{\partial\widehat q_i}
=-\frac{3\widetilde\Lambda^2}{\mu^4m}\det\widehat M.
$$

Thus the determinant is a nonzero perturbative class and becomes zero in the effective scalar ring at fixed $\widetilde\Lambda^2\ne0$.

The explicit odd operator whose $Q$-variation equals the conjugate determinant is given in the [determinant calculation](16-determinant-lift.md#np-a3c60c764266802e-1).

Allowing the inverse instanton weight defines a map out of the scalar ring. Its kernel consists of the scalar classes that disappear. For any size of light meson in this magnetic presentation, the full kernel and the coefficients that annihilate it are

$$
\ker(\text{inversion})=(\det\widehat M)\,\mathbb C[\widehat M],
\qquad
\operatorname{Ann}(\det\widehat M)=
(\widetilde\Lambda^{4-N_c},\widehat q,\widehat{\widetilde q}).
$$

Here $\operatorname{Ann}$ denotes all ring elements whose product with the given class is zero. The two quark entries in the ideal stand for all their components. Every class in this kernel is removed at the first power of the instanton weight. [Full scalar calculation](16-determinant-lift.md#determinant-full-torsion)

<a id="np-summary-147"></a>

### 2.14. The microscopic operator that can pair with the determinant

We construct an incoming operator in the full specified perturbative complex.

The matter calculation writes its full perturbative supercharge as $\mathbf Q_{\mathrm{pert}}$.

Here $\beta^M,\beta^q$ are the odd partners recording the meson and quark F-equations.

A bar denotes complex conjugation.

The expression $A$ contracts three light rows and all four meson columns:

$$
A=\overline M^{\bar i}{}_{\bar{\tilde\jmath}}
\beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}
-\overline q_{\bar i}\beta^{q\,\bar i},
\qquad h=q_4\widetilde q^{\tilde4}.
$$

Let $v,t$ be the heavy conjugate quark column and dual row, so $\bar h=tv$.

The ghost $\mathsf c_\gamma$ encodes infinitesimal gauge transformations of these scalars; $g$ is the source gauge coupling.

The derivatives $\partial_{+\dot\alpha}$ are the retained holomorphic spacetime derivatives; repeated dotted indices use the source spinor pairing.

With $a,b$ abstract color slots, the correction has polynomial numerator

$$
\bar hK=
t_a(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b
(\partial_{+\dot\alpha}v)^b
-(\partial_{+\dot\alpha}t)_a
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b v^b.
$$

The exact triangle coefficient in this convention is $g^2/(64\pi^2)$.

At every perturbative order,

$$
\mathbf Q_{\mathrm{pert}}
\left(\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right)=0,
\qquad
\left[\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right]\ne0
\quad(m\ne0).
$$

Matching this microscopic operator to the effective instanton interaction would give the one-instanton image

$$
-\frac{3\overline{\widetilde\Lambda}^{\,2}}{\bar\mu^3}
\det\overline{\widehat M}.
$$

That microscopic equality remains conditional on matching the local composite operators and contact terms, meaning the extra terms produced when insertion points coincide.

At $m=0$, let $\beta^{M^4{}_{\tilde4}}$ be the heavy-meson odd partner and let $\operatorname{Tr}_{\mathrm{fund}}$ trace color matrices in the fundamental representation.

The explicit closed operator

$$
C=\bar\mu\beta^{M^4{}_{\tilde4}}
\operatorname{Tr}_{\mathrm{fund}}
[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)]
+\frac{i}{g}\bar hK
$$

satisfies

$$
[C]\ne0\quad\text{over }\mathbb C[[\hbar]],
\qquad
\left[\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right]_{m=0}
=\frac{3\hbar g^2}{64\pi^2}[C].
$$

The formal series ring $\mathbb C[[\hbar]]$ admits nonnegative powers of $\hbar$.

The incoming class survives on the zero-mass fiber. Independent choices of the two meson flavor labels give sixteen independent copies over the formal series ring. Each component obeys

$$
\hbar^k[C]\ne0\qquad(k\ge0).
$$

For microscopic matching, let $D$ be the coefficient of one instanton and zero loops in a prescription admitting nonnegative loop powers. Evaluate $D A$ with all fermions and derivatives zero, arbitrary light mesons and heavy quark contraction $-\bar\mu\bar m$. The proposed determinant image at nonzero mass requires

$$
\operatorname{ev}_{\bar m}(D A)
=\frac{3}{\bar\mu^4\bar m}\det\overline{\widehat M}.
$$

This fixes the residue any such prescription must reproduce. [All-loop obstruction](16-determinant-lift.md#massless-all-powers) · [Mass-pole test](16-determinant-lift.md#mass-pole-ward)

There is a further distinction when the mass and loop parameter vary together. Let $R=\mathbb C[[\bar m,\hbar]]$, let $H_R$ be the microscopic perturbative cohomology over this ring, and let $Z$ be the quadratic ghost trace appearing in $C$. Exactly the series with zero constant term annihilate $[Z]$. The polynomial incoming class gives a nonzero element in the kernel of

$$
H_R/(\bar m,\hbar)H_R\ \longrightarrow\
H\!\left(\left.\mathbf Q_{\rm pert}\right|_{\bar m=\hbar=0}\right),
\qquad \operatorname{Ann}_R[Z]=(\bar m,\hbar).
$$

Thus setting both parameters to zero after taking cohomology retains a class that becomes a tree variation when they are set to zero in the differential. [Joint specialization and its explicit second differential](16-determinant-lift.md#mass-loop-base-change)

<a id="np-summary-162"></a>

### 2.15. Vacuum constraints and nonpolynomial superpotentials

The exact SQCD dynamics of [Intriligator–Seiberg](https://arxiv.org/abs/hep-th/9509066) removes magnetic vacua with excessive meson rank.

For an invertible electric mass matrix $m$, the gaugino condensate $S$ instead labels $N_c$ isolated massive vacua:

$$
\operatorname{rank}\langle M\rangle\le N_c
\quad\text{in the massless theory},
$$

$$
S_k=e^{2\pi ik/N_c}
(\Lambda^{3N_c-N_f}\det m)^{1/N_c},
\qquad
\langle M\rangle_k=S_km^{-1},
\qquad W_{\mathrm{vac},k}=N_cS_k,
\quad k=1,\ldots,N_c.
$$

Fractional powers retain distinct condensate values.

For $N_f<N_c$, the polynomial relation and its superpotential value are

$$
S^{N_c-N_f}\det M=\Lambda^{3N_c-N_f},
\qquad W_{\mathrm{eff}}=(N_c-N_f)S.
$$

For $SU(2)$, eliminating a flavor of mass $m_2$ from the two-flavor constrained description gives the Affleck–Dine–Seiberg (ADS) superpotential.

Here $M$ is the remaining meson and $\Lambda_L$ the daughter scale:

$$
W_{\mathrm{eff}}=\frac{\Lambda_L^5}{M},
\qquad \Lambda_L^5=m_2\Lambda^4,
\qquad M\ne0.
$$

It has no finite stationary vacuum.

Adding a second nonzero mass gives

$$
W_{\mathrm{eff}}=\frac{\Lambda_L^5}{M}+m_1M,
\qquad M^2=\frac{\Lambda_L^5}{m_1},
\qquad W_{\mathrm{vac},\pm}=\pm2\sqrt{m_1m_2\Lambda^4}.
$$

The pole records the field-elimination domain $M\ne0$.

Retaining the condensate $S$ gives a logarithmic superpotential; eliminating it selects a root of the displayed polynomial relation.

[The full effective calculation](17-nonpolynomial-superpotentials.md#np-section-15) keeps these domains and all roots explicit.

For massless $SU(2)$ with two flavors, the four-by-four antisymmetric meson $V$ obeys

$$
\operatorname{Pf}V=V^{12}V^{34}-V^{13}V^{24}+V^{14}V^{23}=\Lambda^4.
$$

This fixes the ordering convention for the Pfaffian polynomial; its square is the determinant.

The resulting vacuum space is smooth and has five complex dimensions.

With three flavors, $V$ is six by six and

$$
W=-\frac{\operatorname{Pf}V}{\Lambda^3},
\qquad \operatorname{rank}V\le2.
$$

All fifteen meson multiplets remain massless at $V=0$.

These known interactions implement the electric composite relations and retain a continuous vacuum space.

The symplectic and orthogonal examples retain their source normalizations.

The compact group $SP(N_c)$ preserves a nondegenerate alternating color form; its fundamental has complex dimension $2N_c$, and $2N_f$ matter multiplets form an antisymmetric meson matrix $M$.

The three adjacent [Intriligator–Pouliot](https://arxiv.org/abs/hep-th/9505006) results for $SP(2)$ are

$$
W_{N_f=4}=-\frac{\operatorname{Pf}M}{2\Lambda_{2,4}^5},
\qquad
\operatorname{Pf}M\big|_{N_f=3}=2\Lambda_{2,3}^6,
\qquad
W_{N_f=2}=\frac{2\Lambda_{2,2}^7}{\operatorname{Pf}M}.
$$

The subscripts on $\Lambda$ specify the color and flavor counts.

The group $SO(N_c)$ consists of determinant-one transformations preserving a positive symmetric form on the real color space.

Its matter scalars transform as complex color vectors, and their meson is symmetric.

For $N_f=N_c-1$, the known instanton term and condensate give

$$
W_{\mathrm{inst}}=-\frac{\det M}{64\Lambda_{N_c,N_c-1}^{2N_c-5}},
\qquad
W_{\mathrm{cond}}=\epsilon\frac{\det M}{64\Lambda_{N_c,N_c-1}^{2N_c-5}},
\qquad \epsilon\in\{1,-1\}.
$$

The condensate choice $\epsilon=1$ retains a continuous family of vacua by exact cancellation.

The [orthogonal analysis](18-sp-and-so.md#np-7f93d3422761f50f-1) includes the other condensate choices, Coulomb singularities, and the distinct allowed electric and magnetic line charges.

The confined $SP(1)=SU(2)$ theory with three flavor pairs has another magnetic description. Hold the fifteen entries of its antisymmetric meson $V$ fixed while turning off that description's instanton weight. Its full scalar kernel is generated by the fifteen Pfaffians of the four-by-four matrices obtained by deleting two rows and columns of $V$. At zero weight the two coordinate choices give spaces of dimensions $15$ and $9$. At nonzero weight the superpotentials agree after rescaling the confined baryons to the light magnetic quarks. In the holomorphic Pfaffian effective model, one variation replaces one odd partner by two meson letters. Among targets with two scalar letters and $p-1$ odd partners and no derivatives, the number of independent images is

$$
\binom{15}{p}\quad(1\le p\le15),\qquad \binom{15}{2}=105 .
$$

[Coordinate comparison](18-sp-and-so.md#pfaffian-full-torsion) · [Fermionic images](18-sp-and-so.md#pfaffian-jet-images)

<a id="np-summary-183"></a>

### 2.16. Fermionic interactions and monopole-class removal

A **higher F-term** is a supersymmetric interaction whose integral over two anticommuting supersymmetry coordinates contains additional fermion factors.

An ordinary scalar superpotential is the case with zero such pairs.

In the smooth Higgs-branch effective descriptions of the [Beasley–Witten](https://arxiv.org/abs/hep-th/0409149) $SU(2)$ examples, $p=N_f-N_c+1$ counts those pairs.

Here scalar expectation values break the gauge group and the heavy fields have been eliminated.

Integrating out successive massive flavors gives

$$
N_f:4\longrightarrow3\longrightarrow2\longrightarrow1,
\qquad
p:3\longrightarrow2\longrightarrow1\longrightarrow0.
$$

The corresponding higher F-term effects are an eight-fermion interaction, a six-fermion interaction, a deformation of the meson constraint, and the ADS superpotential.

At $N_f=3$, the six-fermion term arises on a smooth Higgs patch of the meson theory with $W=-\operatorname{Pf}V/\Lambda^3$; both descriptions therefore belong to the same theory.

The eight-fermion coefficient contains the instanton factor $\Lambda^2$; its absolute numerical normalization remains unspecified in the source.

Thus supersymmetric instanton data include interactions beyond an ordinary function of scalar fields.

After one spatial direction is compactified to a circle, a three-dimensional **monopole operator** is inserted by prescribing magnetic flux through a small sphere around its insertion point.

For $SU(N_c)$ with $N_f=N_c\ge2$, let $Y$ be the monopole operator and $\eta=\Lambda^{2N_c}$ the finite-circle instanton parameter.

The known [circle superpotential](https://arxiv.org/abs/1305.3924) is

$$
W=Y(B\widetilde B-\det M+\eta).
$$

Our explicit F-equation relation gives the following complete scalar result:

$$
\eta=0:\quad [Y^k]\ne0\quad(k\ge1),
\qquad
\eta\ne0:\quad [Y]=0,\quad \det M-B\widetilde B=\eta.
$$

Here brackets are taken in the stated effective scalar ring.

At zero circle parameter, the starting three-dimensional theory already includes its other nonperturbative interactions.

**Parameter torsion** means that multiplication by a parameter makes a nonzero class vanish.

With formal coefficients in $\mathbb C[[\bar\eta]]$, the explicit conjugate relation is

$$
\bar\eta[\bar Y]=0.
$$

The inverse $1/\bar\eta$ is absent from that formal ring.

In a model with fixed nonzero $\bar\eta$, it is a scalar coefficient, and the [explicit primitive](19-higher-f-terms-and-monopoles.md#np-66ee4aae6c4219a5-1) gives $[\bar Y]=0$.

This is also the distinction behind the determinant and ADS examples.

The entire circle-parameter kernel is generated by $Y$. In each positive power $Y^j$ the allowed coefficients are meson polynomials modulo all $(N_c-1)$-rowed minors. For a general multiplier interaction $W=Y(F+\eta)$, where $F$ is a polynomial of the other scalar fields, the equation $\eta=-F$ acts on these coefficients modulo the derivatives of $F$. Its powers determine the exact order of each removal. In the local polynomial model $F=x^5+y^5+x^3y^3$ near $x=y=0$ one obtains, for every $j\ge1$,

$$
\eta[Y^j]=\tfrac15[Y^j x^3y^3]\ne0,\qquad
\eta^2[Y^j]=0 .
$$

The parameter enters the interaction linearly. Its second power first kills the displayed class. [Proof and explicit primitive](19-higher-f-terms-and-monopoles.md#multiplier-pages)

For a finite complex of free modules over $\mathbb C[[\eta]]$, express the incoming variations in a basis of the kernel of the outgoing variation. Row and column operations invertible over this ring put the matrix in diagonal form. An entry $\eta^r$ times an invertible series pairs one source and one target at order $r$. The minimum exponents $v_j$ among nonzero $j$-rowed minors determine every such order, in increasing order, by $r_j=v_j-v_{j-1}$ with $v_0=0$. [Finite-sector calculation](23-operator-data.md#complete-tests)

Two known two-dimensional benchmarks show different outcomes.

In the minimal $(0,2)$ model with target $\mathbb{CP}^1$, [Yagi](https://arxiv.org/abs/1001.0118) finds an odd operator $\theta$ whose anticommutator with the supercharge is a nonzero constant:

$$
\{Q,\theta\}=\Lambda^2e^{-t_0}O_1,
\qquad O_1\text{ is a nonzero constant operator}.
$$

Here $e^{-t_0}$ is the instanton weight.

At fixed nonzero weight, every closed local operator is exact.

In the [Hori–Vafa](https://arxiv.org/abs/hep-th/0002222) $(2,2)$ model with target $\mathbb{CP}^{N-1}$, the scalar $\sigma$ in the gauge field-strength multiplet instead obeys

$$
\sigma^N=\Lambda^N,
\qquad
\sigma=\Lambda e^{2\pi i\ell/N},
\quad\ell=0,\ldots,N-1.
$$

The labels $(0,2)$ and $(2,2)$ count left- and right-moving supersymmetries; $\mathbb{CP}^{N-1}$ is the space of complex lines in $\mathbb C^N$.

The second model retains $N$ vacua.

The three-dimensional [Aharony examples](19-higher-f-terms-and-monopoles.md#np-aa0e649d863ebf40-1) likewise remove specified magnetic monopole classes while retaining the stated meson and electric-monopole branches.

<a id="np-summary-200"></a>

### 2.17. A scalar invisible to vacuum values that survives every perturbative order

Take $SU(5)$ SQCD with four flavors and a traceless adjoint scalar $X$, meaning a matrix transformed by conjugation under the gauge group.

With nonzero coupling $s_0$ and vanishing adjoint mass $s_1$,

$$
W=\frac{s_0}{3}\operatorname{Tr}X^3,
\qquad
F(X)=X^2-\frac{\operatorname{Tr}X^2}{5}I.
$$

Here $I$ is the identity matrix, and $(F)$ denotes all polynomial combinations of the entries of $F(X)$.

We determine the complete gauge-invariant polynomial quotient of the adjoint F-equations:

$$
\left(\mathbb C[\mathfrak{sl}_5]/(F)\right)^{SU(5)}
=\frac{\mathbb C[\operatorname{Tr}X^2]}
{((\operatorname{Tr}X^2)^3)},
\qquad [(\operatorname{Tr}X^2)^2]\ne0.
$$

Here $\mathfrak{sl}_5$ is the space of traceless complex five-by-five matrices, and $\mathbb C[\mathfrak{sl}_5]$ denotes their polynomial functions.

The quotient has the three independent classes $1,[\operatorname{Tr}X^2],[(\operatorname{Tr}X^2)^2]$.

The last class is **nilpotent**: a positive power vanishes although the class itself is nonzero.

Every ordinary vacuum gives $\operatorname{Tr}X^2=0$.

For the conjugate scalar letter $\gamma$, our full perturbative calculation gives

$$
\mathbf Q_{\mathrm{pert}}(\operatorname{Tr}\gamma^2)^2=0,
\qquad
\hbar^r[(\operatorname{Tr}\gamma^2)^2]\ne0
\quad(r\ge0)
\quad\text{over }\mathbb C[[\hbar]].
$$

The [Kutasov–Schwimmer–Seiberg](https://arxiv.org/abs/hep-th/9510222) magnetic trace presentation predicts the stronger relation $(\operatorname{Tr}X^2)^2=0$ at the mesonic origin, where the meson expectation values vanish.

Extending it to the full quantum operator ring requires the [specified operator matching](20-cubic-adjoint-su5.md#np-1073436cbe3d7b68-1), including the dynamical mesons.

That extension remains conditional.

<a id="np-summary-209"></a>

### 2.18. Extended supersymmetry: corrected geometry and protected scalars

In four-dimensional $\mathcal N=2$ gauge theory, eight supercharges constrain the low-energy couplings through a holomorphic function $F$, the **prepotential**.

The electric and magnetic periods $a,a_D$ are contour integrals on the Seiberg–Witten curve, with $a_D=\partial_aF$.

They determine the coupling and the exact mass bound for charged supersymmetric particles.

In the pure $SU(2)$ source normalization with strong scale one, $F_k$ is the coefficient of $a^{2-4k}$, and the checked instanton coefficients are

$$
F_1=\frac1{8\pi i},
\qquad F_2=\frac5{256\pi i}.
$$

These are known [Seiberg–Witten](https://arxiv.org/abs/hep-th/9407087) data, reproduced through the period integral.

Adding an adjoint mass $m$ selects two vacua from the Coulomb branch, the continuous vacuum family with massless Abelian gauge fields.

With the source scale restored,

$$
u_{\mathrm{vac}}=\pm\Lambda^2,
\qquad W_{\mathrm{vac},\pm}=\pm m\Lambda^2,
\qquad M\widetilde M=i\sqrt2m\Lambda
\quad\text{in the monopole frame}.
$$

Here $M,\widetilde M$ are the light monopole fields, and $u$ is the Coulomb coordinate.

The second vacuum uses a dyon, a particle carrying both electric and magnetic charge.

The [SQCD curve calculations](21-seiberg-witten.md#np-section-19) also retain the exact flavor-dependent singularities, coupling-coordinate corrections, and the Argyres–Seiberg description of $SU(3)$ with six matter pairs by an $SU(2)$ gauge theory coupled to an interacting sector with $E_6$ global symmetry.

A global symmetry transforms physical states while preserving their dynamics; gauge transformations identify equivalent descriptions of a state.

The single remaining Coulomb coordinate $\widetilde v$ gives

$$
\widetilde y^2=\widetilde x^3-\widetilde v^4.
$$

The variables $\widetilde x,\widetilde y$ are coordinates on that auxiliary curve.

In the holomorphic-topological twist of pure $\mathcal N=2$ theory, local operators have an integer grading called ghost number. Each basic odd field has ghost number one. Write $\partial$ for translation in the holomorphic direction. For a simple gauge algebra, all tree classes of ghost number two form one derivative tower, with one class at each integer spin $s\ge2$. For $SU(2)$ write $X=-[\delta_{ij} b^i\partial c^j]$ for its lowest class and $Y=[\epsilon_{ijk}b^ib^jb^k]$ for the invariant three-fermion class. Here $b,\partial c$ are the two lowest odd fields, and $\delta_{ij},\epsilon_{ijk}$ are the invariant color tensors in Khan's basis. A bulk correction preserving translation and spin and raising ghost number by one satisfies

$$
d(\partial^k X)=d(\partial^k Y)=0\qquad(k\ge0).
$$

Every derivative in these two towers also survives every incoming differential under the same formal deformation hypotheses. A proposed map $\partial X\mapsto Y$ fails this test. The tower calculation is independent of the conjecture that $X,Y$ generate the complete perturbative algebra. [Bulk Ward calculation](21-seiberg-witten.md#bulk-translation) · [Khan's proposed algebra](https://arxiv.org/abs/2604.03500)

In $\mathcal N=4$ $SU(N)$ SYM, sixteen supercharges accompany three complex adjoint scalars $X,Y,Z$.

Their scalar relations are

$$
[X,Y]=[Y,Z]=[Z,X]=0.
$$

The [invariant-theory input](https://arxiv.org/abs/math/0602660) identifies every nonzero gauge-invariant scalar class with a polynomial that has a nonzero value on some simultaneous diagonal vacuum.

Those vacua survive in the exact theory.

Consequently every such nonzero scalar class remains nonzero under the exact supercharge.

Total holomorphic translations preserve this protection for the half-BPS primaries made of one scalar $Z$. At degree $p\ge2$ in $Z$ and total derivative order $k\ge0$, the dimension of their translated subspace is

$$
(k+1)[t^p]\prod_{\substack{j=2,\ldots,N\\\mathrm{invariant\ degrees}}}(1-t^j)^{-1}.
$$

The bracket extracts the coefficient of $t^p$. For $SU(2)$ the single primary $\operatorname{Tr}Z^2$ therefore gives $k+1$ independent spinful states at $k>0$. The exact two-point functions detect every linear combination in that subspace. [Descendant calculation](22-n4-sym.md#n4-descendant-detector)

A **long multiplet** is a supersymmetry family whose scaling dimension is free to change with the coupling.

Let $\mathcal O$ be the Konishi scalar.

Its anomalous dimension $\gamma_K$ changes the scaling power from the classical value two.

For a positive dilation $t$ of separated positions,

$$
\langle\mathcal O(tx)\mathcal O(0)\rangle=t^{-2(2+\gamma_K)}\langle\mathcal O(x)\mathcal O(0)\rangle,
\qquad t>0,\quad x\ne0.
$$

Use the coupling $g$ and vacuum angle $\theta$ of [Alday–Korchemsky](https://arxiv.org/abs/1605.06346), with

$$
\tau=\frac\theta{2\pi}+\frac{4\pi i}{g^2}.
$$

For positive real argument $z$, the Euler gamma function and the group coefficient $\kappa_N$ are

$$
\Gamma(z)=\int_0^\infty t^{z-1}e^{-t}\,dt,
\qquad z>0,
\qquad
\kappa_N=\frac{2\Gamma(N-\tfrac12)}{\sqrt\pi\Gamma(N-1)},
\qquad \kappa_2=1.
$$

The leading one-instanton and anti-instanton correction is

$$
\left.\gamma_K\right|_{\mathrm{leading\ one\ instanton\ and\ anti\ instanton}}
=-\frac{27\kappa_N}{5(N^2-1)}
\left(\frac{g^2}{4\pi^2}\right)^2
(e^{2\pi i\tau}+e^{-2\pi i\bar\tau}).
$$

These are genuine instanton corrections.

The known removal of classes in the $SO(7)$–$USp(6)$ comparison instead occurs at first perturbative quantum order, as calculated in the [Konishi-pairing analysis](22-n4-sym.md#np-87f770bad4f64cbb-1).

For a finite-dimensional space of perturbative survivors with positive inner product, our positivity calculation excludes a leading energy shift made entirely of nonzero Fourier modes in $\theta$.

An instanton correction to the supercharge remains compatible with positivity: if its first coefficient is $d$, then

$$
\left.\{Q_{\mathrm{eff}},Q_{\mathrm{eff}}^\dagger\}\right|_{\mathrm{lead}}
=e^{-16\pi^2/g^2}\{d,d^\dagger\}\ge0.
$$

The dagger denotes the adjoint in that inner product.

A fully perturbatively surviving $\mathcal N=4$ operator with a calculated nonzero microscopic instanton image remains open beyond the pure-scalar sector.
