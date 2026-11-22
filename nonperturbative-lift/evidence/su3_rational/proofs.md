# SU(3) first partner: rational chain calculation

The normalization is the original `checks/su3` component-engine normalization: $C_0=\operatorname{tr}(\partial_1c\,\partial_2c)$, with $c$ in the Chevalley basis and $b$ in the trace-dual basis. The coefficients below are the exact rational lifts of those checks. The original article and checks remain unchanged.

The selected source representative is exactly

$$
a_1=-\frac12\operatorname{tr}\bigl(b\,\partial_1c\,\partial_2c-b\,\partial_2c\,\partial_1c\bigr).
$$

The degree labels $N_b,N_c$ count the $b,c$ letters. The weights $(n_1,n_2)$ count holomorphic derivatives, with each $b$ contributing $(1,1)$. Every $c$ is differentiated. All displayed operators are gauge and rotation invariant.

**Proposition 1.** The rational lifts of Stage A satisfy

$$
\begin{aligned}
Q_0x_3&=C_0^3,& Q_0a_1&=0,&Q_1a_1&=0,\\
Q_0w_1&=C_0a_1,&Q_1w_1+Q_0w_0&=\theta(a_1).
\end{aligned}
$$

> **Proof.** The source polynomials have respectively 996, 48, 424 and 42 monomials. Their coefficients lift to rationals with denominators at most 18. `verify_stageA_Q.py` evaluates every displayed identity over $\mathbb Q$, checks the two gauge raising operators and the rotation raising operator, and verifies the multiplication homotopy on $b$, $\partial_1b$ and $\partial_2b$. Every residual is zero; the coefficient-level data are `stageA_Q.pkl` and `stageA_Q_certificate.json`. $\square$

For the two Cartan generators $\operatorname{diag}(1,-1,0)$ and $\operatorname{diag}(0,1,-1)$, write $c^1,c^2$ for the coefficients of $c$ and $b_1,b_2$ for the trace-dual coefficients of $b$. The coefficient of $b_1\partial_2c^1\partial_1c^1$ in $a_1$ is $1/3$. Cartan restriction annihilates every $Q_0$ boundary, hence $[a_1]\ne0$. Its outgoing higher maps require negative $N_b$. An incoming map of order $r\ge1$ would require $1+r$ copies of $b$ and one differentiated $c$ at $(2,2)$, which has minimum total weight $2(1+r)+1>4$. Thus $a_1$ represents a nonzero full perturbative class.

**Proposition 2.** The exact rational Stage B and C1 witnesses obey

$$
Q_0\eta(b_A)=3C_0^2\theta(b_A)-Q_1(x_3b_A),
\qquad
Q_0n_1=x_3a_1-C_0^2w_1=m_1.
$$

Here $A$ is the original component-engine basis label, used in the component identity for the derivation $\eta$.

> **Proof.** `solve_rational_witness.py` uses the original rational trace-product and rotation-projector descriptors. It selects 44 independent $Q_0$ images for $n_1$ and 73 for $\eta$, solves the resulting square systems over $\mathbb Q$, forms the full operators, and verifies the full polynomial equations. The resulting operators have 104210 and 82664 monomials. Their gauge and rotation tests pass, and reduction modulo $2147483647$ reproduces the original modular witnesses. The source coefficients exceed the range of unique reconstruction from that single prime; the exact linear solves supply the rational lift. The full witnesses and coefficients are in `n1_Q.pkl`, `eta_Q.pkl`, `n1_Q_certificate.json`, and `eta_Q_certificate.json`. $\square$

**Proposition 3.** The representative in the source's equation `eq-su3-m0` is an exact rational perturbative cocycle with nonzero $Q_0$ cohomology class.

> **Proof.** With the source's components
>
> $$
> m_0=\eta(a_1)-C_0^2w_0-2C_0\theta(w_1),
> $$
>
> exact rational polynomial evaluation gives
>
> $$
> Q_1m_1+Q_0m_0=0,
> \qquad
> Q_0(m_0-Q_1n_1)=0.
> $$
>
> The component counts are
>
> $$
> \#m_1=20228,\qquad \#m_0=95176,
> \qquad \#Q_1n_1=98586,
> \qquad \#(m_0-Q_1n_1)=98592.
> $$
>
> Here $\#$ counts nonzero canonical component monomials. The reduced representative contains seven differentiated $c$'s and zero $b$'s at $(5,5)$. Every $Q_r$, $r\ge1$, lowers $N_b$ by $r$, so its action on this representative vanishes. Consequently
>
> $$
> Q_{\mathrm{pert}}(m_0-Q_1n_1)=0.
> $$
>
> Its restriction to the Cartan subalgebra $\mathfrak t$ is
>
> $$
> \left.(m_0-Q_1n_1)\right|_{\mathfrak t}
> =40\,\partial_1c^1\,\partial_2c^1\,
> \partial_1c^2\,\partial_2c^2\,
> \left.\operatorname{tr}\bigl((\partial_1^2c)(\partial_1\partial_2c)(\partial_2^2c)\bigr)\right|_{\mathfrak t}.
> $$
>
> In the fixed Cartan basis, the cubic trace has coefficients $+1$ for the ordered basis-label triples $(1,1,2),(1,2,1),(2,1,1)$, coefficients $-1$ for $(1,2,2),(2,1,2),(2,2,1)$, and zero for $(1,1,1),(2,2,2)$. `verify_cartan_trace_Q.py` verifies the complete factorization over $\mathbb Q$. The Cartan restriction of $Q_1n_1$ is exactly zero. Thus $[m_0-Q_1n_1]\ne0$ in $H(Q_0)$. The full calculation is reproduced by `verify_partner_Q.py`; `partner_Q.pkl` contains the complete polynomials. $\square$

**Proposition 4.** The rational target has dimension one, and every incoming higher differential of order $r\ge2$ has zero source.

> **Proof.** The exact gauge-invariant dimension and modular rank lower bounds in the source census give
>
> $$
> \dim H(Q_0)^{(0,7)}
> \le968-669-298=1.
> $$
>
> Proposition 3 supplies a rational nonzero class, so equality holds and both subtracted $Q_0$ ranks are exact over $\mathbb Q$. For the incoming order-$r$ source $(N_b,N_c)=(r,6)$, six differentiated ghosts require
>
> $$
> 10-2r\ge6\quad\Longrightarrow\quad r\le2.
> $$
>
> For $r=2$, the gauge-invariant source dimension is 45. The independently reseeded computation `incoming_d2_Q0_rank.out` gives rank 45 for $Q_0:V_2^5\to V_2^6$, with zero failed invariance tests. The target has no outgoing $Q_0$ because seven differentiated ghosts would need more than six available derivatives. Therefore
>
> $$
> \dim H(Q_0)^{(2,6)}=45-45-0=0.
> $$
>
> Every later page at this source vanishes. $\square$

**Proposition 5.** The incoming $d_1:H(Q_0)^{(1,6)}\to H(Q_0)^{(0,7)}$ vanishes over $\mathbb Q$ at $(5,5)$.

> **Proof.** Take the coefficient of the canonical Cartan monomial
>
> $$
> \partial_2c^1\,\partial_2^2c^1\,\partial_1c^1\,
> \partial_1\partial_2c^1\,\partial_2c^2\,
> \partial_1c^2\,\partial_1^2c^2.
> $$
>
> It evaluates to 40 on the representative of Proposition 3 and vanishes on every $Q_0$ boundary. Since the target cohomology is one dimensional, this coefficient detects its entire cohomology.
>
> The scalar gauge-invariant source dimension is the exact character difference
>
> $$
> \dim V_1^6(5,5)_{\text{rotation spin }0}=1002-837=165.
> $$
>
> `verify_incoming_d1_Q.py` constructs 165 rational trace-product scalar operators. A minor of their coefficient sketches has
>
> $$
> \det\equiv469534650\pmod{2147483647}\ne0.
> $$
>
> Therefore they form a rational basis of the entire source. The script constructs a rational linear functional on the $Q_0$ target from 44 independent $Q_0$ images and verifies, on all 165 basis vectors, that its value on $Q_0$ equals the displayed Cartan coefficient of $Q_1$. Every equality is checked over $\mathbb Q$. Thus the Cartan coefficient of $Q_1$ vanishes on the entire $Q_0$ kernel. The target detector is injective on cohomology, so the incoming $d_1$ is zero. The basis descriptors, exact functional, minor and its coefficients are recorded in `incoming_d1_Q_certificate.json`. $\square$

**Theorem 6.** For the source's nilpotent perturbative differential, with its gauge/rotation equivariance and degree shifts, the SU(3) partner has

$$
\dim_{\mathbb Q} E_\infty^{(0,7)}(5,5)=1,
\qquad
[m_0-Q_1n_1]\ne0.
$$

> **Proof.** Proposition 3 gives an exact rational full perturbative cocycle, and Proposition 4 gives a one-dimensional $E_1$ target. Every outgoing map requires negative $N_b$. Proposition 5 removes the incoming $d_1$, and Proposition 4 removes every incoming higher differential. Finite filtration at the fixed weight identifies the surviving class with a nonzero perturbative cohomology class. $\square$

The selected source Massey representative is therefore exactly $\hbar[m_0-Q_1n_1]$ in these conventions. The nonperturbative differential and the full representative-ambiguity calculation are separate from the closure and nonexactness theorem proved here.

**Proposition 7.** The complete perturbative source of the $C_0^2$ indeterminacy in the partner degree vanishes over $\mathbb Q$. Consequently the partner remains nonzero modulo $C_0^2E_\infty$.

> **Proof.** Multiplication by $C_0^2$ adds $(2,2)$ to the weights and four to $N_c$. Its possible source for the partner has $(n_1,n_2)=(3,3)$ and $N_c=3$. Derivative counting gives
>
> $$
> 6-2N_b\ge3\quad\Longrightarrow\quad N_b\in\{0,1\}.
> $$
>
> All dimensions below are for the full gauge-invariant spaces, with every rotation representation included. Independently reseeded modular rank lower bounds give
>
> $$
> \dim H(Q_0)^{(0,3)}\le19-6-13=0,
> \qquad
> \dim H(Q_0)^{(1,3)}\le24-12-8=4.
> $$
>
> The rational cocycle $\partial_1\partial_2a_1$ has nonzero Cartan restriction. Thus
>
> $$
> 1\le\dim H(Q_0)^{(1,2)}\le19-6-12=1,
> $$
>
> fixing the subtracted rank $\operatorname{rank}_{\mathbb Q}Q_0^{(1,2)}=12$. The rank $\operatorname{rank}_{\mathbb Q}Q_0^{(1,3)}=8$ saturates its eight-dimensional target; $\operatorname{rank}_{\mathbb Q}Q_0^{(0,3)}=13$ follows from the acyclicity just displayed. The map $Q_0^{(2,2)}$ has zero target. The two-step block-rank identity therefore gives rational bounds
>
> $$
> \begin{aligned}
> \operatorname{rank}d_1^{(2,2)}&\ge14-0-12=2,\\
> \operatorname{rank}d_1^{(1,3)}&\ge23-8-13=2.
> \end{aligned}
> $$
>
> Since $d_1^2=0$,
>
> $$
> 0\le\dim E_2^{(1,3)}\le4-2-2=0.
> $$
>
> Both permitted $N_b$ sectors vanish on $E_2$, hence on $E_\infty$. Any well-defined quantum multiplication by $C_0^2$ therefore has zero image in the partner degree, including all its loop corrections. Together with Theorem 6,
>
> $$
> \dim_{\mathbb Q}\left(E_\infty/C_0^2E_\infty\right)^{(0,7)}(5,5)=1.
> $$
>
> The source of a possible $C_0^2$ multiple representing $a_1$ would have $N_c=2-4=-2$, so $[a_1]$ also lies outside that image. `verify_C0_squared_image.py` and `C0_squared_image_certificate.json` record the exact rational class and all seven independently reseeded rank computations. $\square$

The cubic trace in Proposition 3 is the polarized symmetric tensor of the primitive invariant $\operatorname{tr}(X^3)$ for an even Lie-algebra element $X$, evaluated on the three distinct odd jets in the displayed order. For $X=\operatorname{diag}(x^1,-x^1+x^2,-x^2)$,

$$
\operatorname{tr}(X^3)=3x^1x^2(x^1-x^2).
$$

Its polarization gives precisely the six signed coefficients in Proposition 3. This is a finite-rank invariant-theory check on the selected perturbative representative; the nonperturbative differential is a separate input.

**Proposition 8.** The SU(3) partner is annihilated by $C_0$ in perturbative cohomology, with a rational pure-ghost $Q_0$ primitive.

> **Proof.** Since both factors contain no $b$, the quantum product is the ordinary product. Direct rational multiplication gives 278026 terms and
>
> $$
> Q_0\bigl(C_0(m_0-Q_1n_1)\bigr)=0.
> $$
>
> The product lies in the rotation-scalar sector $(N_b,N_c)=(0,9)$ at $(6,6)$. Its exact scalar-space dimension and independently computed rank lower bounds are
>
> $$
> \dim V_0^9(6,6)_{\text{rotation spin }0}=267,
> \quad\operatorname{rank}Q_0^{(0,8)}\ge216,
> \quad\operatorname{rank}Q_0^{(0,9)}\ge51.
> $$
>
> These are valid rational lower bounds from the gauge and rotation projectors at prime $2147483647$, with zero failed invariance checks. The chain condition gives
>
> $$
> 0\le\dim H(Q_0)^{(0,9)}(6,6)_{\text{rotation spin }0}
> \le267-216-51=0.
> $$
>
> Thus
>
> $$
> C_0(m_0-Q_1n_1)\in Q_0V_0^8(6,6).
> $$
>
> Every higher $Q_r$ vanishes on the pure-ghost primitive because it has $N_b=0$. Therefore $C_0[m_0-Q_1n_1]=0$ in full perturbative cohomology. This argument uses no division by $\hbar$. The exact product, character counts and ranks are saved in `c0_product_Q.pkl`, `c0_product_dims.json` and `c0_product_ranks.json`. $\square$

**Proposition 9.** The candidate SU(3) $C_2$ sector $(N_b,N_c)=(2,2)$ at $(3,3)$ is completely lifted by $d_1$.

> **Proof.** Its gauge-invariant chain space has dimension four, incoming $Q_0$ rank at least two and zero outgoing $Q_0$ target. Proposition 7 gives the exact rational incoming $d_1$ rank at this source:
>
> $$
> \dim H(Q_0)^{(2,2)}\le4-2=2,
> \qquad
> \operatorname{rank}d_1^{(2,2)}\ge14-0-12=2.
> $$
>
> Therefore the kernel is zero and $E_2^{(2,2)}(3,3)=0$. In particular this sector supplies no $E_\infty$ generator named $C_2$. $\square$

**Proposition 10.** The exact quantum product $a_1\star a_1$ represents a nonzero full perturbative class at weight $(4,4)$, with leading sector $(N_b,N_c)=(2,4)$.

> **Proof.** The original product-forest kernels give
>
> $$
> a_1\star a_1=a_1^2+\hbar\mu_1(a_1,a_1)+\hbar^2\mu_2(a_1,a_1).
> $$
>
> The respective exact rational component polynomials contain 822, 8704 and 5216 terms. The one-loop and two-loop components agree coefficientwise with the independently implemented optimized source formulas. Direct full-polynomial evaluation gives
>
> $$
> Q_0(a_1^2)=0,
> \qquad Q_1(a_1^2)+Q_0\mu_1(a_1,a_1)=0,
> $$
>
> $$
> Q_2(a_1^2)+Q_1\mu_1(a_1,a_1)+Q_0\mu_2(a_1,a_1)=0.
> $$
>
> Each component at order $\hbar^n$ has $2-n$ copies of $b$. Since $Q_r$ removes $r$ copies of $b$, all remaining closure equations vanish termwise. At leading order, the gauge-invariant target chain dimension is 20, incoming $Q_0$ rank is at least 19, and the outgoing target is absent. The Cartan restriction of $a_1^2$ has coefficient $-2/3$ on each of the three canonically ordered monomials with bosonic factors $b_1^2,b_1b_2,b_2^2$ and ghost factors $\partial_2c^1\partial_1c^1\partial_2c^2\partial_1c^2$. Here $b_1,b_2$ are the coefficients of the trace-dual Cartan basis and $c^1,c^2$ are the coefficients of its Chevalley-dual basis. Thus
>
> $$
> 1\le\dim H(Q_0)^{(2,4)}(4,4)\le20-19=1.
> $$
>
> An incoming $d_r$, for $r\ge1$, would require a source with $(N_b,N_c)=(2+r,3)$. Its minimum total weight is $2(2+r)+3=7+2r>8$, so no such source exists. Closure of the full product excludes outgoing differentials. Hence $[a_1\star a_1]\ne0$ in $E_\infty$.
>
> Exact Cartan comparison in the one-dimensional classical target additionally gives the following equalities in $H(Q_0)$:
>
> $$
> [a_1^2]= -\frac16[\operatorname{tr}(b^2)C_0^2]
> $$
>
> $$
> =-\frac16\left[C_0\operatorname{tr}\!\left(
> b^2\partial_1c\partial_2c-b^2\partial_2c\partial_1c
> +b\partial_1c\,b\partial_2c\right)\right].
> $$
>
> These are cohomology equalities, not literal polynomial identities. The two differences have 921 and 1044 terms, respectively, and are $Q_0$-exact by the same one-dimensional target and nonzero Cartan detector. The exact products and residual certificates are `C1_square_Q.pkl`, `C1_square_Q_certificate.json`, `C1_star_Q.pkl` and `C1_star_Q_certificate.json`. $\square$

**Theorem 11.** At weight $(4,4)$, the full gauge-invariant sector with leading $(N_b,N_c)=(3,2)$ has

$$
\dim_{\mathbb Q}E_\infty^{(3,2)}(4,4)=1.
$$

An exact full cocycle is stored as $C_3$ in `C3_Q.pkl`. Its normalization is specified by the following trace identities; no prior normalization of a Claude $C_3$ is assumed.

The dimension statement concerns the associated-graded sector with leading $N_b=3$. Proposition 12 separately classifies the six constructed lifts and proves uniqueness in their complete rotation-scalar lift block. Non-scalar lower-filtration sectors are outside that statement.

> **Proof.** The exact coefficient polynomials at orders $1,\hbar,\hbar^2,\hbar^3$ have 434, 2434, 96 and zero terms. Independent readback from the saved trace descriptors verifies
>
> $$
> Q_0[\hbar^0]C_3=0,
> \qquad Q_0[\hbar]C_3+Q_1[\hbar^0]C_3=0,
> $$
>
> $$
> Q_0[\hbar^2]C_3+Q_1[\hbar]C_3+Q_2[\hbar^0]C_3=0,
> $$
>
> $$
> Q_1[\hbar^2]C_3+Q_2[\hbar]C_3+Q_3[\hbar^0]C_3=0.
> $$
>
> These checks use the full tree differential, including undifferentiated ghosts. Gauge and both rotation-generator checks vanish. Every higher equation vanishes by $b$-degree. In the unchanged $C_0,a_1$ trace convention, the literal leading coefficient is
>
> $$
> [\hbar^0]C_3=
> -\frac{111}{112}\operatorname{tr}(b^2)a_1
> +\frac{57}{56}\operatorname{tr}(b^3)C_0
> $$
>
> $$
> \hspace{2em}
> +\frac{15}{14}\left[
> \operatorname{tr}(b^2\partial_2c)\operatorname{tr}(b\partial_1c)
> -\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b\partial_2c)\right].
> $$
>
> Its nonzero Cartan restriction has 16 monomials. A direct rational primitive check gives the stronger literal identity
>
> $$
> [\hbar^0]C_3-
> \frac9{112}\operatorname{tr}(b^2)a_1
> +\frac3{56}\operatorname{tr}(b^3)C_0
> $$
>
> $$
> =Q_0\left\{\frac{15}{14}\operatorname{tr}\!\left[
> b^2\bigl((\partial_2b)(\partial_1c)-(\partial_1b)(\partial_2c)\bigr)\right]\right\}.
> $$
>
> The primitive has 256 terms and the exact difference has 354 terms. The last nonzero coefficient is also a short literal trace:
>
> $$
> [\hbar^2]C_3=\frac14\operatorname{tr}\!\left\{b\left[
> (\partial_2^3c)(\partial_1^3c)
> -3(\partial_1\partial_2^2c)(\partial_1^2\partial_2c)
> +3(\partial_1^2\partial_2c)(\partial_1\partial_2^2c)
> -(\partial_1^3c)(\partial_2^3c)
> \right]\right\}.
> $$
>
> Exact character counts give full leading chain dimension six. Direct sparse rational elimination on the exhaustive trace contractions gives three independent incoming $Q_0$ images and two independent Cartan restrictions of $Q_1$ images. Since Cartan restriction annihilates $Q_0$ boundaries, the latter supply a rational rank-two lower bound for $d_1$. Consequently
>
> $$
> \dim E_2^{(3,2)}(4,4)\le6-3-2=1.
> $$
>
> An incoming $d_r$ with $r\ge1$ would require $(N_b,N_c)=(3+r,1)$, whose minimum total weight is $7+2r>8$. The displayed full cocycle has nonzero leading Cartan restriction, so it cannot be a tree boundary and has no incoming differential. It supplies the matching lower bound. Thus the full $E_\infty$ dimension is exactly one over $\mathbb Q$. The independent witnesses are `C3_full_gauge_Q_certificate.json`, `C3_ranks_Q_certificate.json`, `C3_Q_certificate.json` and `C3_leading_relation_Q_certificate.json`. $\square$

**Proposition 12.** The complete rotation-scalar $C_3$ lift block at $(4,4)$, $N_c=2$, and $N_b+\text{loop order}=3$ has one-dimensional cohomology over $\mathbb Q$. All six saved lifts represent that same direction.

> **Proof.** The source weights and character calculation give scalar chain dimensions $3,15,8,0$ at $N_b=3,2,1,0$. Exact sparse rational elimination of the full differential gives rank18. The complete ghost-number-one source dimensions are $4,5,1,0$, and its full differential has rational rank7. Therefore
>
> $$
> \dim H(Q_{\rm pert})_{\text{specified scalar lift block}}
> =26-18-7=1.
> $$
>
> The relative leading Cartan coefficients of the six stored lifts are $1,1,-1,-10/27,2/3,-2/3$. After subtracting these multiples of the first lift, all differences are exact. The five nontrivial primitives have $(N_b=3,N_b=2)$ component term counts $(340,410),(340,410),(340,410),(340,278),(340,0)$, with zero lower components. Every saved primitive was independently re-evaluated under the full rational $Q$ and reproduces its complete difference coefficientwise.
>
> These checks are `C3_lift_differences_Q_certificate.json` and `C3_scalar_complex_Q_certificate.json`. The block assigns $\hbar^{3-N_b}$ to each component. It does not classify non-scalar sectors or other homogeneous $\hbar$ blocks. $\square$

**Proposition 13.** For the leading scalar sectors $(n_1,n_2)=(k+1,k+1)$, $(N_b,N_c)=(k,2)$, the exact $E_2$ dimensions for $k=1,\ldots,6$ are

| $k$ | Leading scalar chain dimension | Incoming $Q_0$ rank | $\operatorname{rank}(\mathrm{Cartan}\circ Q_1)$ | $\dim E_2$ |
|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 | 1 |
| 2 | 3 | 1 | 2 | 0 |
| 3 | 3 | 1 | 1 | 1 |
| 4 | 5 | 2 | 2 | 1 |
| 5 | 5 | 2 | 3 | 0 |
| 6 | 7 | 3 | 3 | 1 |

> **Proof.** Every leading $b$ is undifferentiated and the two ghosts are first jets. The scalar chain is the invariant part of the tensor product of the $k$th symmetric power and the second symmetric power of the adjoint representation. Exact integer characters give the table's chain dimensions. The Cayley--Hamilton identity
>
> $$
> b^3=\frac12\operatorname{tr}(b^2)b+\frac13\operatorname{tr}(b^3)\mathbf1
> $$
>
> reduces every bare-$b$ gap in a trace to lengths zero, one or two. The resulting trace covariants attain the full character dimensions by direct rational polynomial elimination. The incoming scalar source has one differentiated $b$ and one first-jet ghost; its complete reduced trace basis gives the displayed exact $Q_0$ ranks. The leading Cartan ranks equal the chain dimensions minus these ranks, so Cartan restriction is injective on the leading classical cohomology.
>
> Every leading vector is $Q_0$-closed because a target with $N_c=3$ would require total weight at least $2k+3>2k+2$. Cartan restriction of its $Q_1$ image supplies the displayed rational lower bound for $d_1$. Thus $k=2,5$ are killed. For $k=1,3$, the already constructed full cocycles provide the matching surviving direction. For $k=4,6$, direct rational first-lift equations give
>
> $$
> Q_0(\text{saved primitive})=Q_1(\text{saved leading coefficient}).
> $$
>
> The respective leading, $Q_1$-image and primitive term counts are $(864,11320,10270)$ and $(3214,68100,75826)$. These are full polynomial equalities, not tests of Cartan restrictions alone. An incoming $d_1$ would need $(N_b,N_c)=(k+1,1)$, whose minimum total weight exceeds $2(k+1)$. Hence each surviving candidate gives exactly one $E_2$ class.
>
> The certificates are `leading_tower_Q_certificate.json`, `tower_k4_first_lift_Q_certificate.json` and `tower_k6_first_lift_Q_certificate.json`. Compact, coefficientwise verified leading Cartan formulas are in `tower_cartan_trace_formulas.md`. The proposition asserts no higher survival for $k=4,6$ and no uncomputed congruence-pattern theorem. $\square$

**Theorem 14.** At weight $(5,5)$, the leading scalar $SU(3)$ sector $(N_b,N_c)=(4,2)$ has one-dimensional $E_\infty$ over $\mathbb Q$. Its complete rotation-scalar homogeneous lift block $N_b+\text{loop order}=4$ also has one-dimensional full perturbative cohomology.

> **Proof.** The saved representative `k4_full/C4_Q.pkl` has the following literal leading coefficient, checked coefficientwise using the original $a_1$:
>
> $$
> [\hbar^0]C_4=\operatorname{tr}(b^3)a_1
> -\frac32\operatorname{tr}(b^2\partial_1c\,b^2\partial_2c)
> -\frac32\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c).
> $$
>
> The coefficients at orders $1,\hbar,\hbar^2,\hbar^3,\hbar^4$ contain $864,10728,14686,0,0$ component monomials. All $124$ scalar correction coordinates were included simultaneously. The exact full-map column solve has zero residual. A fresh process reconstructed every nonzero correction coefficient from its original trace descriptor and recomputed the full equations
>
> $$
> Q_0[\hbar^0]C_4=0,
> \qquad Q_0[\hbar]C_4+Q_1[\hbar^0]C_4=0,
> $$
>
> $$
> Q_0[\hbar^2]C_4+Q_1[\hbar]C_4+Q_2[\hbar^0]C_4=0,
> $$
>
> $$
> Q_1[\hbar^2]C_4+Q_2[\hbar]C_4+Q_3[\hbar^0]C_4=0,
> $$
>
> $$
> Q_2[\hbar^2]C_4+Q_3[\hbar]C_4+Q_4[\hbar^0]C_4=0.
> $$
>
> These checks include undifferentiated ghosts in $Q_0$, both rotation generators, and gauge invariance. Higher equations vanish by $b$-degree. The nonzero leading Cartan restriction has $20$ terms, and Proposition13 proves that it spans the one-dimensional leading scalar $E_2$. Any incoming $d_r$ with $r\ge1$ would require $(N_b,N_c)=(4+r,1)$, whose minimum total weight $9+2r$ exceeds $10$. Thus the full cocycle gives the matching lower bound for $E_\infty$.
>
> The complete ghost-number-two scalar chain dimensions at $N_b=4,3,2,1,0$ are $5,43,67,13,1$, totaling $129$. The complete ghost-number-one dimensions are $6,19,9,1,0$, totaling $35$. Each trace basis attains its independently computed exact integer-character dimension.
>
> A nonzero rational minor proves rank at least $95$ for the full correction map. The two independent leading images detected by $\mathrm{Cartan}\circ Q_1$ cannot lie in this correction image: its $N_b=3$ projection is a $Q_0$ image, annihilated by Cartan restriction. A second nonzero rational minor proves incoming full-$Q$ rank at least $31$. Consequently
>
> $$
> 1\le\dim H(Q_{\rm pert})_{\text{specified scalar lift block}}
> \le129-(95+2)-31=1.
> $$
>
> Both ranks therefore saturate: the outgoing rank is $97$ and the incoming rank is $31$. All lifts with the same leading class differ by a full-$Q$ boundary within this homogeneous scalar block. No non-scalar lower sector, other homogeneous $\hbar$ block, or nonperturbative differential is classified here.
>
> The exact certificates are `k4_full/C4_Q_certificate.json`, `k4_full/independent_full_Q_readback.json`, `k4_full/scalar_ambiguity_Q_certificate.json`, and `k4_full/incoming_rank_Q_certificate.json`. The $k=6$ second lift is treated separately in Theorem17. $\square$

**Theorem 15.** For every $k\ge0$, the leading rotation-scalar tree cohomology at $(N_b,N_c)=(k,2)$ and weight $(k+1,k+1)$ is the degree-$k$ part of a free module over $\mathbb Q[\operatorname{tr}(b^2),\operatorname{tr}(b^3)]$, with generators $C_0,a_1,\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ of degrees $0,1,2$.

> **Proof.** Invariant trace generation and the source Cayley--Hamilton identity reduce every leading scalar invariant to a combination of six single-trace and three double-trace covariants with bare-$b$ gaps at most two. All nine exact reductions, including full-$Q_0$ primitives for the differences, are recorded in `leading_module_Q.pkl` and displayed in `leading_module_and_radial_kernel.md`. Their maximum degree is four. Multiplication by the $Q_0$-closed pure-$b$ traces extends these identities to every degree.
>
> The three resulting Cartan tensors have source-coordinate determinant $-3b_1b_2(b_1+b_2)$, where $b_1,b_2$ are the original trace-dual Cartan coefficients. Since this polynomial is nonzero and every $Q_0$ boundary has zero Cartan restriction, no polynomial module relation remains. Thus, with $t$ counting $b$-degree,
>
> $$
> \frac{1+t+t^2}{(1-t^2)(1-t^3)}=\frac1{(1-t)(1-t^2)},
> \qquad
> \dim H(Q_0)^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
> =\left\lfloor\frac k2\right\rfloor+1.
> $$
>
> These finite exact reductions establish the global module statement without assuming injectivity from generic Cartan localization. $\square$

**Theorem 16.** For every $k\ge0$,

$$
\dim_{\mathbb Q}E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}
0,&k\equiv2\pmod3,\\
1,&k\not\equiv2\pmod3.
\end{cases}
$$

> **Proof.** The source $Q_1$ has differential order two. Its eighteen values on the three module generators multiplied by the pure-$b$ trace monomials of ordinary polynomial degree at most two determine its action on arbitrary polynomial module coordinates. After Cartan restriction, every second-coordinate-derivative coefficient vanishes. The third module coordinate is forced to zero by its exclusive second-ghost-jet term. The remaining two independent target operators have a source-coordinate minor $-216b_1$.
>
> The $C_0$ and $a_1$ polynomial coefficients therefore satisfy
>
> $$
> \frac{\partial(\text{coefficient of }C_0)}{\partial\operatorname{tr}(b^2)}
> +\frac{\operatorname{tr}(b^2)}4
> \frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)}=0,
> $$
>
> $$
> \frac{\partial(\text{coefficient of }C_0)}{\partial\operatorname{tr}(b^3)}
> +\frac23\frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)}=0.
> $$
>
> Their compatibility equates the second derivative in $\operatorname{tr}(b^2)$ of the $a_1$ coefficient with $3\operatorname{tr}(b^2)/8$ times its second derivative in $\operatorname{tr}(b^3)$. At weighted degree $k-1$, the coefficient recurrence in `leading_module_and_radial_kernel.md` has no nonzero solution for $k=2\pmod3$ and exactly one solution direction otherwise. The $C_0$ coefficient is then uniquely fixed for $k>0$ by weighted Euler differentiation; $k=0$ gives $C_0$ itself.
>
> This Cartan-kernel calculation is promoted to full $E_2$ by five exact primitive identities: the three pure-Casimir brackets and the two first-derivative relations in `radial_first_primitives_Q.pkl`. The order-two chain rule expresses the full $Q_1$ image of any kernel solution as $Q_0$ of the corresponding polynomial derivative combination of those five primitives. Independent readback of all five identities and direct checks on a $k=7$ kernel solution and an off-kernel polynomial pass in `../su3_transgression/verify_universal_first_lift.json`.
>
> Incoming $d_1$ sources would have $(N_b,N_c)=(k+1,1)$ and minimum total weight $2k+3>2k+2$, so none exist. This proves the stated exact $E_2$ dimensions. No all-$k$ $E_\infty$ or nonperturbative survival is asserted. $\square$

**Theorem 17.** At weight $(7,7)$, the leading rotation-scalar $SU(3)$ sector $(N_b,N_c)=(6,2)$ has one-dimensional $E_3$ over $\mathbb Q$.

> **Proof.** The saved leading coefficient represents the following tree class, where square brackets denote $Q_0$ cohomology:
>
> $$
> [\text{leading coefficient}]
> =\left[3\operatorname{tr}(b^2)\operatorname{tr}(b^3)a_1
> -\left(\frac{\operatorname{tr}(b^2)^3}{4}
> +\operatorname{tr}(b^3)^2\right)C_0\right].
> $$
>
> Their Cartan restrictions agree coefficientwise. The complete leading tree module of Theorem15 makes this an equality of tree classes. The original leading polynomial is retained without rescaling.
>
> The complete first- and second-correction scalar spaces have dimensions $118$ and $594$. The first space has full-$Q_0$ rank $56$ and kernel dimension $62$. The original fixed first correction gives a $1314$-term Cartan obstruction; the Cartan $Q_1$ image of all closed first-correction freedom has rank $8$ and cancels it. The subsequent joint solve retains all $118+594=712$ correction coordinates.
>
> A nonzero $487\times487$ rational determinant after a fixed linear projection selects independent map columns. The resulting rational coefficients satisfy the complete unprojected polynomial equations. A fresh process then reconstructs every nonzero correction from its original trace descriptor and independently verifies
>
> $$
> Q_0(\text{leading coefficient})=0,
> $$
>
> $$
> Q_0(\text{first correction})
> +Q_1(\text{leading coefficient})=0,
> $$
>
> $$
> Q_0(\text{second correction})
> +Q_1(\text{first correction})
> +Q_2(\text{leading coefficient})=0.
> $$
>
> The three source coefficients contain $3214,76234,355880$ component monomials. The full-$Q_0$ evaluations include bare ghosts; gauge invariance and both rotation generators were also checked on all three coefficients.
>
> Theorem16 gives a one-dimensional leading $E_2$, and the displayed second lift makes its outgoing $d_2$ zero. An incoming $d_r$ with $r\ge1$ would require $(N_b,N_c)=(6+r,1)$ and minimum total weight $13+2r>14$. Hence
>
> $$
> \dim_{\mathbb Q}E_3^{(6,2)}(7,7)_{\mathrm{scalar}}=1.
> $$
>
> The witnesses and fresh readback are `k6_two_loop/joint_two_loop_Q.pkl`, `k6_two_loop/joint_two_loop_Q_certificate.json`, and `k6_two_loop/independent_joint_two_loop_readback.json`. No third- or higher-loop equation, full-$k=6$ lift, or uniqueness in the complete $1848$-coordinate lower scalar block is asserted. This is an $E_3$ result, not an $E_\infty$ or physical nonperturbative result. $\square$
