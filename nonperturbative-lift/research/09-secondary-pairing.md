<a id="np-section-07"></a>

## 9. From multiplicative extension to the nonperturbative differential

In pure $SU(N)$ SYM, the gaugino bilinear $C_0$ satisfies $[C_0]^N=0$ in perturbative cohomology.

The pure-ghost primitive $x_N$ satisfies

$$
Q_0x_N=C_0^N,
\qquad Q_rx_N=0\quad(r\ge1).
$$

Here $Q_r$ is the $r$-loop supercharge, and $Q_{\rm pert}$ is the full perturbative sum.

Write $L=C_0\star{}$ for quantum left multiplication and $h$ for an odd multiplication homotopy.

A multiplication homotopy gives a primitive for repeated multiplication on a closed operator.

Explicitly, $L\mathcal O=C_0\star\mathcal O$, and $L^N\mathcal O$ means $N$ successive applications of $L$.

The operator $h$ changes Grassmann parity.

Its perturbative identity is $Q_{\rm pert}h+hQ_{\rm pert}=L^N$.

Thus $Q_{\rm pert}\mathcal O=0$ gives $Q_{\rm pert}(h\mathcal O)=L^N\mathcal O$.

The instanton scale is $\Lambda^{3N}$; $\nu$ denotes its scalar coefficient in the source normalization.

The theorem below assumes the exact product and homotopy identities

$$
[Q,L]=0,
\qquad
\{Q,h\}=L^N-\nu\Lambda^{3N}\mathbf1,
\qquad \nu\ne0.
$$

Here $Q$ is the exact supercharge; at $\Lambda=0$, $L,h$ reduce to the corresponding perturbative operators.

This hypothesis specifies the compatibility of the local product, composite representative, and chain homotopy.

The brackets here are operator commutators: $[Q,L]=QL-LQ$ and $\{Q,h\}=Qh+hQ$; $\mathbf1$ is the identity map.

For an operator closed under the full $Q$ $\mathcal O$, the hypothesis makes $(L^N-\nu\Lambda^{3N})\mathcal O$ exact, with primitive $h\mathcal O$.

It simultaneously specifies multiplication and its primitive at one insertion point.

<a id="np-b082a68ee6a1c259-1"></a>

### The full iterated associator

Let $H_{A,B}$ be the source associator homotopy.

For the closed pure-ghost multipliers $A,B$ used here, the associator homotopy compares two ways to multiply an input operator: first form the ordinary product $AB$ and multiply by it, or apply quantum multiplication by $B$ and then by $A$.

The difference is a supercharge variation.

The next identity states this comparison as an equality of operators acting on every allowed input.

On the closed pure-$c$ inputs used here, it satisfies

$$
\{Q_{\rm pert},H_{A,B}\}
=(AB)\star{}-A\star(B\star{}).
$$

**Proposition 9.1.**

A perturbative multiplication homotopy is

$$
h_{\rm diag}=x_N\star{}
-\sum_{\substack{j=1,\ldots,N-1\\\mathrm{successive\ multiplication\ homotopies}}}
H_{C_0^j,C_0}L^{N-j-1},
\qquad
\{Q_{\rm pert},h_{\rm diag}\}=L^N.
$$

> **Proof.**
>
> The quantum Leibniz identity gives
>
> $$
> \{Q_{\rm pert},x_N\star{}\}=C_0^N\star{},
> \qquad [Q_{\rm pert},L]=0.
> $$
>
> The successive homotopies telescope:
>
> $$
> \begin{aligned}
> &\left\{Q_{\rm pert},
> \sum_{\substack{j=1,\ldots,N-1\\\mathrm{successive\ multiplication\ homotopies}}}
> H_{C_0^j,C_0}L^{N-j-1}\right\}\\
> &=\sum_{\substack{j=1,\ldots,N-1\\\mathrm{successive\ multiplication\ homotopies}}}
> \left[(C_0^{j+1}\star{})L^{N-j-1}-(C_0^j\star{})L^{N-j}\right]\\
> &=C_0^N\star{}-L^N.
> \end{aligned}
> $$
>
> Subtracting these identities proves the assertion.
>
> $\square$

For $SU(3)$, both correction terms occur:

$$
h_{\rm diag}=x_3\star{}-H_{C_0^2,C_0}-H_{C_0,C_0}L.
$$

The difference between the power of the quantum multiplication operator and multiplication by the ordinary power is explicit.

On the gauge-invariant input $\operatorname{tr}b^2$,

$$
[\hbar]\left[(C_0^2)\star\operatorname{tr}b^2
-C_0\star(C_0\star\operatorname{tr}b^2)\right]\ne0.
$$

The exact component polynomial contains $2064$ nonzero monomials.

Here $[\hbar]$ extracts the coefficient linear in the loop parameter.

<a id="np-90f763e6956a69ff-1"></a>

### Secondary classes and instanton pairing

Let $E_\infty^{\rm pert}$ denote the stable perturbative cohomology page.

For a class $a$ satisfying $C_0a=0$, choose a cocycle representative $a$ and a primitive $w$ with

$$
Q_{\rm pert}a=0,
\qquad Q_{\rm pert}w=La.
$$

The equation $Q_{\rm pert}w=La$ records how multiplication by $C_0$ becomes zero on the class of $a$: it supplies the specific primitive $w$.

Applying the homotopy to $a$ gives $ha$; applying the remaining $N-1$ multiplications to $w$ gives $L^{N-1}w$.

Their supercharge variations agree.

Their difference is therefore a new closed operator, called the secondary representative.

The associated Massey-type secondary class is

$$
\langle C_0^{N-1},C_0,a\rangle
=[ha-L^{N-1}w]
\quad\text{in}\quad
E_\infty^{\rm pert}/C_0^{N-1}E_\infty^{\rm pert},
$$

In this expression, $h,L$ are their perturbative values.

Brackets denote cohomology classes; the denominator records the ambiguity in the choice of $w$.

> **Proof: closure and independence of representatives.**
>
> $$
> Q_{\rm pert}(ha-L^{N-1}w)
> =L^Na-hQ_{\rm pert}a-L^{N-1}La=0.
> $$
>
> The replacement $w\mapsto w+y$, with $Q_{\rm pert}y=0$, changes the representative by $-L^{N-1}y$.
>
> The replacements $a\mapsto a+Q_{\rm pert}u$ and $w\mapsto w+Lu$ change it by
>
> $$
> hQ_{\rm pert}u-L^Nu=-Q_{\rm pert}(hu).
> $$
>
> Thus the quotient class is independent of both choices.
>
> $\square$

Let $Q_{1\text{-inst}}$ denote the first differential of the instanton-number filtration: the action induced on perturbative cohomology by the $\Lambda^{3N}$ coefficient of the exact $Q$.

The instanton-number filtration groups the exact supercharge by powers of $\Lambda^{3N}$.

On a perturbatively closed operator, its first nonzero instanton correction can obstruct exact closure.

Passing this correction to perturbative cohomology defines $Q_{1\text{-inst}}$.

A nonzero image identifies a pair removed at that instanton order.

In the theorem below the image is $-\nu a$ when $Q_{1\text{-inst}}a=0$; removal requires $\nu a$ to be nonzero in the displayed quotient.

<a id="np-93111f11b210145b-1"></a>

**Conditional Theorem 9.2.**

Under the exact-product hypothesis above,

$$
Q_{1\text{-inst}}\langle C_0^{N-1},C_0,a\rangle
+\langle C_0^{N-1},C_0,Q_{1\text{-inst}}a\rangle
=-\nu a
\quad\bmod C_0^{N-1}E_\infty^{\rm pert}.
$$

> **Proof.**
>
> Extend $a,w$ to formal representatives $\widetilde a,\widetilde w$.
>
> In the quotient modulo $\Lambda^{6N}$, write
>
> $$
> Q\widetilde a=\Lambda^{3N}e,
> \qquad
> L\widetilde a=Q\widetilde w+\Lambda^{3N}u.
> $$
>
> The zero-instanton coefficient of $e$ represents $Q_{1\text{-inst}}a$.
>
> Applying $Q$ and using $[Q,L]=0$ gives
>
> $$
> Q_{\rm pert}u=Le.
> $$
>
> Thus $C_0$ also annihilates $Q_{1\text{-inst}}a$, and
>
> $$
> \begin{aligned}
> Q(h\widetilde a-L^{N-1}\widetilde w)
> &=(L^N-\nu\Lambda^{3N})\widetilde a
> -h\Lambda^{3N}e-L^{N-1}(L\widetilde a-\Lambda^{3N}u)\\
> &=-\Lambda^{3N}(\nu\widetilde a+he-L^{N-1}u).
> \end{aligned}
> $$
>
> Taking the cohomology class of the $\Lambda^{3N}$ coefficient proves the assertion.
>
> $\square$

If the charges force $Q_{1\text{-inst}}a=0$ and $a$ is nonzero in this quotient, then

$$
Q_{1\text{-inst}}\langle C_0^{N-1},C_0,a\rangle=-\nu a.
$$

The perturbative annihilation and indeterminacy conditions hold for $C_2$ in $SU(2)$ and $C_1$ in $SU(3)$.

Their nonzero secondary representatives supply incoming classes.

The physical application of this identity assumes the exact-product hypothesis.

The relation $[C_0]^N=\nu\Lambda^{3N}$ is a multiplicative extension.

The primitive $x_N$ has the nonzero perturbative differential $Q_0x_N=C_0^N$.

The theorem converts this multiplicative input into a differential on cohomology classes.

<a id="np-6570a44df50d9296-1"></a>

### The formal loop parameter and matric secondary products

The secondary operation retains its loop degree.

For $SU(2)$, write the full $C_2$ cocycle as $a_2+\hbar a_1+\hbar^2a_0$ and $w=w_2+\hbar w_1+\hbar^2w_0$.

Choose $n_2$ with $Q_0n_2=x_2a_2-C_0w_2$.

Then

$$
ha-Lw-Q_{\rm pert}n_2
=\hbar M_1+\hbar^2M_0,
$$

$$
M_1=x_2a_1+\eta(a_2)-C_0w_1-\theta(w_2)-Q_1n_2.
$$

For $SU(3)$, $C_1=a_1$ and $w=w_1+\hbar w_0$.

Removing the leading boundary with $Q_0n_1=x_3a_1-C_0^2w_1$ gives the secondary representative

$$
\hbar\bigl[\eta(a_1)-C_0^2w_0-2C_0\theta(w_1)-Q_1n_1\bigr].
$$

The right-hand sides $-\nu[C_2]$ and $-\nu[C_1]$ in the conditional NP theorem correspond to these full secondary classes.

Over $\mathbb Q[[\hbar]]$, $\hbar$ is a nonunit.

Expressing the NP image of a normalized partner requires a specified localization or physical loop convention.

**Definition 9.5 (Parameter torsion).**

A cohomology class has $\hbar$-torsion if multiplication by a positive power of $\hbar$ gives zero.

Localizing the coefficient ring at $\hbar$ sends every such class to zero.

A unit is a coefficient with a multiplicative inverse in the chosen ring.

A formal series in $\mathbb Q[[\hbar]]$ is invertible exactly when its constant coefficient is nonzero.

Multiplication by $\hbar$ shifts every coefficient to a positive power, so $\hbar^{-1}$ lies outside this ring.

Localizing at $\hbar$ means adjoining that inverse.

Accordingly, the equations for $\hbar M_1+\hbar^2M_0$ retain their displayed loop factors until a coefficient ring or physical specialization is chosen.

We retain the displayed $\hbar$ factors.

The same exact-product matching determines the loop dependence of $\nu$.

Since $C_0C_4\ne0$, the partner of $C_4$ uses a matric secondary product.

In this subsection, retain the source operator notation $N=C_2\star{}$: $N$ is a multiplication operator, and the gauge group is fixed to $SU(2)$.

Choose $v$ and an odd operator $\kappa$ satisfying

$$
Q_{\rm pert}v=LC_4+\frac13NC_2,
\qquad
\{Q_{\rm pert},\kappa\}=LN-NL.
$$

The first identity uses the established product relation $[C_2]^2=-3[C_0][C_4]$.

The matric secondary representative is

$$
M_4=hC_4-Lv+\frac13(Nw+\kappa C_2).
$$

> **Proof: closure.**
>
> Preserving operator order gives
>
> $$
> \begin{aligned}
> Q_{\rm pert}M_4
> &=L^2C_4-L\left(LC_4+\frac13NC_2\right)
> +\frac13\left(NLC_2+(LN-NL)C_2\right)\\
> &=0.
> \end{aligned}
> $$
>
> The commutativity correction $\kappa$ cancels $LN-NL$ exactly.
>
> $\square$

A full rational cocycle exists in this weight and field-count block, and the specified scalar sector is one-dimensional.

Identifying the matric formula with that cocycle coefficient by coefficient requires matching the higher operations and homotopy choices.

If $L,N,h,\kappa$ and the chosen primitives admit compatible exact continuations, and the obstruction classes vanish in the source charge sectors, the same calculation at one-instanton order gives

$$
Q_{1\text{-inst}}[M_4]=-\nu[C_4],
\qquad
Q_{1\text{-inst}}[C_0M_4]=-\nu[C_0C_4].
$$

The original higher-tower calculation supplies finite-field and leading-component evidence for

$$
[C_2][C_4]=-\frac53[C_0][C_6],
$$

It also gives a leading $b^5J$ representative in the incoming sector for $C_6$.

Here

$$
J=(\partial_2c)(\partial_2^2c)(\partial_1c)
(\partial_1\partial_2c)(\partial_1^2c)
$$

uses the one-dimensional Cartan ghost.

The rational survival of this leading sector and the computed corrections appear above.

The full $C_6$ product identity and complete matric matching enter the higher-tower conjecture.

With the source rescaling $D_0=-C_0$, $D_j=(2j-1)C_{2j}$ for $j\ge1$, the first exact relation is

$$
[D_1]^2=[D_0][D_2].
$$

**Conjecture 9.3: the full tower.**

In the source normalization with leading components satisfying $C_0a_{2i+2j}=a_{2i}a_{2j}$,

$$
[D_i][D_j]=[D_0][D_{i+j}],
\qquad i,j\ge1.
$$

The predicted ratio of the Cartan coefficients of the partner and the leading $C_{2j}$ is

$$
-2(2j+1)(2j+3).
$$

The first three finite calculations give $-30,-70,-126$.

Extending them requires full quantum products and compatible higher homotopies.

Assume the conjectural exact ring relations and the exact-product hypothesis, with $\nu$ invertible in the chosen coefficient ring.

After inverting $\Lambda$, $[C_0]$ is invertible; once $[C_2]=0$, the ring relations successively make every remaining $[C_{2j}]$ vanish.

<a id="np-17abaa189a90b8dc-1"></a>

### Coherence and finite homotopies

A contraction between the kernel and quotient requires an additional coherence identity.

Direct calculation gives

$$
L(ha-L^{N-1}w)=[L,h]a-Q_{\rm pert}(hw).
$$

Thus the secondary operation takes values in $\ker L$ when $[L,h]$ acts trivially on the relevant class.

The strict commutator of the common-root associator is nonzero on general inputs.

The explicit differences satisfy the following identities, where $a_2$ is the leading component of $C_2$ and $\mu_1$ is the one-loop product coefficient:

$$
[C_0,H^{(1)}_{C_0,C_0}]a_2=Q_0u,
$$

$$
\eta_{\rm algebraic}(a_2)
-\bigl(\mu_1(x_2,a_2)-H^{(1)}_{C_0,C_0}(a_2)\bigr)
=Q_0u'.
$$

The exact rational primitives $u,u'$ contain $2092,1298$ component monomials, respectively.

The corresponding difference for the leading representative of $C_4$ has an explicit primitive with $10598$ monomials.

For $SU(3)$, the difference between the full iterated expression and the algebraic $\eta$ is a pure-ghost $Q_0$ boundary.

Every higher $Q_r$ annihilates pure-ghost polynomials, so this difference is also a $Q_{\rm pert}$ boundary.

The partner coefficient in the source Cartan basis remains $40$.

A second construction solves for $h$ in the algebra of operators with pure-$c$ coefficients and derivatives with respect to $b$ jets.

This algebra is graded commutative and preserved by $[Q_{\rm pert},-]$.

> **Proof.**
>
> The input and output field counts of an $n$-loop fan are
>
> $$
> (n,1;0,2),\qquad(n+1,0;1,1),
> $$
>
> The four entries count input $b$, input $c$, output $b$, and output $c$, respectively.
>
> Hence $Q_n$ differentiates at most one $c$, and its coefficient is at most linear in $b$.
>
> Commuting it with pure-$c$ multiplication or a $b$ derivative produces another $b$-differential operator with pure-$c$ coefficients.
>
> The $b$ derivatives commute with one another and annihilate these coefficients.
>
> Operator composition is therefore graded commutative.
>
> For an odd operator $h$ in this algebra,
>
> $$
> [L,h]=0,
> \qquad h^2=-h^2=0.
> $$
>
> The last equality uses characteristic zero.
>
> $\square$

For $SU(2)$, write $[\hbar^0]h=x_2\cdot{}$, $[\hbar]h=\eta$, and $L=C_0\cdot{}+\hbar\theta+\hbar^2L^{(2)}+\cdots$.

The solved recursive equations include

$$
\{Q_0,\eta\}=2C_0\theta-\{Q_1,x_2\cdot{}\},
$$

$$
\{Q_0,[\hbar^2]h\}
=2C_0L^{(2)}+\theta^2
-\{Q_1,\eta\}-\{Q_2,x_2\cdot{}\},
$$

$$
\begin{aligned}
\{Q_0,[\hbar^3]h\}
={}&2C_0L^{(3)}+\theta L^{(2)}+L^{(2)}\theta\\
&-\{Q_1,[\hbar^2]h\}-\{Q_2,\eta\}-\{Q_3,x_2\cdot{}\},
\end{aligned}
$$

$$
\begin{aligned}
\{Q_0,[\hbar^4]h\}
={}&2C_0L^{(4)}+\theta L^{(3)}+L^{(3)}\theta+(L^{(2)})^2\\
&-\{Q_1,[\hbar^3]h\}-\{Q_2,[\hbar^2]h\}
-\{Q_3,\eta\}-\{Q_4,x_2\cdot{}\}.
\end{aligned}
$$

| Coefficient | Solved total input jet degree | Number of nonzero component coefficients |
|---|---:|---:|
| $[\hbar]h$ | $0\le m\le6$ | $84$ |
| $[\hbar^2]h$ | $0\le m\le4$ | $324$ |
| $[\hbar^3]h$ | $0\le m\le2$ | $163$ |
| $[\hbar^4]h$ | $m=0$ | $11$ |

Here $m$ is the total derivative degree of the differentiated $b$-jet slots.

The fourth-order target contains $20525$ monomials.

The full rotation-scalar source dimension and the rank of its $Q_0$ image are both $14$.

The resulting coefficient contains $3159$ monomials; an independent scalar solve gives the same result.

The full $Q_0$ equation, both rotation generators, and ghost spectators have been checked in all $76$ jet sectors.

The homotopy identity therefore holds through fourth order on every polynomial in undifferentiated $b$, with arbitrary ghost spectators.

Extending this construction to every jet and loop order, and identifying it with the exact nonperturbative homotopy, remain open problems.

The finite projector that removes the lower-spin ambiguity of a primitive also has a closed form.

Let $J_+,J_-$ be the source rotation raising and lowering operators, $n$ the operator arity, and $m$ the total input jet degree.

On weight-zero coefficient maps, the invariant projector is

$$
\prod_{\substack{j\in\{1,\ldots,n+2+m\}\\\mathrm{rotation\ spins}}}
\left(1-\frac{[J_-,[J_+,-]]}{j(j+1)}\right).
$$

> **Proof.**
>
> On the weight-zero line of the spin-$j$ representation, $[J_-,[J_+,-]]$ has eigenvalue $j(j+1)$.
>
> The displayed polynomial is therefore the identity on spin zero and vanishes on every allowed nonzero integral spin.
>
> Since $Q_0$ commutes with rotations, applying this projector preserves the tree equation when its right-hand side is equivariant.
>
> $\square$

<a id="np-317a4949e0efe6ef-1"></a>

### Finite ghost-contraction operators

Proposition 9.4.

On the full jet algebra of $\mathfrak{sl}_N$, consider a ghost differential operator of finite jet support and positive differential order.

If its coefficients are polynomials in finitely many ghost jets and it graded-commutes with $Q_0$, its positive-order part vanishes.

Its multiplication part is $Q_0$-closed.

> **Proof.**
>
> Choose a nonzero coefficient of highest differential order.
>
> Choose an integer $t$ greater than twice the first-direction derivative order of every coefficient jet and differentiation jet.
>
> Write the source structure constants as $f^i{}_{jk}$, where $i,j,k$ label Lie algebra components.
>
> For each active jet $(r,s)$, the tree differential contains
>
> $$
> \left.Q_0(\partial_1^{r+t}\partial_2^s c^i)\right|_
> {(\partial_1^r\partial_2^s c)(\partial_1^tc)}
> =\binom{r+t}{r}f^i{}_{jk}
> (\partial_1^r\partial_2^s c^j)(\partial_1^tc^k).
> $$
>
> In the commutator, contracting the first factor produces the new differentiation slot $(r+t,s)$ and the fresh coefficient jet $\partial_1^tc^k$.
>
> The $Q_0$ variation of an original coefficient preserves the original differentiation slots, so its contribution has a distinct slot from this new term.
>
> Extracting each fresh-jet coefficient forces the color coefficient to commute with all of $\mathfrak{sl}_N$, whose center is zero.
>
> The highest-order part vanishes.
>
> Descending in differential order proves the assertion.
>
> $\square$

For gauge-invariant test polynomials in the reduced positive-ghost-jet algebra, pair each active slot with a distinct fresh spectator jet in an invariant trace.

The same argument applies.

The positive-jet hypothesis is essential: with the bare ghost retained, $[Q_0,\partial/\partial c^i]$ is a constant-gauge derivation and vanishes on gauge invariants.

The proposition concerns finite-support operators.

Ordinary field transformations are locally finite sums over all jets.
