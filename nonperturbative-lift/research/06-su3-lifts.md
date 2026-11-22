<a id="np-section-04"></a>

## 6. SU(3): first partner, products, and higher lifts

<a id="np-3e7812cbaa09bd65-1"></a>

### The complete $C_1$ partner

For the preceding $a_1$, the source primitives satisfy

$$
Q_0x_3=C_0^3,\qquad Q_0w_1=C_0a_1,
\qquad Q_1w_1+Q_0w_0=\theta(a_1),
$$

$$
Q_0n_1=x_3a_1-C_0^2w_1.
$$

Here $\theta$ is the one-loop correction to multiplication by $C_0$, and $\eta$ is the coefficient of the source odd multiplication homotopy.

The source defines

$$
m_0=\eta(a_1)-C_0^2w_0-2C_0\theta(w_1).
$$

**Proposition P5.**

The expression $m_0-Q_1n_1$ is a complete rational perturbative cocycle, with

$$
\dim_{\mathbb Q}E_\infty^{(0,7)}(5,5)=1,
\qquad [m_0-Q_1n_1]\ne0.
$$

> **Proof.**
>
> The full component calculation gives
>
> $$
> Q_1(x_3a_1-C_0^2w_1)+Q_0m_0=0,
> \qquad Q_0(m_0-Q_1n_1)=0.
> $$
>
> The final representative contains $98592$ terms, seven differentiated ghosts, and zero copies of $b$.
>
> Consequently, every $Q_r$ with $r\ge1$ annihilates it.
>
> For Chevalley Cartan coefficients $c^1,c^2$, its Cartan restriction is
>
> $$
> 40(\partial_1c^1)(\partial_2c^1)(\partial_1c^2)(\partial_2c^2)
> \left.\operatorname{tr}\bigl((\partial_1^2c)(\partial_1\partial_2c)(\partial_2^2c)\bigr)\right|_{\mathfrak t}.
> $$
>
> Here $\mathfrak t$ is the Cartan subalgebra; the odd factors retain the displayed order.
>
> The cubic trace has coefficient $+1$ on basis-label triples $(1,1,2),(1,2,1),(2,1,1)$ and $-1$ on $(1,2,2),(2,1,2),(2,2,1)$; the other two coefficients vanish.
>
> It is therefore nonzero, while every tree boundary has zero Cartan restriction.
>
> The character dimension and exact rank bounds give
>
> $$
> 1\le\dim H(Q_0)^{(0,7)}\le968-669-298=1.
> $$
>
> The possible higher incoming source $(2,6)$ has tree cohomology dimension $45-45=0$.
>
> The complete scalar source of incoming $d_1$ has dimension $1002-837=165$.
>
> On its $165$ rational basis vectors, the target Cartan detector after the first-loop differential equals a rational functional composed with $Q_0$.
>
> The detector therefore vanishes on the $Q_1$ image of the full $Q_0$ kernel.
>
> The target cohomology is one-dimensional, so incoming $d_1=0$.
>
> $\square$

In the source convention, the secondary product is $\hbar[m_0-Q_1n_1]$.

Its $C_0^2$ indeterminacy image vanishes.

The possible source has weight $(3,3)$ and $N_c=3$; its $N_b=0$ tree cohomology vanishes, and its $N_b=1$ part satisfies

$$
0\le\dim E_2^{(1,3)}\le4-2-2=0.
$$

Multiplication by $C_0$ also annihilates the partner.

Both factors have zero $b$ count, so their quantum product equals their ordinary product.

The scalar complex at weight $(6,6)$ and $(N_b,N_c)=(0,9)$ satisfies

$$
0\le\dim H(Q_0)^{(0,9)}\le267-216-51=0.
$$

Every higher $Q_r$ annihilates its pure-ghost primitive.

Hence

$$
[C_0][m_0-Q_1n_1]=0.
$$

<a id="np-e2537f1f2c55de35-1"></a>

### $C_1^2$ and complete $C_3,C_4$ representatives

Represent $C_1$ by $a_1$ as above.

Its full product is

$$
a_1\star a_1=a_1^2+\hbar\mu_1(a_1,a_1)+\hbar^2\mu_2(a_1,a_1),
$$

The three components contain $822,8704,5216$ terms, respectively.

The complete closure equations are

$$
Q_0a_1^2=0,
\qquad Q_1a_1^2+Q_0\mu_1(a_1,a_1)=0,
$$

$$
Q_2a_1^2+Q_1\mu_1(a_1,a_1)+Q_0\mu_2(a_1,a_1)=0.
$$

The leading target satisfies $1\le\dim H(Q_0)^{(2,4)}(4,4)\le20-19=1$.

Cartan restriction gives the nonzero lower bound; incoming sources have minimum total weight $7+2r>8$.

Thus

$$
[C_1]^2\ne0.
$$

Tree cohomology also has the exact relation

$$
[a_1^2]= -\frac16[\operatorname{tr}(b^2)C_0^2].
$$

<a id="np-93544cd67af7d9c7-1"></a>

The full cocycles $C_3$ and $C_4$ each have three nonzero coefficients.

A representative is complete when its full differential vanishes at every loop order.

For the present $C_3,C_4$, the nonzero coefficients occur at orders $0,1,2$; higher equations still act on those coefficients.

$$
C_k=[\hbar^0]C_k+\hbar[\hbar]C_k+\hbar^2[\hbar^2]C_k,
\qquad k\in\{3,4\}.
$$

> **Proof of the finite completeness test.**
>
> The coefficient $[\hbar^s]C_k$ has $N_b=k-s$.
>
> Collecting powers of $\hbar$ gives
>
> $$
> \begin{aligned}
> [\hbar^0](Q_{\mathrm{pert}}C_k)
> &=Q_0([\hbar^0]C_k),\\
> [\hbar](Q_{\mathrm{pert}}C_k)
> &=Q_0([\hbar]C_k)+Q_1([\hbar^0]C_k),\\
> [\hbar^2](Q_{\mathrm{pert}}C_k)
> &=Q_0([\hbar^2]C_k)+Q_1([\hbar]C_k)+Q_2([\hbar^0]C_k),\\
> [\hbar^3](Q_{\mathrm{pert}}C_k)
> &=Q_1([\hbar^2]C_k)+Q_2([\hbar]C_k)+Q_3([\hbar^0]C_k).
> \end{aligned}
> $$
>
> For $C_4$, the additional allowed equation is
>
> $$
> [\hbar^4](Q_{\mathrm{pert}}C_4)
> =Q_2([\hbar^2]C_4)+Q_3([\hbar]C_4)+Q_4([\hbar^0]C_4)=0.
> $$
>
> The trace tables and rational certificates satisfy every displayed equation with zero right-hand side.
>
> Terms above order $k$ vanish because $Q_r$ lowers the $b$ count by $r$.
>
> Thus the three nonzero coefficients give full representatives.
>
> In the specified scalar blocks, the dimensions of the spaces of classes are
>
> $$
> \begin{array}{c|c}
> N_b+\text{loop order}=3&26-18-7=1\\
> N_b+\text{loop order}=4&129-97-31=1.
> \end{array}
> $$
>
> Each row subtracts the outgoing differential rank and incoming differential rank from the number of coefficients available in the block.
>
> Therefore a fixed leading class has a unique lift modulo full perturbative boundaries in that block.
>
> $\square$

Their explicit leading terms are

$$
\begin{aligned}
[\hbar^0]C_3={}&-\frac{111}{112}\operatorname{tr}(b^2)a_1
+\frac{57}{56}\operatorname{tr}(b^3)C_0\\
&+\frac{15}{14}\bigl[
\operatorname{tr}(b^2\partial_2c)\operatorname{tr}(b\partial_1c)
-\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b\partial_2c)\bigr],
\end{aligned}
$$

$$
[\hbar^0]C_4=\operatorname{tr}(b^3)a_1
-\frac32\operatorname{tr}(b^2\partial_1c\,b^2\partial_2c)
-\frac32\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c).
$$

The last coefficient of $C_3$ has a short trace expression:

$$
\begin{aligned}
[\hbar^2]C_3=\frac14\operatorname{tr}\Bigl\{b\bigl[&
(\partial_2^3c)(\partial_1^3c)
-3(\partial_1\partial_2^2c)(\partial_1^2\partial_2c)\\
&+3(\partial_1^2\partial_2c)(\partial_1\partial_2^2c)
-(\partial_1^3c)(\partial_2^3c)\bigr]\Bigr\}.
\end{aligned}
$$

The complete trace tables give the remaining coefficients.

The component term counts are

$$
C_3:(434,2434,96),\qquad C_4:(864,10728,14686).
$$

Coefficientwise evaluation of $Q_{\mathrm{pert}}C_3=Q_{\mathrm{pert}}C_4=0$ includes full $Q_0$, every permitted higher $Q_r$, gauge invariance, and both rotation actions.

The corresponding leading sectors have one-dimensional $E_\infty$.

In the fixed homogeneous loop blocks $N_b+\text{loop order}=3$ and $4$, the full scalar cohomology dimensions are

$$
26-18-7=1,
\qquad
129-97-31=1
$$

respectively.

The three numbers are the chain dimension, the outgoing full-differential rank, and the incoming full-differential rank.

Thus all scalar lifts of a given leading class differ by a full perturbative boundary within the specified block.

<a id="np-0dbe41c29df07e6c-1"></a>

### The second lift at $k=6$

The leading tree class is

$$
\left[3\operatorname{tr}(b^2)\operatorname{tr}(b^3)a_1
-\left(\frac{\operatorname{tr}(b^2)^3}{4}+\operatorname{tr}(b^3)^2\right)C_0\right].
$$

The first-correction space has dimension $118$, full-$Q_0$ rank $56$, and kernel dimension $62$.

The second obstruction for the initial fixed first correction has $1314$ Cartan terms.

The Cartan $Q_1$ image of the closed first-correction freedom has rank $8$ and cancels this obstruction.

The full joint solve retains $118+594=712$ correction coordinates and gives

$$
Q_0(\text{leading})=0,
$$

$$
Q_0(\text{first correction})+Q_1(\text{leading})=0,
$$

$$
Q_0(\text{second correction})+Q_1(\text{first correction})+Q_2(\text{leading})=0.
$$

The three components contain $3214,76234,355880$ terms, respectively.

The leading $E_2$ is one-dimensional, and the incoming source has minimum total weight $13+2r>14$.

Hence

$$
\dim_{\mathbb Q}E_3^{(6,2)}(7,7)_{\mathrm{scalar}}=1.
$$

This establishes $E_3$ survival; the third and higher corrections determine survival on later pages.
