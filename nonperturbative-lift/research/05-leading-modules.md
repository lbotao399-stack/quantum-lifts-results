<a id="np-section-03"></a>

## 5. Finite-rank leading classes: SU(2), SU(3), and SU(4)

Consider the leading sector

$$
(N_b,N_c)=(k,2),\qquad(n_1,n_2)=(k+1,k+1),\qquad k\ge0.
$$

Here “leading” specifies the coefficient before loop corrections.

A tree class is a solution of $Q_0\mathcal O=0$ modulo additions of $Q_0A$.

Multiplication by the invariant traces $\operatorname{tr}(b^2),\ldots,\operatorname{tr}(b^N)$ makes the leading classes into a module: its coefficients are polynomials in these traces.

“Free of finite rank” means that every class has a unique expression as such a polynomial combination of a fixed finite list of classes.

The ranks below are $1,3,6$ for $SU(2),SU(3),SU(4)$, respectively.

A complete tree module lists all generators and all relations in this leading scalar sector for every $k$.

The word “scalar” selects expressions unchanged by rotations of the two derivative directions.

Loop corrections are tested after this tree classification.

Both ghosts are first jets, and every $b$ is undifferentiated.

The source rotation generators $J_+,J_-$ act on nonnegative jet orders $m,n$ by

$$
J_+(\partial_1^m\partial_2^n b)=n\partial_1^{m+1}\partial_2^{n-1}b,
\qquad
J_-(\partial_1^m\partial_2^n b)=m\partial_1^{m-1}\partial_2^{n+1}b,
$$

and act on $c$ by the same rule.

A term with zero prefactor is defined to vanish.

The rotation-scalar projector for two first jets is $1-\frac12J_-J_+$.

<a id="np-6c4b9916fd1b9a28-1"></a>

### SU(2): first lifts in every degree

Cayley–Hamilton identity

$$
b^2=\frac12\operatorname{tr}(b^2)\mathbf1
$$

reduces each $b$ gap in a trace to length $0,1$.

Full-$Q_0$ primitives reduce the scalar projections of the four remaining trace seeds to the free module

$$
H(Q_0)_{\mathrm{leading\ scalar}}
=\mathbb Q[\operatorname{tr}(b^2)]\,[C_0].
$$

> **Example: what the rank-one module counts.**
>
> The generator $C_0$ has zero $b$-degree, and each factor $\operatorname{tr}(b^2)$ adds two.
>
> At degrees $k=0,2,4,6$, representatives are
>
> $$
> C_0,\qquad \operatorname{tr}(b^2)C_0,\qquad
> \operatorname{tr}(b^2)^2C_0,\qquad
> \operatorname{tr}(b^2)^3C_0.
> $$
>
> Every polynomial coefficient is a linear combination of these powers.
>
> There is one class in each even degree and zero classes in each odd degree.
>
> The source representative has $[\hbar^0]C_2=-\operatorname{tr}(b^2)C_0$, fixing its sign within this module.
>
> $\square$

**Proposition P4.**

For every $k\ge0$,

$$
\dim E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}1,&k\text{ is even},\\0,&k\text{ is odd}.\end{cases}
$$

> **Proof.**
>
> Cartan restriction is nonzero on each $\operatorname{tr}(b^2)^mC_0$ and annihilates every tree boundary, proving freeness.
>
> The first obstruction obeys two full identities:
>
> $$
> Q_1\bigl(\operatorname{tr}(b^2)C_0\bigr)
> =Q_0([\hbar]C_2),
> $$
>
> $$
> Q_1\bigl(\operatorname{tr}(b^2)^2\bigr)
> -2\operatorname{tr}(b^2)Q_1\operatorname{tr}(b^2)
> =16Q_0\left[
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}((\partial_1b)b\partial_2b)\right].
> $$
>
> Since $Q_1$ has differential order two, for integers $m\ge2$,
>
> $$
> \begin{aligned}
> Q_1\bigl(\operatorname{tr}(b^2)^mC_0\bigr)
> =Q_0\Bigg[&m\operatorname{tr}(b^2)^{m-1}[\hbar]C_2\\
> &+8m(m-1)\operatorname{tr}(b^2)^{m-2}C_0
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}((\partial_1b)b\partial_2b)\Bigg].
> \end{aligned}
> $$
>
> The cases $m=0,1$ follow from $Q_1C_0=0$ and the first identity, respectively.
>
> An incoming $d_1$ source would contain $k+1$ copies of $b$ and one differentiated ghost, with minimum total weight $2k+3>2k+2$.
>
> All leading classes therefore survive the first differential.
>
> $\square$

This proves the all-degree $E_2$ statement.

Higher-loop survival is determined by the subsequent closure equations.

<a id="np-96250f1e20680d7c-1"></a>

### SU(3): tree module and radial equations

The first class uses the source representative

$$
a_1=-\frac12\operatorname{tr}(b\partial_1c\partial_2c-b\partial_2c\partial_1c).
$$

It represents $C_1$ and satisfies $Q_0a_1=Q_1a_1=0$.

The Cayley–Hamilton identity is

$$
b^3=\frac12\operatorname{tr}(b^2)b+\frac13\operatorname{tr}(b^3)\mathbf1.
$$

Full tree primitives reduce the nine trace seeds to three free generators over the Casimir ring $\mathbb Q[\operatorname{tr}(b^2),\operatorname{tr}(b^3)]$:

$$
C_0,\qquad a_1,\qquad
\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
$$

With $t$ counting $b$-degree, the Hilbert series and degree dimensions are

$$
\frac{1+t+t^2}{(1-t^2)(1-t^3)}
=\frac1{(1-t)(1-t^2)},
\qquad
\dim H(Q_0)^{(k,2)}_{\mathrm{scalar}}
=\left\lfloor\frac k2\right\rfloor+1.
$$

> **Example: the two tree classes at degree two.**
>
> The degrees of the three generators are $0,1,2$.
>
> At $k=2$, the degree-zero generator takes the coefficient $\operatorname{tr}(b^2)$, and the degree-two generator takes a constant.
>
> The degree-one generator would require a degree-one trace, which vanishes for traceless matrices.
>
> Hence a basis is
>
> $$
> [\operatorname{tr}(b^2)C_0],\qquad
> [\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)],
> \qquad \dim H(Q_0)^{(2,2)}_{\mathrm{scalar}}=2.
> $$
>
> The radial equations below set the double-trace coefficient to zero.
>
> The remaining coefficient of $C_0$ is a multiple of $\operatorname{tr}(b^2)$, and its derivative must vanish:
>
> $$
> \frac{\partial(\text{coefficient of }C_0)}{\partial\operatorname{tr}(b^2)}=0,
> \qquad
> \dim E_2^{(2,2)}(3,3)_{\mathrm{scalar}}=0.
> $$
>
> This is an explicit degree where the complete tree module has two classes and the first loop removes both.
>
> $\square$

> **Proof of freeness.**
>
> Trace reduction establishes spanning in the full complex.
>
> In the source trace-dual Cartan coordinates $b_1,b_2$, the coefficient determinant of the three symmetric ghost tensors is
>
> $$
> \det\begin{pmatrix}
> 2&-(b_1+2b_2)/3&b_1^2\\
> -1&(-b_1+b_2)/3&b_1b_2\\
> 2&(2b_1+b_2)/3&b_2^2
> \end{pmatrix}
> =-3b_1b_2(b_1+b_2)\ne0.
> $$
>
> Cartan restriction annihilates every tree boundary, so the three generators are independent over the polynomial ring.
>
> $\square$

The Cartan $Q_1$ kernel forces the coefficient of the third generator to vanish.

The remaining two polynomial coefficients satisfy

$$
\frac{\partial(\text{coefficient of }C_0)}{\partial\operatorname{tr}(b^2)}
+\frac{\operatorname{tr}(b^2)}4
\frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)}=0,
$$

$$
\frac{\partial(\text{coefficient of }C_0)}{\partial\operatorname{tr}(b^3)}
+\frac23
\frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)}=0.
$$

These are ordinary derivatives with respect to the two independent Casimir polynomial coordinates.

> **Proof of the all-degree kernel.**
>
> The coefficient minor of two target tensors is $-216b_1\ne0$, giving both equations coefficient by coefficient.
>
> Equality of mixed derivatives gives
>
> $$
> \frac{\partial^2(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)^2}
> =\frac{3\operatorname{tr}(b^2)}8
> \frac{\partial^2(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)^2}.
> $$
>
> Expand in Casimir monomials of weighted degree $k-1$, with exponents $r,s\ge0$ satisfying $2r+3s=k-1$.
>
> Adjacent coefficients obey the recurrence
>
> $$
> \begin{aligned}
> &r(r-1)\,[\text{coefficient of }\operatorname{tr}(b^2)^r\operatorname{tr}(b^3)^s]\\
> &\qquad=\frac38(s+2)(s+1)
> [\text{coefficient of }\operatorname{tr}(b^2)^{r-3}\operatorname{tr}(b^3)^{s+2}].
> \end{aligned}
> $$
>
> Coefficients with negative exponents vanish.
>
> For $k\equiv2\pmod3$, the smallest possible $r$ is $2$, and the endpoint equation annihilates the entire recurrence chain.
>
> At $k=2$, the weighted-degree-one monomial space is empty.
>
> In the remaining degrees, the smallest $r$ is $0$ or $1$.
>
> One coefficient is free, and the recurrence determines all others uniquely.
>
> For $k>0$, the weighted Euler identity gives the other coordinate:
>
> $$
> \text{coefficient of }C_0
> =-\frac1k\left[
> \frac{\operatorname{tr}(b^2)^2}{2}
> \frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)}
> +2\operatorname{tr}(b^3)
> \frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)}
> \right].
> $$
>
> At $k=0$, the solution is $C_0$.
>
> $\square$

Lifting the Cartan solutions to full $E_2$ uses five full-$Q_0$ primitives: three pure-Casimir second-order brackets, together with

$$
Q_1\bigl(\operatorname{tr}(b^2)a_1\bigr)
-\frac23Q_1\bigl(\operatorname{tr}(b^3)C_0\bigr),
$$

$$
Q_1\bigl(\operatorname{tr}(b^3)a_1\bigr)
-\frac{\operatorname{tr}(b^2)}4
Q_1\bigl(\operatorname{tr}(b^2)C_0\bigr).
$$

The exact primitives of the latter two expressions contain $2546,10744$ terms, respectively.

The second-order chain rule expresses the full $Q_1$ image of every polynomial solution as the $Q_0$ image of a combination of these five primitives.

Repeated-coordinate second derivatives carry $1/2$, and mixed derivatives carry $1$.

The minimum-weight bound excludes incoming $d_1$, yielding

$$
\dim_{\mathbb Q}E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}
0,&k\equiv2\pmod3,\\
1,&k\not\equiv2\pmod3.
\end{cases}
$$

<a id="np-7d31a48e28a59909-1"></a>

### SU(4): complete tree module and Cartan classification

The Cayley–Hamilton identity becomes

$$
b^4=\frac12\operatorname{tr}(b^2)b^2+\frac13\operatorname{tr}(b^3)b
+\left[\frac14\operatorname{tr}(b^4)-\frac18\operatorname{tr}(b^2)^2\right]\mathbf1.
$$

Full-$Q_0$ reductions of ten single-trace and six double-trace seeds give six free generators over $\mathbb Q[\operatorname{tr}(b^2),\operatorname{tr}(b^3),\operatorname{tr}(b^4)]$:

$$
\begin{gathered}
C_0,\quad a_1,\quad
\operatorname{tr}(b^2\partial_1c\partial_2c),\quad
\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c),\\
\operatorname{tr}(b^3\partial_1c\partial_2c),\quad
\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c),
\end{gathered}
$$

Each expression is projected to its scalar part as above.

Their degrees are $0,1,2,2,3,4$.

At the source Cartan point $(b_1,b_2,b_3)=(1,2,4)$, the determinant of the six symmetric ghost tensors is $1016064\ne0$, proving independence over the polynomial ring.

Consequently,

$$
\frac{1+t+2t^2+t^3+t^4}{(1-t^2)(1-t^3)(1-t^4)}
=\frac1{(1-t)(1-t^2)^2},
$$

$$
\dim H(Q_0)^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\binom{\lfloor k/2\rfloor+2}{2}.
$$

> **Example: the complete SU(4) tree module at degree two.**
>
> The six generator degrees are $0,1,2,2,3,4$.
>
> At $k=2$, the degree-zero generator contributes once through $\operatorname{tr}(b^2)$, and each degree-two generator contributes once through a constant.
>
> Thus
>
> $$
> \begin{gathered}
> [\operatorname{tr}(b^2)C_0],\qquad
> \left[\left(1-\frac12J_-J_+\right)
> \operatorname{tr}(b^2\partial_1c\partial_2c)\right],\\
> [\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)],
> \qquad
> \dim H(Q_0)^{(2,2)}(3,3)_{\mathrm{scalar}}=1+1+1=3.
> \end{gathered}
> $$
>
> The Cartan test leaves the following one-dimensional combination:
>
> $$
> \frac12\operatorname{tr}(b^2)C_0
> -4\left(1-\frac12J_-J_+\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
> +\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
> $$
>
> The first correction displayed below solves the full equation for this combination.
>
> Its second correction vanishes, and the remaining allowed equations also vanish.
>
> This establishes a complete perturbative representative at $k=2$.
>
> $\square$

The Cartan $Q_1$ equations reduce the six polynomial coefficients to three.

The degree-two combination is

$$
\left(1-\frac12J_-J_+\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
-\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
$$

Expand its polynomial coefficient in characteristic-polynomial coordinates:

$$
\left(-\frac{\operatorname{tr}(b^2)}2\right)^r
\left(-\frac{\operatorname{tr}(b^3)}3\right)^s
\left(\frac{\operatorname{tr}(b^2)^2}{8}-\frac{\operatorname{tr}(b^4)}4\right)^t,
\qquad2r+3s+4t=k-2.
$$

Here $r,s,t$ are nonnegative monomial exponents.

> **Proof of the Cartan classification.**
>
> Comparison of coefficients gives
>
> $$
> 4(r+2)(r+1)[\text{coefficient at }(r+2,s,t)]
> =-(t+1)(2r+s+1)[\text{coefficient at }(r,s,t+1)],
> $$
>
> $$
> (r+1)(t+1)[\text{coefficient at }(r+1,s,t+1)]
> =(s+2)(s+1)[\text{coefficient at }(r,s+2,t)],
> $$
>
> $$
> \begin{aligned}
> 0={}&4(r+1)(s+1)[\text{coefficient at }(r+1,s+1,t)]\\
> &+2(s+1)(t+1)[\text{coefficient at }(r-1,s+1,t+1)]\\
> &+(t+2)(t+1)[\text{coefficient at }(r,s-1,t+2)].
> \end{aligned}
> $$
>
> The first equation connects all $(r,t)$ at fixed $s$; the second connects chains whose $s$ values differ by $2$.
>
> The solution space therefore has dimension at most one.
>
> At $k-2=5+4t$, the endpoint of the third equation is $4[\text{coefficient at }(1,1,t)]=0$, annihilating the entire chain.
>
> At $k-2=1$, the monomial space is empty.
>
> The characteristic residue in the next section supplies nonzero solutions in the remaining degrees.
>
> Hence
>
> $$
> \dim\ker(\mathrm{Cartan}\circ Q_1)\big|_{H(Q_0)^{(k,2)}_{\mathrm{scalar}}}
> =\begin{cases}0,&k\equiv3\pmod4,\\1,&k\not\equiv3\pmod4.\end{cases}
> $$
>
> $\square$

In sectors $k\equiv3\pmod4$, the full $E_2$ vanishes.

In the remaining sectors, the full $E_2$ dimension is at most one.

Full-$Q_0$ primitives have been constructed at $k=2,4,5$.

The $k=6$ primitive is the remaining finite step in upgrading the Cartan theorem to an all-degree full-$E_2$ theorem.

At $k=2$, a complete cocycle has leading term

$$
\frac12\operatorname{tr}(b^2)C_0
-4\left(1-\frac12J_-J_+\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
+\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c),
$$

and first correction

$$
\begin{aligned}
\left(1-\frac{J_-J_+}{6}\right)\left(1-\frac{J_-J_+}{2}\right)\Bigg[&
-\frac56\operatorname{tr}(b\partial_2^2c\partial_1^2c)\\
&+\frac76\operatorname{tr}((\partial_2b)(\partial_2c)(\partial_1^2c))\\
&+\frac76\operatorname{tr}((\partial_2b)(\partial_1^2c)(\partial_2c))\\
&+\frac{23}{6}\operatorname{tr}((\partial_2^2b)(\partial_1c)(\partial_1c))\Bigg].
\end{aligned}
$$

The second correction vanishes.

The full closure equations are

$$
Q_0(\text{first correction})=-Q_1(\text{leading term}),
\qquad
Q_1(\text{first correction})=-Q_2(\text{leading term}).
$$

The $b$ count makes higher equations vanish.

The nonzero leading Cartan class and incoming weight bound establish full perturbative survival.

The first lifts at $k=4,5$ satisfy full rational polynomial identities; their higher-loop closure equations determine subsequent survival.
