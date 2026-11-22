# SU(4): complete leading tree module and all-degree Cartan kernel

All source traces, fields and loop normalization are unchanged. The sector has $N_b=k$, $N_c=2$ and weight $(k+1,k+1)$ in the source gauge-invariant rotation-scalar complex with positive-jet ghosts. Here $k$ is a nonnegative integer. Its leading terms contain only bare $b$ and first jets of $c$.

The source generators are

$$
C_0=\operatorname{tr}(\partial_1c\partial_2c),\qquad
a_1=-\frac12\operatorname{tr}(b\partial_1c\partial_2c-b\partial_2c\partial_1c).
$$

**Proposition 1.** The leading tree module is free over $\mathbb Q[\operatorname{tr}(b^2),\operatorname{tr}(b^3),\operatorname{tr}(b^4)]$, with the six scalar-projected generators

$$
C_0,\quad a_1,\quad
\operatorname{tr}(b^2\partial_1c\partial_2c),\quad
\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c),\quad
\operatorname{tr}(b^3\partial_1c\partial_2c),\quad
\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c).
$$

> **Proof.** Invariant contractions are products of matrix traces. With two ghosts, either both occupy one trace or each occupies a separate trace. All remaining traces contain only $b$. The source Cayley--Hamilton identity is
>
> $$
> b^4=\frac12\operatorname{tr}(b^2)b^2+\frac13\operatorname{tr}(b^3)b+
> \left[\frac14\operatorname{tr}(b^4)-\frac18\operatorname{tr}(b^2)^2\right]\mathbf1.
> $$
>
> Every bare-$b$ gap therefore has length $0,1,2,3$. Graded cyclicity and scalar projection identify the two gap orders. Tracelessness removes a separate trace with one ghost and no $b$. There remain ten single-trace and six double-trace seeds.
>
> In the following table, apply the source scalar projector $1-J_-J_+/2$ to every entry; $J_\pm$ are the source rotation generators on the two jet directions. Each seed minus its displayed representative equals full $Q_0$ of the displayed primitive.
>
> | Seed | Representative modulo $Q_0$ | Primitive |
> |:---|:---|:---|
> | $\operatorname{tr}(\partial_1c\partial_2c)$ | $C_0$ | $0$ |
> | $\operatorname{tr}(\partial_1cb\partial_2c)$ | $-a_1$ | $0$ |
> | $\operatorname{tr}(\partial_1cb^{2}\partial_2c)$ | $\operatorname{tr}(b^2\partial_1c\partial_2c)$ | $0$ |
> | $\operatorname{tr}(\partial_1cb^{3}\partial_2c)$ | $\operatorname{tr}(b^3\partial_1c\partial_2c)$ | $0$ |
> | $\operatorname{tr}(b\partial_1cb\partial_2c)$ | $\operatorname{tr}(b^2\partial_1c\partial_2c)$ | $\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1cb^{2}\partial_2c)$ | $\operatorname{tr}(b^3\partial_1c\partial_2c)$ | $\operatorname{tr}((\partial_2b)b^{2}\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1cb^{3}\partial_2c)$ | $\frac{1}{4}\operatorname{tr}(b^4)C_0-\frac{1}{3}\operatorname{tr}(b^3)a_1+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^2\partial_1c\partial_2c)-\frac{1}{8}\operatorname{tr}(b^2)^{2}C_0$ | $\operatorname{tr}((\partial_2b)b^{3}\partial_1c)$ |
> | $\operatorname{tr}(b^{2}\partial_1cb^{2}\partial_2c)$ | $\frac{1}{4}\operatorname{tr}(b^4)C_0-\frac{1}{3}\operatorname{tr}(b^3)a_1+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^2\partial_1c\partial_2c)-\frac{1}{8}\operatorname{tr}(b^2)^{2}C_0$ | $\operatorname{tr}(b(\partial_2b)b^{2}\partial_1c)+\operatorname{tr}((\partial_2b)b^{3}\partial_1c)$ |
> | $\operatorname{tr}(b^{2}\partial_1cb^{3}\partial_2c)$ | $-\frac{1}{4}\operatorname{tr}(b^4)a_1+\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}(b^2\partial_1c\partial_2c)+\frac{1}{8}\operatorname{tr}(b^2)^{2}a_1+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^3\partial_1c\partial_2c)$ | $\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b^{2}\partial_1c)+\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}((\partial_2b)b\partial_1c)+\operatorname{tr}(b(\partial_2b)b^{3}\partial_1c)$ |
> | $\operatorname{tr}(b^{3}\partial_1cb^{3}\partial_2c)$ | $-\frac{1}{6}\operatorname{tr}(b^2)\operatorname{tr}(b^3)a_1+\frac{1}{4}\operatorname{tr}(b^4)\operatorname{tr}(b^2\partial_1c\partial_2c)+\frac{1}{8}\operatorname{tr}(b^2)\operatorname{tr}(b^4)C_0+\frac{1}{8}\operatorname{tr}(b^2)^{2}\operatorname{tr}(b^2\partial_1c\partial_2c)+\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}(b^3\partial_1c\partial_2c)-\frac{1}{16}\operatorname{tr}(b^2)^{3}C_0$ | $\operatorname{tr}(b^{2}(\partial_2b)b^{3}\partial_1c)+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b^{3}\partial_1c)+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b(\partial_2b)b^{2}\partial_1c)+\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}((\partial_2b)b^{2}\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $0$ |
> | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b^{2}\partial_2c)$ | $-\frac{1}{6}\operatorname{tr}(b^3)C_0+\frac{1}{2}\operatorname{tr}(b^2)a_1+2\operatorname{tr}(b^3\partial_1c\partial_2c)$ | $\operatorname{tr}((\partial_2b)b^{2}\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b^{3}\partial_2c)$ | $\frac{1}{2}\operatorname{tr}(b^4)C_0-\frac{1}{2}\operatorname{tr}(b^3)a_1+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^2\partial_1c\partial_2c)+\frac{1}{4}\operatorname{tr}(b^2)\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)-\frac{1}{4}\operatorname{tr}(b^2)^{2}C_0-\frac{1}{2}\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c)$ | $\frac{3}{2}\operatorname{tr}((\partial_2b)b^{3}\partial_1c)+\frac{1}{2}\operatorname{tr}(b(\partial_2b)b^{2}\partial_1c)-\frac{1}{4}\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b^{2}\partial_1c)\operatorname{tr}(b^{2}\partial_2c)$ | $\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c)$ | $0$ |
> | $\operatorname{tr}(b^{2}\partial_1c)\operatorname{tr}(b^{3}\partial_2c)$ | $-\frac{1}{2}\operatorname{tr}(b^4)a_1+\frac{1}{2}\operatorname{tr}(b^3)\operatorname{tr}(b^2\partial_1c\partial_2c)-\frac{1}{6}\operatorname{tr}(b^3)\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)+\frac{1}{4}\operatorname{tr}(b^2)^{2}a_1+\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^3\partial_1c\partial_2c)$ | $\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b^{2}\partial_1c)+\operatorname{tr}(b(\partial_2b)b^{3}\partial_1c)+\frac{1}{2}\operatorname{tr}(b^3)\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b^{3}\partial_1c)\operatorname{tr}(b^{3}\partial_2c)$ | $\frac{1}{9}\operatorname{tr}(b^3)^{2}C_0-\frac{2}{3}\operatorname{tr}(b^2)\operatorname{tr}(b^3)a_1+\operatorname{tr}(b^4)\operatorname{tr}(b^2\partial_1c\partial_2c)-\frac{1}{4}\operatorname{tr}(b^4)\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)+\frac{1}{4}\operatorname{tr}(b^2)\operatorname{tr}(b^4)C_0-\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}(b^3\partial_1c\partial_2c)-\frac{1}{2}\operatorname{tr}(b^2)\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c)-\frac{1}{8}\operatorname{tr}(b^2)^{3}C_0+\frac{1}{8}\operatorname{tr}(b^2)^{2}\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\operatorname{tr}(b^{2}(\partial_2b)b^{3}\partial_1c)+\frac{1}{4}\operatorname{tr}(b^4)\operatorname{tr}((\partial_2b)b\partial_1c)+\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b^{3}\partial_1c)+\frac{1}{3}\operatorname{tr}(b^3)\operatorname{tr}((\partial_2b)b^{2}\partial_1c)+\operatorname{tr}(b^2)\operatorname{tr}(b(\partial_2b)b^{2}\partial_1c)-\frac{1}{8}\operatorname{tr}(b^2)^{2}\operatorname{tr}((\partial_2b)b\partial_1c)$ |

>
> These sixteen identities were checked with bare ghosts retained. Multiplication by each Casimir commutes with $Q_0$, so they prove spanning in every degree. On the source Cartan coordinates $(b_1,b_2,b_3)=(1,2,4)$, the six symmetric ghost coefficient tensors have determinant $1016064$. Cartan restriction annihilates every $Q_0$ boundary; this nonzero determinant proves independence over the Casimir polynomial ring.
>
> Let $t$ count bare-$b$ degree. The six generator degrees are $0,1,2,2,3,4$, hence
>
> $$
> \frac{1+t+2t^2+t^3+t^4}{(1-t^2)(1-t^3)(1-t^4)}
> =\frac1{(1-t)(1-t^2)^2},
> $$
>
> $$
> \dim H(Q_0)^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
> =\binom{\lfloor k/2\rfloor+2}{2}.
> $$
>
> Independently, the symmetric square of the standard $S_4$ representation has characters $6,2,2,0,0$ on the identity, transposition, double-transposition, three-cycle and four-cycle classes. Their sizes are $1,6,3,8,6$. Its Molien series is
>
> $$
> \frac1{24}\left[
> \frac6{(1-t)^3}+\frac{12}{(1-t)(1-t^2)}+
> \frac6{(1-t)(1+t)^2}\right]
> =\frac1{(1-t)(1-t^2)^2}.
> $$
>
> $\square$

**Proposition 2.** The Cartan restriction of the one-loop map has all-degree kernel

$$
\dim\ker(\mathrm{Cartan}\circ Q_1)
\big|_{H(Q_0)^{(k,2)}_{\mathrm{scalar}}}
=\begin{cases}0,&k\equiv3\pmod4,\\1,&k\not\equiv3\pmod4.\end{cases}
$$

> **Proof.** The six zeroth images and eighteen first-Casimir images are saved in `leading_radial.pkl`. The second-Casimir-derivative terms vanish on Cartan by the full-$Q_0$-exact universal Casimir bracket. Thus these twenty-four evaluations determine the entire polynomial map.
>
> Its second-jet-ghost sector has a $3\times3$ coefficient minor $-169728$ at the same source Cartan point. It forces the degree-three and degree-four generator coefficients to vanish, and combines the two degree-two generators only as
>
> $$
> \left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
> -\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
> $$
>
> Call its polynomial multiplier the quadratic-trace coefficient. The remaining six tensor operators are independent; `radial_equations.pkl/json` records their nonzero coefficient minor and all exact reduction identities. Comparing mixed derivatives of the two other coefficients gives three recurrence equations for the quadratic-trace coefficient.
>
> Expand that coefficient in the three actual coefficients of the characteristic polynomial:
>
> $$
> \left(-\frac{\operatorname{tr}(b^2)}2\right)^r
> \left(-\frac{\operatorname{tr}(b^3)}3\right)^s
> \left(\frac{\operatorname{tr}(b^2)^2}{8}-\frac{\operatorname{tr}(b^4)}4\right)^t,
> \qquad 2r+3s+4t=k-2.
> $$
>
> Here $r,s,t$ are nonnegative monomial exponents. A coefficient with a negative exponent is zero. Direct comparison gives
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
> At fixed $s$, the first equation connects every allowed exponent pair $(r,t)$ and leaves at most one coefficient. Allowed $s$ values differ by two. The second equation connects each consecutive pair of those chains, with nonzero coefficients, so the entire kernel has dimension at most one. For $k-2=5+4t$, the third equation at $(r,s)=(0,0)$ reads
>
> $$
> 4[\text{coefficient at }(1,1,t)]=0.
> $$
>
> This kills the minimal-$s$ chain and therefore every chain. At $k-2=1$, no monomial exists. These are precisely $k\equiv3\pmod4$.
>
> For existence in all other degrees, $z$ is the auxiliary Laurent variable and $[z^{-1}]$ extracts its coefficient at infinity. The original characteristic polynomial is
>
> $$
> \det(z\mathbf1-b)=z^4-\frac12\operatorname{tr}(b^2)z^2
> -\frac13\operatorname{tr}(b^3)z
> +\frac18\operatorname{tr}(b^2)^2-\frac14\operatorname{tr}(b^4).
> $$
>
> A representative is
>
> $$
> \begin{aligned}
> [z^{-1}]\det(z\mathbf1-b)^{(k-3)/4}\Bigg[&
> \left(z^2-\frac14\operatorname{tr}(b^2)\right)C_0-z a_1\\
> &+\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
> -\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)\Bigg].
> \end{aligned}
> $$
>
> The exact finite operator identity in `period_identity.pkl` makes its Cartan $Q_1$ image a total Laurent derivative. At the recurrence level, the first and third displayed equations follow respectively from
>
> $$
> \frac{k-3}{4}[z^{-1}]\frac{d}{dz}
> \left[z\det(z\mathbf1-b)^{(k-7)/4}\right]=0,
> $$
>
> $$
> \frac{k-3}{4}[z^{-1}]\frac{d}{dz}
> \det(z\mathbf1-b)^{(k-7)/4}=0;
> $$
>
> the second follows from $z^2\cdot1=z\cdot z$ in the differentiated determinant coefficients. For $k\ge2$ outside the excluded congruence class, the exponent $(k-3)/4$ is nonintegral and at least one monomial of weight $k-2$ exists. Its generalized multinomial coefficient is nonzero. The two remaining polynomial coefficients are uniquely reconstructed from their gradients by weighted Euler division by $k$ and $k-1$. At $k=0,1$, the representatives are $C_0$ and $-a_1$.
>
> The recurrence and total-derivative arguments prove all degrees; the saved finite checks are readbacks. $\square$

An incoming $d_1$ is absent by weight. Consequently this theorem gives zero $E_2$ in every $k\equiv3\pmod4$ sector and an upper bound of one in every other sector. Full $Q_0$ primitives are still required to turn all the latter candidates into an $E_2$ theorem.

**Proposition 3.** At $k=2$, the explicit full perturbative cocycle has leading term

$$
\frac12\operatorname{tr}(b^2)C_0
-4\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
+\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
$$

> **Proof.** The first correction is the following literal trace polynomial, with the source rotation generators acting on its four derivative slots:
>
> $$
> \begin{aligned}
> \left(1-\frac{J_-J_+}{6}\right)\left(1-\frac{J_-J_+}{2}\right)\Bigg[&
> -\frac56\operatorname{tr}(b\partial_2^2c\partial_1^2c)\\
> &+\frac76\operatorname{tr}((\partial_2b)(\partial_2c)(\partial_1^2c))\\
> &+\frac76\operatorname{tr}((\partial_2b)(\partial_1^2c)(\partial_2c))\\
> &+\frac{23}6\operatorname{tr}((\partial_2^2b)(\partial_1c)(\partial_1c))\Bigg].
> \end{aligned}
> $$
>
> `k2_full.pkl` stores its full rational reconstruction. The leading and first correction have $721$ and $960$ component monomials; the second correction is zero. Direct evaluation gives
>
> $$
> Q_0(\text{first correction})=-Q_1(\text{leading term}),\qquad
> Q_1(\text{first correction})=-Q_2(\text{leading term}).
> $$
>
> Higher terms vanish by the fan input-$b$ count. The correction domain has dimensions $4+1$ and the joint map has rank $4$. Both rotation equations pass. The Cartan leading class is nonzero. An incoming differential of order $r\ge1$ would require $N_b=2+r,N_c=1$, of minimum total weight $5+2r>6$, so the class survives. A fresh process reconstructed all five trace sources and recomputed every equation with the original Wick evaluator, skipping only channels excluded by the exact input-$b$ count. All checks passed in $349.23$ seconds; omitting the first correction trace produced a $2016$-term nonzero residual (`k2_independent_readback.json`). $\square$

At $k=4$, the leading expression

$$
\begin{aligned}
&\operatorname{tr}(b^2)\left[
\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
-\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)\right]\\
&+\left[\frac12\operatorname{tr}(b^4)-\frac5{16}\operatorname{tr}(b^2)^2\right]C_0
-\frac23\operatorname{tr}(b^3)a_1
\end{aligned}
$$

has $7712$ terms. Its $162588$-term full $Q_1$ image equals full $Q_0$ of an $81401$-term rational primitive selected from $35$ independent trace images (`k4_first/first_lift.pkl/json`). Thus its scalar $E_2$ dimension is one. Higher-loop survival is open.

No physical instanton operator or nonperturbative coefficient is inferred from these perturbative statements.
