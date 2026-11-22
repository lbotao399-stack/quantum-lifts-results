# SU(4): finite gates for a universal first lift

Status at 07:13: the $k=2$ full cocycle and the $k=4,5$ positive first-lift primitives are certified and independently reconstructed. The $k=6$ primitive remains under construction. The following construction remains conditional on that final full rational identity and independent readback. Its algebraic signs have been independently reviewed by the source lane, and `verify_first_lift_chain_rule.py/json` independently checks the complete formal chain rule, all three first-derivative reductions, and all six remaining coefficients against `radial_equations.json`.

The quadratic-trace combination means exactly

$$
\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
-\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c).
$$

The four required leading expressions, in order $k=2,4,5,6$, are:

$$
\frac12\operatorname{tr}(b^2)C_0-4(\text{quadratic-trace combination}),
$$

$$
\operatorname{tr}(b^2)(\text{quadratic-trace combination})+
\left[\frac12\operatorname{tr}(b^4)-\frac5{16}\operatorname{tr}(b^2)^2\right]C_0
-\frac23\operatorname{tr}(b^3)a_1,
$$

$$
\operatorname{tr}(b^3)(\text{quadratic-trace combination})+
\left[\frac3{16}\operatorname{tr}(b^2)^2-\frac34\operatorname{tr}(b^4)\right]a_1,
$$

$$
\begin{aligned}
&\left[3\operatorname{tr}(b^2)^2-8\operatorname{tr}(b^4)\right](\text{quadratic-trace combination})\\
&+\left[\operatorname{tr}(b^2)\operatorname{tr}(b^4)
-\frac{11}{24}\operatorname{tr}(b^2)^3-\frac49\operatorname{tr}(b^3)^2\right]C_0
+\frac43\operatorname{tr}(b^2)\operatorname{tr}(b^3)a_1.
\end{aligned}
$$

Here a **positive primitive** satisfies $Q_0(\text{primitive})=Q_1(\text{leading term})$; the actual first correction is its negative. The saved $k=2$ artifact instead records the actual first correction, so its sign is retained separately.

The source zeroth images of $C_0$ and $a_1$ under $Q_1$ vanish; the quadratic-trace combination has the nonzero zeroth image retained below. These are literal source-operator statements.

> **Conditional proof.** Divide the actual first correction at $k=2$ by four. Its full $Q_0$ image is
>
> $$
> Q_1(\text{quadratic-trace combination})
> -\frac18 Q_1(\operatorname{tr}(b^2)C_0).
> $$
>
> Every pure-Casimir bracket has the already proved full-$Q_0$ primitive
>
> $$
> 4mn\left(1-\frac{J_-J_+}{2}\right)
> \operatorname{tr}\bigl((\partial_1b^{m-1})b(\partial_2b^{n-1})\bigr).
> $$
>
> This means its $Q_0$ image is $Q_1(\operatorname{tr}(b^m)\operatorname{tr}(b^n))-\operatorname{tr}(b^m)Q_1\operatorname{tr}(b^n)-\operatorname{tr}(b^n)Q_1\operatorname{tr}(b^m)$, in the unchanged source convention. The integer labels $m,n$ specify the two Casimir powers.
>
> For the first radial relation, take the positive $k=4$ primitive, subtract $\operatorname{tr}(b^2)$ times one quarter of the actual $k=2$ correction, and add $5/16$ times the $(2,2)$ Casimir-bracket primitive multiplied by $C_0$.
>
> For the second radial relation, take the positive $k=5$ primitive, subtract $\operatorname{tr}(b^3)$ times one quarter of the actual $k=2$ correction, and subtract $3/16$ times the $(2,2)$ Casimir-bracket primitive multiplied by $a_1$.
>
> For the third radial relation, add the following seven contributions, subtract the positive $k=6$ primitive, and divide by eight:
>
> | Multiplier | Primitive or source factor |
> |:---|:---|
> | $3\operatorname{tr}(b^2)^2-8\operatorname{tr}(b^4)$ | one quarter of the actual $k=2$ correction |
> | $6\operatorname{tr}(b^2)$ | first radial-relation primitive just constructed |
> | $3$ | $(2,2)$ Casimir-bracket primitive times the quadratic-trace combination |
> | $-11\operatorname{tr}(b^2)/8$ | $(2,2)$ Casimir-bracket primitive times $C_0$ |
> | $1$ | $(2,4)$ Casimir-bracket primitive times $C_0$ |
> | $-4/9$ | $(3,3)$ Casimir-bracket primitive times $C_0$ |
> | $4/3$ | $(2,3)$ Casimir-bracket primitive times $a_1$ |
>
> Applying full $Q_0$ to these three constructions gives exactly the three derivative relations in `radial_equations.pkl/json`. All multipliers are $Q_0$-closed. The universal positive primitive for a leading Cartan-kernel polynomial is then the sum of:
>
> 1. Its quadratic-trace coefficient times one quarter of the actual $k=2$ correction.
> 2. Each first Casimir derivative of that coefficient times the corresponding radial-relation primitive.
> 3. Each second Casimir derivative of each of the three polynomial coefficients, multiplied by its source generator and the corresponding Casimir-bracket primitive. The diagonal terms have factor $1/2$; each unordered mixed pair has factor $1$.
>
> The full order-two chain rule shows that the remaining $Q_1$ terms are precisely the six gradient equations already proved to characterize the Cartan kernel. They vanish. Thus the constructed expression has full $Q_0$ image equal to the original full $Q_1$ image, in every degree.
>
> Incoming $d_1$ is absent by the source weights. Once the remaining rational primitive and independent readback are certified, this construction proves $E_2$ dimension one for $k\not\equiv3\pmod4$ and zero otherwise. It says nothing by itself about higher-loop survival. $\square$

The conditional gate is deliberately retained until `k6_first/first_lift.json` and its independent saved-state readback are present and pass. The corresponding k5 files have passed.
