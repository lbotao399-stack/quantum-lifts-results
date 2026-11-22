# SU(2) leading scalar cohomology through the first differential

Use exactly the product-engine fields and $C_0=\operatorname{tr}(\partial_1c\partial_2c)$ of `../rational_c4/full_hidden_product_exact.pkl`. In that source convention, the leading component of $C_2$ is $-\operatorname{tr}(b^2)C_0$, and $[\hbar]C_2$ denotes its stored first-loop component. This note does not change either normalization.

**Proposition.** In the gauge-invariant rotation-scalar sector with $N_b=k$, $N_c=2$ and weight $(k+1,k+1)$,

$$
\dim H(Q_0)^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\dim E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
=\begin{cases}1,&k\text{ even},\\0,&k\text{ odd}.\end{cases}
$$

> **Proof.** All $b$ fields are undifferentiated and both ghosts are first jets. The $SU(2)$ invariant tensors are matrix-trace contractions; products of traces therefore span the invariant polynomials. Cayley–Hamilton gives
>
> $$
> b^2=\frac12\operatorname{tr}(b^2)\mathbf1.
> $$
>
> Every $b$ gap reduces to length zero or one. With two ghosts, the remaining trace seeds are $\operatorname{tr}(\partial_1c\partial_2c)$, $\operatorname{tr}(\partial_1c\,b\partial_2c)$, $\operatorname{tr}(b\partial_1c\,b\partial_2c)$, and $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$.
>
> The source rotation generators act on first jets by $J_+(\partial_2b)=\partial_1b$, $J_-(\partial_1b)=\partial_2b$, $J_+(\partial_1b)=J_-(\partial_2b)=0$, and the same formulas for $c$, with the Leibniz rule. Their scalar projector on these seeds is $1-J_-J_+/2$. Its odd-$b$ seed is zero. For each of the two degree-two seeds, its scalar projection minus, respectively, $\frac12\operatorname{tr}(b^2)C_0$ and $\operatorname{tr}(b^2)C_0$ is
>
> $$
> Q_0\left[
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}((\partial_2b)b\partial_1c)
> \right].
> $$
>
> These are full component identities, including bare ghosts. Multiplication by $Q_0$-closed powers of $\operatorname{tr}(b^2)$ proves spanning of the entire leading tree cohomology by $\mathbb Q[\operatorname{tr}(b^2)]C_0$. Its Cartan restriction is nonzero in every even degree and every tree boundary restricts to zero. Thus this module is free, with the displayed tree dimensions.
>
> Two additional exact identities are sufficient for the first differential:
>
> $$
> Q_1\bigl(\operatorname{tr}(b^2)C_0\bigr)
> =Q_0\bigl([\hbar]C_2\bigr),
> $$
>
> $$
> Q_1\bigl(\operatorname{tr}(b^2)^2\bigr)
> -2\operatorname{tr}(b^2)Q_1\operatorname{tr}(b^2)
> =16Q_0\left[
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}((\partial_1b)b\partial_2b)
> \right].
> $$
>
> Their primitives have $33$ and $6$ component monomials. Because $Q_1$ has differential order two and $Q_1C_0=0$, for each integer exponent $m\ge2$ the exact chain rule gives
>
> $$
> \begin{aligned}
> Q_1\bigl(\operatorname{tr}(b^2)^mC_0\bigr)
> =Q_0\Bigl[&m\operatorname{tr}(b^2)^{m-1}[\hbar]C_2\\
> &+8m(m-1)\operatorname{tr}(b^2)^{m-2}C_0
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}((\partial_1b)b\partial_2b)\Bigr].
> \end{aligned}
> $$
>
> The cases $m=0,1$ follow from $Q_1C_0=0$ and the first identity. Thus every leading tree class has zero outgoing $d_1$. An incoming $d_1$ would need $k+1$ copies of $b$ and one differentiated ghost; its minimum total weight $2k+3$ exceeds $2k+2$. There is no incoming map, proving the stated $E_2$ dimensions. $\square$

`all_degree_first_lift.pkl` saves the finite identities. The verifier also directly checks the full polynomial equality at powers $m=1,2,3,4,6$; these are checks of the universal chain-rule proof, not its range. Higher-loop survival of the entire tower is not asserted.
