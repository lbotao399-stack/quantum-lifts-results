<a id="np-section-05"></a>

## 7. Characteristic residues and twisted periods

Let $z$ be the auxiliary Laurent variable.

The operation $[z^{-1}]$ extracts the $z^{-1}$ coefficient after expansion at infinity.

[Feigin–Silantyev, Theorem 4.8, Corollary 4.9, and §5](https://arxiv.org/abs/1110.1946) classify type-$A$ invariant polynomial twisted periods.

The following characteristic-residue family belongs to that classification.

A complete family of leading SU(3) Cartan representatives is

$$
[z^{-1}]\det(z\mathbf1-b)^{(k-2)/3}(a_1-zC_0).
$$

For SU(4), the representatives are

$$
\begin{aligned}
[z^{-1}]\det(z\mathbf1-b)^{(k-3)/4}\Bigg[&
\left(z^2-\frac14\operatorname{tr}(b^2)\right)C_0-z a_1\\
&+\left(1-\frac12J_-J_+\right)\operatorname{tr}(b^2\partial_1c\partial_2c)
-\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)\Bigg].
\end{aligned}
$$

These formulas specify Cartan representatives.

The preceding section states the degree ranges in which they lift to the full complex.

A uniform expression starts from the scalar period of degree $k+2$:

$$
[z^{-1}]\det(z\mathbf1-b)^{(k+1)/N}.
$$

Take its ordinary Hessian and contract its two slots, in order, with $\partial_1c,\partial_2c$, retaining this order of odd factors.

The source trace pairing fixes the Hessian contraction, with the additional $1/2$ absent.

> **Proof of the Hessian bridge.**
>
> Direct differentiation of the determinant on the Cartan subalgebra gives
>
> $$
> \begin{aligned}
> [z^{-1}]\det(z\mathbf1-b)^{(k+1)/N}\Bigg[&
> \frac{(k+1)^2}{N^2}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\bigr)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_2c\bigr)\\
> &-\frac{k+1}{N}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-2}\partial_1c\partial_2c\bigr)\Bigg].
> \end{aligned}
> $$
>
> A Laurent total derivative has zero $z^{-1}$ coefficient.
>
> Apply this fact to
>
> $$
> \frac{d}{dz}\left[
> \det(z\mathbf1-b)^{(k+1)/N}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\partial_2c\bigr)\right]
> $$
>
> to rewrite the preceding expression as
>
> $$
> \begin{aligned}
> \frac{(k+1)^2}{N^2}[z^{-1}]\det(z\mathbf1-b)^{(k-N+1)/N}\Bigg[&
> \det(z\mathbf1-b)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\bigr)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_2c\bigr)\\
> &-\frac{d}{dz}\det(z\mathbf1-b)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\partial_2c\bigr)\Bigg].
> \end{aligned}
> $$
>
> The poles of the bracket cancel at every simple eigenvalue, leaving a polynomial.
>
> At $N=3$ it is $3(a_1-zC_0)$; at $N=4$ it is $-4$ times the preceding SU(4) bracket.
>
> Thus the ordinary Hessian equals the saved SU(3) and SU(4) generating representatives multiplied by
>
> $$
> \frac{(k+1)^2}{3},\qquad -\frac{(k+1)^2}{4}
> $$
>
> respectively.
>
> A homogeneous polynomial of degree $k+2\ge2$ with zero Hessian is zero.
>
> Together with the complete Cartan-kernel classifications, this proves the correspondence for these two gauge groups.
>
> $\square$

**Proposition P6.**

For every $N\ge2$, the scalar period residue and its Hessian vanish exactly when $k\equiv N-1\pmod N$.

> **Proof.**
>
> For $k=N-1+Nm$, $m\ge0$, the exponent $(k+1)/N=m+1$ is an integer.
>
> The integrand is a polynomial in $z$, with zero $z^{-1}$ coefficient.
>
> In the remaining degrees, the exponent is nonintegral.
>
> For $N\ge3$, the characteristic-coefficient weights include $2,3$, and every $k+2\ge2$ is a nonnegative integer combination of these weights.
>
> The corresponding generalized multinomial coefficient is a product of nonintegral falling factors divided by positive integer factorials, hence is nonzero.
>
> Distinct algebraically independent characteristic monomials cannot cancel.
>
> For $N=2$, the remaining $k$ are even, and the $(k+2)/2$ power of the quadratic coefficient gives a nonzero term.
>
> Hessian injectivity gives the same vanishing pattern.
>
> $\square$

This proves the period-side statement for arbitrary $N$.

Its correspondence with the full Cartan $Q_1$ kernel is established for $SU(3),SU(4)$.

The arbitrary-$N$ identification with the quantum complex remains a research problem.
