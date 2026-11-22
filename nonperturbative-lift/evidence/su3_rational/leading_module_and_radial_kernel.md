# SU(3) leading tree module and all-degree Cartan kernel

The source fields, traces and loop parameter are unchanged. The leading sector has $N_b=k$, $N_c=2$, and weight $(k+1,k+1)$. Both ghosts are first jets and every $b$ is undifferentiated. Throughout this note, cohomology and scalar mean the source gauge-invariant rotation-scalar complex. The source expressions are

$$
C_0=\operatorname{tr}(\partial_1c\partial_2c),
\qquad
a_1=-\frac12\operatorname{tr}(b\partial_1c\partial_2c-b\partial_2c\partial_1c).
$$

## Finite proof of the complete tree module

> **Proof.**
>
> Invariant tensor contractions are products of traces. With exactly two ghosts, their traces either share a trace or occupy two traces. The other traces contain only $b$. The source Cayley--Hamilton identity
>
> $$
> b^3=\frac12\operatorname{tr}(b^2)b+\frac13\operatorname{tr}(b^3)\mathbf1
> $$
>
> reduces each bare-$b$ gap to length $0,1,2$ and all pure-$b$ traces to polynomials in $\operatorname{tr}(b^2),\operatorname{tr}(b^3)$. Tracelessness removes a trace containing a single ghost and no $b$. Graded cyclicity and scalar projection identify the two gap orders. There remain six single-trace and three double-trace generators, all of degree at most four in $b$.
>
> The source rotation generators $J_+,J_-$ act on first jets by $J_+(\partial_2b)=\partial_1b$, $J_-(\partial_1b)=\partial_2b$, and the same formulas for $c$, with the Leibniz rule. Apply $1-J_-J_+/2$ to every seed and every primitive in the following table. The scalar-projected seed minus the displayed representative is exactly $Q_0$ of the scalar-projected primitive.
>
> | Seed | Representative modulo $Q_0$ | Primitive seed |
> |:---|:---|:---|
> | $\operatorname{tr}(\partial_1c\partial_2c)$ | $C_0$ | $0$ |
> | $\operatorname{tr}(\partial_1c\,b\partial_2c)$ | $-a_1$ | $0$ |
> | $\operatorname{tr}(\partial_1c\,b^2\partial_2c)$ | $\frac16\operatorname{tr}(b^2)C_0+\frac13\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $-\frac13\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1c\,b\partial_2c)$ | $\frac16\operatorname{tr}(b^2)C_0+\frac13\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\frac23\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1c\,b^2\partial_2c)$ | $\frac13\operatorname{tr}(b^3)C_0-\frac12\operatorname{tr}(b^2)a_1$ | $\operatorname{tr}((\partial_2b)b^2\partial_1c)$ |
> | $\operatorname{tr}(b^2\partial_1c\,b^2\partial_2c)$ | $\frac1{12}\operatorname{tr}(b^2)^2C_0-\frac13\operatorname{tr}(b^3)a_1+\frac16\operatorname{tr}(b^2)\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\operatorname{tr}(b(\partial_2b)b^2\partial_1c)+\frac13\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b\partial_1c)$ |
> | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $0$ |
> | $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b^2\partial_2c)$ | $\frac12\operatorname{tr}(b^3)C_0-\frac12\operatorname{tr}(b^2)a_1$ | $\operatorname{tr}((\partial_2b)b^2\partial_1c)$ |
> | $\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c)$ | $\frac16\operatorname{tr}(b^2)^2C_0-\operatorname{tr}(b^3)a_1-\frac16\operatorname{tr}(b^2)\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$ | $\operatorname{tr}(b(\partial_2b)b^2\partial_1c)+\frac23\operatorname{tr}(b^2)\operatorname{tr}((\partial_2b)b\partial_1c)$ |
>
> Every row was expanded and verified using full $Q_0$, including bare ghosts. The nonzero primitives have $88,88,256,816,256,812$ component monomials. Since $Q_0\operatorname{tr}(b^2)=Q_0\operatorname{tr}(b^3)=0$, these nine identities prove spanning in every degree.
>
> For independence, $b_1,b_2$ are the coefficients of $b$ in the original trace-dual Cartan basis; $c^1,c^2$ are the coefficients of $c$ in the dual Chevalley Cartan basis. The three generators' symmetric ghost coefficient matrices, read in the slots $(1,1),(1,2),(2,2)$, have determinant
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
> Cartan restriction kills every $Q_0$ boundary. Thus the three generators are free over $\mathbb Q[\operatorname{tr}(b^2),\operatorname{tr}(b^3)]$. No localization is used: the reduction table proves spanning before restriction, and the nonzero polynomial determinant proves independence.
>
> With $t$ counting the number of $b$ fields, the Hilbert series is
>
> $$
> \frac{1+t+t^2}{(1-t^2)(1-t^3)}
> =\frac1{(1-t)(1-t^2)},
> \qquad
> \dim H(Q_0)^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
> =\left\lfloor\frac k2\right\rfloor+1.
> $$
>
> Independently, the Weyl group is $S_3$ and the symmetric square of its two-dimensional Cartan representation has characters $3,1,0$ on the identity, transposition and three-cycle classes. Their class cardinalities are $1,3,2$, giving the same Molien series:
>
> $$
> \frac16\left(\frac3{(1-t)^2}+\frac3{1-t^2}+\frac0{1+t+t^2}\right)
> =\frac1{(1-t)(1-t^2)}.
> $$
>
> $\square$

## Exact radial one-loop calculation

> **Proof.**
>
> Each leading tree class has three polynomial coefficients, multiplying $C_0$, $a_1$, and $\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)$, respectively. These are polynomials in the two independent generators $\operatorname{tr}(b^2),\operatorname{tr}(b^3)$; partial derivatives below are with respect to those actual polynomial coordinates.
>
> The source $Q_1$ selects two input letters, so its differential order is at most two. Its action on arbitrary such coordinates is fixed by its values on each module generator multiplied by
>
> $$
> 1,\quad\operatorname{tr}(b^2),\quad\operatorname{tr}(b^3),\quad
> \operatorname{tr}(b^2)^2,\quad
> \operatorname{tr}(b^2)\operatorname{tr}(b^3),\quad
> \operatorname{tr}(b^3)^2.
> $$
>
> All eighteen images were evaluated over $\mathbb Q$. Their second-Casimir-derivative coefficients vanish after Cartan restriction. The formula reproduces all twenty-four previously saved leading-basis images for $k=1,\ldots,6$ and two independent higher-degree polynomial images.
>
> Only the zeroth operator multiplying the coefficient of the third trace generator contains a second-jet ghost. In particular, the coefficient of
>
> $$
> (\partial_2c^1)(\partial_1\partial_2c^1)(\partial_1c^2)
> $$
>
> is $6(b_1+2b_2)$ times the coefficient of the third trace generator. The nonzero polynomial factor forces that coefficient to vanish.
>
> The remaining two independent operators are the Cartan restrictions of $Q_1(\operatorname{tr}(b^2)C_0)$ and $Q_1(\operatorname{tr}(b^3)C_0)$. Their coefficients of $(\partial_2b_1)(\partial_2c^1)(\partial_1c^1)(\partial_1c^2)$ and $(\partial_2b_2)(\partial_2c^1)(\partial_1c^1)(\partial_1c^2)$ have determinant
>
> $$
> \det\begin{pmatrix}
> -12&6(-b_1+b_2)\\
> -24&6(b_1+2b_2)
> \end{pmatrix}
> =-216b_1\ne0.
> $$
>
> Direct polynomial identities identify the two $a_1$ derivative operators with $2/3$ times the second of these operators and $\operatorname{tr}(b^2)/4$ times the first. Consequently the Cartan kernel is exactly
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
> Equality of the mixed derivatives of the $C_0$ coefficient gives
>
> $$
> \frac{\partial^2(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)^2}
> =\frac{3\operatorname{tr}(b^2)}8
> \frac{\partial^2(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)^2}.
> $$
>
> For $k>0$, the $a_1$ coefficient has weighted degree $k-1$ when the two trace coordinates have degrees $2,3$. The compatible $C_0$ coefficient, of weighted degree $k$, is uniquely
>
> $$
> \text{coefficient of }C_0=-\frac1k\left[
> \frac{\operatorname{tr}(b^2)^2}{2}
> \frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^3)}
> +2\operatorname{tr}(b^3)
> \frac{\partial(\text{coefficient of }a_1)}{\partial\operatorname{tr}(b^2)}
> \right].
> $$
>
> Weighted Euler differentiation and the preceding second-order equation verify both first-order equations directly. For $k=0$, only a constant multiple of $C_0$ remains.
>
> $\square$

## Recurrence and exact all-degree kernel dimension

> **Proof.**
>
> The nonnegative exponents $r,s$ in the Casimir monomials of the $a_1$ coefficient obey $2r+3s=k-1$. Comparing their coefficients gives
>
> $$
> r(r-1)\left(\text{coefficient of }\operatorname{tr}(b^2)^r\operatorname{tr}(b^3)^s\text{ in the }a_1\text{ coordinate}\right)
> =\frac38(s+2)(s+1)\left(\text{coefficient of }\operatorname{tr}(b^2)^{r-3}\operatorname{tr}(b^3)^{s+2}\text{ in the }a_1\text{ coordinate}\right),
> $$
>
> where a coefficient with a negative exponent is zero. At $k=2$, there is no Casimir monomial of weighted degree one, so the $a_1$ coefficient is zero directly. For $k=5,8,\ldots$, the smallest possible exponent $r$ is $2$; its coefficient vanishes and the recurrence kills all coefficients. Together these are exactly the cases $k\equiv2\pmod3$. If the smallest exponent is $0$ or $1$, its second derivative vanishes, its coefficient is free, and every subsequent coefficient is fixed by the recurrence. The recurrence terminates when $s$ is $0$ or $1$; there is no additional endpoint equation. Thus
>
> $$
> \dim\ker(\mathrm{Cartan}\circ Q_1)\big|_{H(Q_0)^{(k,2)}_{\mathrm{scalar}}}
> =\begin{cases}
> 0,&k\equiv2\pmod3,\\
> 1,&k\not\equiv2\pmod3.
> \end{cases}
> $$
>
> The three full pure-Casimir brackets are $Q_0$-exact. Their primitives, with the source scalar projection $1-J_-J_+/2$, are respectively
>
> $$
> 16\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}((\partial_1b)b\partial_2b),
> $$
>
> $$
> 24\left(1-\frac{J_-J_+}{2}\right)\operatorname{tr}((\partial_1b)b^2\partial_2b),
> $$
>
> $$
> 36\left(1-\frac{J_-J_+}{2}\right)
> \left[\operatorname{tr}(b(\partial_1b)b^2\partial_2b)
> +\operatorname{tr}((\partial_1b)b^3\partial_2b)\right].
> $$
>
> The brackets here mean, in order, $Q_1(\operatorname{tr}(b^2)^2)-2\operatorname{tr}(b^2)Q_1\operatorname{tr}(b^2)$, $Q_1(\operatorname{tr}(b^2)\operatorname{tr}(b^3))-\operatorname{tr}(b^2)Q_1\operatorname{tr}(b^3)-\operatorname{tr}(b^3)Q_1\operatorname{tr}(b^2)$, and $Q_1(\operatorname{tr}(b^3)^2)-2\operatorname{tr}(b^3)Q_1\operatorname{tr}(b^3)$.
>
> The two full first-derivative relations
>
> $$
> Q_1(\operatorname{tr}(b^2)a_1)-\frac23Q_1(\operatorname{tr}(b^3)C_0),
> $$
>
> $$
> Q_1(\operatorname{tr}(b^3)a_1)
> -\frac{\operatorname{tr}(b^2)}4Q_1(\operatorname{tr}(b^2)C_0)
> $$
>
> also have exact full-$Q_0$ primitives, containing $2546$ and $10744$ terms. These were extracted from the certified $C_3,C_4$ first corrections with their exact leading-representative changes included.
>
> The order-two chain rule now gives a universal primitive. Multiply the first-derivative primitives by the respective derivatives of the $a_1$ coefficient. Multiply each pure-Casimir bracket primitive by the corresponding second derivative of the $C_0$ coefficient times $C_0$ plus that of the $a_1$ coefficient times $a_1$; include $1/2$ for the equal-coordinate second derivatives and $1$ for the mixed derivative. All multiplier polynomials are $Q_0$-closed. The remaining full-$Q_1$ terms are precisely the two first-order kernel equations above. Hence every Cartan-kernel solution has a full $Q_0$ primitive of its $Q_1$ image.
>
> An incoming $d_1$ would require $N_b=k+1$, $N_c=1$, with minimum total weight $2k+3>2k+2$. Therefore, for every $k\ge0$,
>
> $$
> \dim_{\mathbb Q}E_2^{(k,2)}(k+1,k+1)_{\mathrm{scalar}}
> =\begin{cases}
> 0,&k\equiv2\pmod3,\\
> 1,&k\not\equiv2\pmod3.
> \end{cases}
> $$
>
> Independent readback verified all five full-$Q_0$ primitive identities and four full first-derivative $Q_1$ images. At $k=7$, the direct full $Q_1$ image has $114724$ terms and the universal primitive has $121438$ terms, with zero residual. An off-kernel test also reproduced its $27100$-term nonexact remainder exactly. These are verification cases for the universal chain rule, not the basis of an extrapolation in $k$. Higher-loop survival remains a separate question.
>
> $\square$

## Saved exact evidence

- `leading_module_Q.pkl` and its certificate contain all nine full-$Q_0$ primitive identities and the source Cartan determinant.
- `leading_radial_Q1.pkl` and its certificate contain the eighteen images, exact derivative operators and finite $k=0,\ldots,30$ census. The all-degree theorem comes from the recurrence above.
- `leading_radial_kernel_Q_certificate.json` records the derivative-sector elimination and nonzero minor.
- `radial_first_primitives_Q.pkl` contains full-$Q_0$ primitives for the two first-derivative relations, with source leading-representative corrections included.
- `../su3_transgression/casimir_brackets.pkl` and `../su3_transgression/verify_universal_first_lift.pkl` contain the three pure-Casimir bracket primitives and independent universal first-lift readback.
