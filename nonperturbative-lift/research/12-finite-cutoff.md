<a id="np-section-10"></a>

## 12. The finite-cutoff normal connection and section transport

<a id="np-0b135cc734d7230a-1"></a>

### Projection coefficients and boundary balance

Let $u_j=DH_j$ be the three orientation zero modes.

Their exact finite-ball Gram matrix is

$$
\int_{|y|\le\sqrt{s}}d^4y\,
(u_i)^j_a(u_k)_j{}^a
=\frac{4\pi^2\rho^2s^3}{(s+\rho^2)^3}\delta_{ik}.
$$

At $s=L^2$, the coefficient transferred from the normal connection to the collective orientation is

$$
P_L=\frac{g^2}{4}
\left(h_0+s h'-\frac{s h}{s+\rho^2}\right)
+4(s+\rho^2)k_{\eta^{\dot1}}(s).
$$

The corresponding direction is $P_L\eta^{\dot1}(u_1-iu_2)$.

Its exact large-cutoff limits are

$$
\lim_{L\to\infty}\frac{P_L}{\log(L^2/\rho^2)}
=\frac{g^2\rho^2}{16},
$$

$$
\lim_{L\to\infty}
\left[P_L-\frac{g^2\rho^2}{16}\log(L^2/\rho^2)\right]
=-\frac{7g^2\rho^2}{48}.
$$

This logarithmic transfer comes from the overlap at the spatial boundary.

The bilinear action pairing vanishes because its color contraction is

$$
(-i,-1,0)(1,-i,0)^T=0.
$$

The Hermitian pairing conjugates the first vector and gives $2i$, retaining the logarithmic overlap.

The full first-order gauge-action surface flux vanishes in every channel.

For a weighted composite insertion, the Ward identity retains the paired surface and bulk terms.

The finite-cutoff coordinate completion is the odd shear

$$
\eta^{\dot2}\longmapsto\eta^{\dot2}+4b_HP_Ln^3\eta^{\dot1}.
$$

Its Berezinian is one.

The added orientation and odd-coordinate divergences are, respectively,

$$
ib_HP_L(n^1-in^2)\eta^{\dot1},
\qquad
-ib_HP_L(n^1-in^2)\eta^{\dot1},
$$

Thus the full collective measure has zero divergence.

Write $Q_L$ for this finite-cutoff vector.

It satisfies $Q_L^2=0$ on every retained coordinate.

<a id="np-593328a6baa0c9e0-1"></a>

### Finite transport of the homogeneous section

Changing the homogeneous constant $h_0$ transports the full collective-coordinate section.

The source construction gives the exact coordinate flow $F_{\Delta h_0}$ at fixed $b_H,g,L$.

For example, the finite maps of size and body-frame color are

$$
\rho_{\Delta h_0}^{2}
=\rho^2[1-b_Hg^2\Delta h_0(n^1-in^2)],
$$

$$
n^i_{\Delta h_0}
=\frac{n^i-b_Hg^2\Delta h_0
[1-\tfrac12b_Hg^2\Delta h_0(n^1-in^2)](1,-i,0)^i}
{1-b_Hg^2\Delta h_0(n^1-in^2)}.
$$

The odd coordinates satisfy

$$
\xi_{\alpha,\Delta h_0}=\xi_\alpha,
\qquad
\eta^{\dot1}_{\Delta h_0}
=\frac{\eta^{\dot1}}{1-b_Hg^2\Delta h_0(n^1-in^2)},
$$

$$
\begin{aligned}
\eta^{\dot2}_{\Delta h_0}={}&
[1-b_Hg^2\Delta h_0(n^1-in^2)]\eta^{\dot2}\\
&-\frac{b_Hg^2\Delta h_0n^3
[2-b_Hg^2\Delta h_0(n^1-in^2)]}
{2[1-b_Hg^2\Delta h_0(n^1-in^2)]}\eta^{\dot1}.
\end{aligned}
$$

Integrating the full framing map along this flow gives

$$
F_{\Delta h_0'}\circ F_{\Delta h_0}=F_{\Delta h_0+\Delta h_0'},
\qquad F_{\Delta h_0}^{-1}=F_{-\Delta h_0}.
$$

> **Proof: density and Ward transport.**
>
> The size Jacobian is $1-b_Hg^2\Delta h_0(n^1-in^2)$, the Haar factor is its inverse square, and the odd determinant is one.
>
> Therefore
>
> $$
> F_{\Delta h_0}^*d\mu_0
> =\frac{d\mu_0}{1-b_Hg^2\Delta h_0(n^1-in^2)}.
> $$
>
> The coordinate change conjugates the Ward vector:
>
> $$
> Q_{L,\Delta h_0}
> =F_{\Delta h_0}^*Q_L(F_{\Delta h_0}^*)^{-1},
> \qquad Q_{L,\Delta h_0}^{2}=0.
> $$
>
> Relative to the original measure, the transformed vector has divergence
>
> $$
> \operatorname{div}_{\mu_0}Q_{L,\Delta h_0}
> =\frac{ib_Hg^2\Delta h_0(n^1-in^2)\eta^{\dot1}}
> {4[1-b_Hg^2\Delta h_0(n^1-in^2)]}.
> $$
>
> Applying $Q_{L,\Delta h_0}$ to the logarithm of the pulled-back density gives the negative of this expression.
>
> The full density therefore preserves Ward balance.
>
> $\square$

The flow is a complexified local or formal coordinate map.

The transformed integral uses the pulled-back insertion and measure together with the inverse-image cycle $F_{\Delta h_0}^{-1}\Gamma$.

Its size faces are

$$
\rho^2[1-b_Hg^2\Delta h_0(n^1-in^2)]=\varepsilon^2,
\qquad
\rho^2[1-b_Hg^2\Delta h_0(n^1-in^2)]=R^2.
$$

These exact identities transport the specified collective-coordinate section.

Matching the physical integration cycle and fluctuation/BV measures is a separate physical problem.
