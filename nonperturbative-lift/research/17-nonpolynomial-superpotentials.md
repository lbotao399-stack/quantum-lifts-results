<a id="np-section-15"></a>

## 17. Nonpolynomial superpotential: elimination, branches, and effective cohomology

<a id="np-acc5558e9ee1d6e4-1"></a>

### The ADS pole from eliminating massive fields

Eliminating a massive field means solving its F-equation and substituting the solution into $W$.

A denominator in that solution specifies the region where this description applies.

If the solution is multivalued, each root gives a branch.

The original polynomial equations retain all branches before elimination.

In Seiberg's constrained effective description of $SU(2)$ with $N_f=2$, $X$ is a Lagrange multiplier, and $m_2$ is the second-flavor mass:

$$
W=X\left(M^1{}_{\tilde1}M^2{}_{\tilde2}
-M^1{}_{\tilde2}M^2{}_{\tilde1}-B\widetilde B-\Lambda^4\right)
+m_2M^2{}_{\tilde2}.
$$

> **Proof: rational effective function.**
>
> The heavy equations give
>
> $$
> XM^1{}_{\tilde1}+m_2=0,
> \qquad
> XM^1{}_{\tilde2}=XM^2{}_{\tilde1}=XB=X\widetilde B=0.
> $$
>
> For $m_2\ne0$, both $X$ and the light meson $M=M^1{}_{\tilde1}$ are nonzero, so
>
> $$
> X=-\frac{m_2}{M},
> \qquad M^2{}_{\tilde2}=\frac{\Lambda^4}{M},
> \qquad
> W_{\mathrm{eff}}=\frac{m_2\Lambda^4}{M}
> =\frac{\Lambda_L^5}{M},
> \qquad \Lambda_L^5=m_2\Lambda^4.
> $$
>
> The light F-equation $-\Lambda_L^5/M^2=0$ has no solution at finite nonzero $M$.
>
> The original constrained equations also exclude $M=0$.
>
> $\square$

The pole records the domain of the description after field elimination.

Fractional ADS powers similarly record condensate sheets.

On the mesonic patch for general $N_f<N_c$, the integrated-in source glueball $S$ satisfies

$$
W(S,M)=S\left[
\log\frac{\Lambda^{3N_c-N_f}}{S^{N_c-N_f}\det M}
+N_c-N_f\right],
$$

$$
\frac{\partial W}{\partial S}=0
\quad\Longrightarrow\quad
S^{N_c-N_f}\det M=\Lambda^{3N_c-N_f},
\qquad W_{\mathrm{eff}}=(N_c-N_f)S.
$$

The polynomial equation retains all roots.

Choosing one sheet gives the fractional-power function.

[Seiberg 1994](https://arxiv.org/abs/hep-th/9402044); [Intriligator–Seiberg, §2.3, §4.1](https://arxiv.org/abs/hep-th/9509066)

Clearing denominators in rational F-equations requires retaining their domain.

For example, adding a mass to the $SU(2)$ ADS model gives

$$
\frac{m_1M^2-\Lambda_L^5}{M^2}=0,
\qquad M\ne0.
$$

Multiplying the numerator equation by $M$ and retaining the polynomial equation $M(m_1M^2-\Lambda_L^5)=0$ introduces a spurious $M=0$ component.

Saturation preserves the original open patch within polynomial algebra.

Here the exact identity is

$$
(M(m_1M^2-\Lambda_L^5)):M^\infty
=(m_1M^2-\Lambda_L^5)
\qquad(\Lambda_L^5\ne0).
$$

The ideal $I:M^\infty$ consists of polynomials that enter $I$ after multiplication by some nonnegative power of $M$.

Equivalently, form the ideal in the algebra localized at $M\ne0$, then contract it to the polynomial algebra.

<a id="np-d680ecc8328087d5-1"></a>

### Unit exactness and the coefficient ring

The Koszul complex used here consists of scalar functions and their odd partners, with the supercharge sending each odd partner to its F-equation.

The unit is the constant operator $1$.

If $1$ is a $Q$-variation, multiplication by its odd primitive makes every closed operator exact.

The following formulas state separately the allowed formal coefficients and the model at fixed nonzero coupling.

In the effective ADS Koszul complex, $\gamma^{(0)}=\bar M$ and $\beta^{(0)}$ are an even letter and its odd partner:

$$
Q\beta^{(0)}=\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2},
\qquad
Q((\gamma^{(0)})^2\beta^{(0)})=\bar\Lambda_L^5.
$$

> **Proof: formal torsion and contraction at fixed coupling.**
>
> On the Laurent-polynomial algebra for $\gamma^{(0)}\ne0$,
>
> $$
> \mathbb C[[\bar\Lambda_L^5]][(\gamma^{(0)})^{\pm1}]\beta^{(0)}
> \xrightarrow{Q}
> \mathbb C[[\bar\Lambda_L^5]][(\gamma^{(0)})^{\pm1}].
> $$
>
> The differential multiplies by $\bar\Lambda_L^5/(\gamma^{(0)})^2$.
>
> Its kernel is zero, and its image is the ideal $(\bar\Lambda_L^5)$.
>
> Hence
>
> $$
> H^{-1}=0,
> \qquad H^0=\mathbb C[(\gamma^{(0)})^{\pm1}],
> \qquad \bar\Lambda_L^5 H^0=0.
> $$
>
> Over the coefficient field with the instanton parameter inverted, or in a separately specified model at fixed $\bar\Lambda_L^5\ne0$,
>
> $$
> Q\left(\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}\right)=1,
> $$
>
> $$
> Q\left(\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}\mathcal O\right)
> +\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}Q\mathcal O
> =\mathcal O.
> $$
>
> Every closed cochain is therefore exact.
>
> $\square$

The constrained model admits a certificate polynomial in the mesons.

Its global polynomial certificate is

$$
\begin{aligned}
&\Lambda^4\frac{\partial W}{\partial M^2{}_{\tilde2}}
-(M^1{}_{\tilde1})^2\frac{\partial W}{\partial M^1{}_{\tilde1}}
-M^1{}_{\tilde1}M^1{}_{\tilde2}\frac{\partial W}{\partial M^1{}_{\tilde2}}
-M^1{}_{\tilde1}B\frac{\partial W}{\partial B}
+M^1{}_{\tilde1}X\frac{\partial W}{\partial X}
=m_2\Lambda^4.
\end{aligned}
$$

At fixed nonzero $m_2\Lambda^4$, the unit belongs to the full effective Jacobian ideal.

Comparison with the full microscopic electric jet complex requires a separate physical input.

> **Example 17.1. **
>
> **The mass-deformed ADS algebra.**
>
> Adding $m_1M$ gives
>
> $$
> W_{\mathrm{eff}}=\frac{\Lambda_L^5}{M}+m_1M,
> \qquad M^2=\frac{\Lambda_L^5}{m_1},
> \qquad
> W_{\mathrm{vac},\pm}=\pm2\sqrt{m_1m_2\Lambda^4}.
> $$
>
> At fixed nonzero parameters, the effective cohomology becomes a two-dimensional algebra:
>
> $$
> H^{-1}=0,
> \qquad
> H^0=\frac{\mathbb C[(\gamma^{(0)})^{\pm1}]}
> {((\gamma^{(0)})^2-\bar\Lambda_L^5/\bar m_1)}.
> $$

<a id="np-9330362f2f4aa9a8-1"></a>

### Quantum constraint and s-confining relation

The antisymmetric meson $V^{ij}$ of $SU(2)$ with $N_f=2$ has six independent entries:

$$
\operatorname{Pf}V=V^{12}V^{34}-V^{13}V^{24}+V^{14}V^{23}=\Lambda^4.
$$

Its six partial derivatives are $(V^{34},-V^{24},V^{23},V^{14},-V^{13},V^{12})$, with common zero $V=0$.

For $\Lambda\ne0$, this point lies outside the hypersurface.

The quantum moduli space is therefore smooth, of complex dimension $5$.

For $N_f=3$, $V$ is a six-by-six antisymmetric matrix:

$$
W=-\frac{\operatorname{Pf}V}{\Lambda^3},
\qquad
\frac{\partial W}{\partial V^{12}}
=-\frac{V^{34}V^{56}-V^{35}V^{46}+V^{36}V^{45}}{\Lambda^3}.
$$

The independent derivatives impose the four-by-four Pfaffian minors, hence $\operatorname{rank}V\le2$.

At $V=0$, the first and second derivatives vanish, leaving fifteen massless meson chirals.

The nonperturbative interaction implements the electric composite relations while retaining the moduli space.
