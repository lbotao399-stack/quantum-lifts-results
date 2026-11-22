# A global effective contraction before eliminating the massive meson

Primary input: Seiberg, [hep-th/9402044](https://arxiv.org/pdf/hep-th/9402044), (4.1), (4.8), and the mass-decoupling paragraph following (4.8). Work in its constrained effective model for SU(2) with two flavors, adding the second flavor mass $m_2$:

$$
W=X\left(M^1{}_{\tilde1}M^2{}_{\tilde2}-M^1{}_{\tilde2}M^2{}_{\tilde1}-B\widetilde B-\Lambda^4\right)+m_2M^2{}_{\tilde2}.
$$

$X$ is the source Lagrange multiplier. The following polynomial identity holds on the entire space of these effective coordinates, without inverting any meson:

$$
\begin{aligned}
&\Lambda^4\frac{\partial W}{\partial M^2{}_{\tilde2}}
-(M^1{}_{\tilde1})^2\frac{\partial W}{\partial M^1{}_{\tilde1}}
-M^1{}_{\tilde1}M^1{}_{\tilde2}\frac{\partial W}{\partial M^1{}_{\tilde2}}
-M^1{}_{\tilde1}B\frac{\partial W}{\partial B}
+M^1{}_{\tilde1}X\frac{\partial W}{\partial X}\\
&=\Lambda^4(XM^1{}_{\tilde1}+m_2)
-X(M^1{}_{\tilde1})^2M^2{}_{\tilde2}
+XM^1{}_{\tilde1}M^1{}_{\tilde2}M^2{}_{\tilde1}
+XM^1{}_{\tilde1}B\widetilde B\\
&\quad+XM^1{}_{\tilde1}(M^1{}_{\tilde1}M^2{}_{\tilde2}-M^1{}_{\tilde2}M^2{}_{\tilde1}-B\widetilde B-\Lambda^4)
=m_2\Lambda^4.
\end{aligned}
$$

Each partial derivative is the differential of its odd partner in the effective Koszul complex, with the complex conjugation and overall minus sign fixed by the manuscript's $Q\beta^{(0)}=-\partial\bar W/\partial\gamma^{(0)}$. Complex-conjugating the displayed identity and reversing the coefficients therefore supplies an explicit odd polynomial whose $Q$ image is $\bar m_2\bar\Lambda^4$.

For fixed $m_2\Lambda^4\ne0$, the unit is exact in this entire constrained effective complex; division is by a nonzero numerical parameter, not a meson coordinate. Multiplication by the resulting odd primitive contracts every closed element. At fixed $m_2\ne0$ over formal nonnegative powers of $\Lambda^4$, the same identity gives instanton-parameter torsion instead. Specialization to $\Lambda^4=0$ has a critical point with $M^1{}_{\tilde1}=1$, $X=-m_2$, and the remaining fields zero, so the unit cannot have been exact over that formal ring.

This sharpens the previously localized ADS calculation in `n1_duality.md`: localization at the surviving meson is unnecessary if the constrained model is retained. The model's identification with the complete microscopic electric jet complex is a separate physical input; this polynomial computation constructs no ultraviolet primitive.

Verification: `check_global_ads_identity.py`, `global_ads_identity_readback.json`. All cancellations are exact polynomial identities.
