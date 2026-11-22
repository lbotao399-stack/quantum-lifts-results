# SU(3): the common-root and algebraic partners agree

Use the original SU(3) component-engine normalization and the source polynomials $C_0,x_3,a_1,w_1,w_0,n_1,\theta,\eta$ recorded in `../su3_rational/proofs.md`. The superscript $(1)$ denotes the coefficient of the stored loop parameter. The common-root homotopy has the raw-kernel-to-engine factor $2$, fixed by the associator equation.

**Proposition 1.** The first corrected coefficient of the iterated multiplication homotopy is

$$
\eta_{\rm diag}(a_1)=\mu_1(x_3,a_1)
-H^{(1)}_{C_0^2,C_0}(a_1)
-H^{(1)}_{C_0,C_0}(C_0a_1).
$$

It gives the same class in $H(Q_0)^{(0,7)}(5,5)$ as the source's algebraically solved $\eta(a_1)$.

> **Proof.** The telescoping multiplication identity in `iterated_multiplication.md` gives
>
> $$
> h_{\rm diag}=x_3\star{}-H_{C_0^2,C_0}-H_{C_0,C_0}L,
> \qquad \{Q_{\rm pert},h_{\rm diag}\}=L^3.
> $$
>
> Expanding its coefficient of the loop parameter gives the displayed formula. Exact rational evaluation yields
>
> $$
> \begin{aligned}
> \#\mu_1(x_3,a_1)&=94536,\\
> \#H^{(1)}_{C_0^2,C_0}(a_1)&=30284,\\
> \#H^{(1)}_{C_0,C_0}(C_0a_1)&=17430,\\
> \#\eta(a_1)&=95100,\\
> \#\bigl(\eta(a_1)-\eta_{\rm diag}(a_1)\bigr)&=95786.
> \end{aligned}
> $$
>
> Here $\#$ counts canonical component monomials. The full tree differential, including undifferentiated ghosts, and the Cartan restriction satisfy
>
> $$
> Q_0\bigl(\eta(a_1)-\eta_{\rm diag}(a_1)\bigr)=0,
> \qquad
> \left.\bigl(\eta(a_1)-\eta_{\rm diag}(a_1)\bigr)\right|_{\mathfrak t}=0.
> $$
>
> The complete rational target cohomology has dimension one, and its nonzero class has nonzero Cartan restriction (`../su3_rational/proofs.md`, Propositions 3–4). Hence the displayed difference is $Q_0$-exact. $\square$

**Corollary 2.** Replacing the source's algebraic $\eta$ by the common-root coefficient preserves the full SU(3) partner class and its coefficient $40$ in the source Cartan basis.

> **Proof.** The reduced source partner is
>
> $$
> \eta(a_1)-C_0^2w_0-2C_0\theta(w_1)-Q_1n_1.
> $$
>
> Proposition 1 changes it by a $Q_0$ boundary in the pure-ghost sector. Every higher $Q_r$ annihilates a pure-ghost primitive. Thus this is also a $Q_{\rm pert}$ boundary. The same-degree source of multiplication by $C_0^2$ has
>
> $$
> E_\infty^{N_c=3}(3,3)=0
> $$
>
> for every permitted $N_b$, by `../su3_rational/C0_squared_image_certificate.json`. The Massey quotient therefore preserves the nonzero partner class. $\square$

Evidence: `su3_common_eta_cartan.py` checks the survivor-projected operations against full operations before projecting. `su3_common_eta_full.py` evaluates the full difference and its full $Q_0$ differential. Saved polynomials are `su3_eta_choice_exact.pkl`; the readback is `su3_eta_choice_exact.json`. The primitive exists by the complete cohomology detector argument; its coefficients have not been constructed.
