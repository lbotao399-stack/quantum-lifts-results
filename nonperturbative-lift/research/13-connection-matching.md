<a id="np-section-11"></a>

## 13. The complete connection contribution and the open physical coefficient

The physical source polynomial of the original $C_2$ partner is

$$
\hbar M_1+\hbar^2M_0,
$$

Here $M_1$ has one $b$ and five ghosts, while $M_0$ has zero $b$ and five ghosts.

Their source expansions contain $597$ and $973$ monomials, respectively.

This section evaluates the coefficient of $b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2$ in the normal connection variation.

The connection appears inside each covariant derivative of a field.

Varying it produces commutators with the field and, for higher jets, derivatives of the connection variation.

The symbol $\delta_R$ applies this replacement to every covariant jet in $M_1$, with the source symmetrization weights.

The subscript $b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2$ extracts the coefficient of that ordered monomial: second order in the even curvature parameter and top degree in the six indicated fermion parameters.

The normal connection is

$$
R_a=(Q_{\rm physical}-Q_{\rm old})\delta(gV)_a
-P_L\eta^{\dot1}(u_1-iu_2)_a,
\qquad
R_\alpha=\sigma^a_{\alpha\dot1}R_a.
$$

The first term measures the difference between the physical supersymmetry variation of the connection response and its original collective-coordinate variation.

The subtraction removes its projected orientation-zero-mode component.

Thus $R$ is the remaining normal variation of the field, while $P_L$ records the motion transferred to the collective orientation at finite spatial radius $L$.

<a id="np-361905ea0e0c5763-1"></a>

### Source symmetrization and the angular integral

The physical dictionary symmetrizes the derivative slots and gaugino slot together with total weight one.

For the ghost jet $c_{rt}$, a connection jet of degree $(a,b)$ has source weight

$$
\frac{r+t-a-b}{r+t}\binom ra\binom tb,
\qquad 1\le a+b<r+t.
$$

Here $r,t,a,b$ are derivative counts.

In particular,

$$
\delta_Rc_{20}=[R_1,c_{10}],
\qquad
\delta_Rc_{30}=2[R_1,c_{20}]+[D_1R_1,c_{10}],
\qquad
\delta_R(D_\alpha b)=[R_\alpha,b].
$$

Acting on all $597$ terms gives $480$ nonzero terms from the $b$ variation and $1262$ from the ghost variation.

The two numerators contain $11190$ and $21225$ monomials, respectively, with common denominator $(s+\rho^2)^{16}$.

The exact angular average is

$$
\begin{aligned}
&\frac1{2\pi^2}\int_{S^3}d\Omega\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}\\
={}&\frac{4\rho^8s^3(51s-463\rho^2)}{9(s+\rho^2)^{12}}
(4P_L-g^2h_0)\\
&+\frac{32\rho^6s^4(80s+1077\rho^2)}{9(s+\rho^2)^{11}}k_1
+\frac{32\rho^6s^4(161s+130\rho^2)}{9(s+\rho^2)^{10}}k_1'\\
&-\frac{16\rho^6s^4(179s^2-798\rho^2s-2971\rho^4)}{27(s+\rho^2)^{11}}k_3
-\frac{16\rho^6s^5(259s+197\rho^2)}{27(s+\rho^2)^{10}}k_3'\\
&-\frac{32\rho^6s^3(140s^2+211\rho^2s-926\rho^4)}{9(s+\rho^2)^{11}}k_{\eta^{\dot1}}
-\frac{32\rho^6s^4(140s+109\rho^2)}{9(s+\rho^2)^{10}}k_{\eta^{\dot1}}'.
\end{aligned}
$$

A prime denotes $d/ds$.

An independent calculation with full exterior multiplication agrees with the cofactor reduction.

The two normal-connection columns built from external odd parameters vanish in this top component: the associated one-forms, two-forms, and remaining ghost forms lie in one five-dimensional annihilator space, so their degree-six wedge vanishes.

<a id="np-c209f4b75ff549bb-1"></a>

### The radial adjoint and finite boundary terms

The position measure is $\pi^2s\,ds$.

In each channel, let $A,B$ be the coefficients of $k,k'$ after including this measure.

The rational adjoint $z$ satisfies

$$
4sz''-4lz'
+\left[\frac{4c}{s+\rho^2}-\frac{8s}{(s+\rho^2)^2}\right]z=A-B'.
$$

> **Proof: the full finite radial identity.**
>
> Integrate $Bk'$ once by parts and $z(4sk''+4(l+2)k')$ twice to obtain
>
> $$
> \int_0^{L^2}(Ak+Bk')ds
> =\int_0^{L^2}zf\,ds
> +[Bk-4s(zk'-z'k)-4(l+1)zk]_0^{L^2}.
> $$
>
> Regularity makes the origin faces vanish; the upper faces remain at finite $L$.
>
> The exact beta and logarithmic moments of the three rational adjoints with their forcing terms give the integral below.
>
> $\square$

For the Green representative with upper endpoint at infinity and fixed finite $P_L$,

$$
\int_{\mathbb R^4}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}
=-\pi^2\left[
\frac{3364P_L}{10395\rho^4}
+\frac{724901g^2}{27941760\rho^2}\right].
$$

A change of the homogeneous constant $h_0$ cancels pointwise:

$$
\partial_{h_0}k_1=\partial_{h_0}k_{\eta^{\dot1}}
=-\frac{g^2\rho^2}{16(s+\rho^2)^2},
\qquad
\partial_{h_0}k_3=\partial_{h_0}P_L=0.
$$

An admissible scalar homogeneous change acts simultaneously as

$$
k_{\eta^{\dot1}}\mapsto k_{\eta^{\dot1}}+\frac{C}{s+\rho^2},
\qquad P_L\mapsto P_L+4C.
$$

Its additional upper boundary face is $13456\pi^2C/(10395\rho^4)$, which cancels the change in the displayed bulk functional.

Using the actual finite-ball $P_L$ gives the exact cutoff coefficients of this connection sector:

$$
\lim_{L\to\infty}
\frac{\rho^2}{\pi^2g^2\log(L^2/\rho^2)}
\int_{|y|\le L}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}
=-\frac{841}{41580},
$$

$$
\begin{aligned}
\lim_{L\to\infty}\bigg[&
\int_{|y|\le L}d^4y\,
[\delta_RM_1]_{b_H^2\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}e_1e_2}\\
&+\frac{841\pi^2g^2}{41580\rho^2}\log(L^2/\rho^2)\bigg]
=\frac{197929\pi^2g^2}{9313920\rho^2}.
\end{aligned}
$$

These values evaluate the connection insertion in the fixed source convention.

The full partner retains its overall $\hbar$ factor.

<a id="np-4c2ddf4548d3357c-1"></a>

### Equations for the full microscopic matching

The physical image of the source field $b_H$ includes an antifield:

$$
\Phi(b_H)=\frac{f_{++}^3}{2ig}
-\frac{\lambda_{3,\mathrm{old}}^{\ddagger-}}g.
$$

Here $\Phi$ is the source-to-physical insertion dictionary, and $\lambda^{\ddagger-}$ is the source BV antifield.

The elementary variations cancel exactly:

$$
Q\left(\frac{f_{++}^3}{2ig}\right)
=\frac i g\nabla_{+\dot\alpha}\bar\lambda^{3\dot\alpha},
\qquad
Q\left(-\frac{\lambda_{3,\mathrm{old}}^{\ddagger-}}g\right)
=-\frac i g\nabla_{+\dot\alpha}\bar\lambda^{3\dot\alpha}.
$$

Thus the BV completion cancels the elementary-$Qb$ contribution obtained from curvature alone.

The completed elementary variation contributes zero to the connection calculation.

The lowest remaining field insertion is

$$
(Q_{\rm physical}-Q_L)c_\alpha,
$$

This is the even normal gaugino variation with two fermions.

Its evaluation requires the collective-coordinate and body-color dependence of the complete Green family.

At the same background order, the matching includes the second-order response of $M_0$, the remaining normal and BV completion, action and measure responses, and composite/contact terms.

For example, the background expansion of the physical action $S$ gives

$$
\begin{aligned}
\frac{\partial^2}{\partial b_H^2}
\left(e^{-S}\Phi(\hbar M_1+\hbar^2M_0)\right)
=e^{-S}\bigg[&
\frac{\partial^2\Phi(\hbar M_1+\hbar^2M_0)}{\partial b_H^2}
-2\frac{\partial S}{\partial b_H}
\frac{\partial\Phi(\hbar M_1+\hbar^2M_0)}{\partial b_H}\\
&+\left(\left(\frac{\partial S}{\partial b_H}\right)^2
-\frac{\partial^2S}{\partial b_H^2}\right)\Phi(\hbar M_1+\hbar^2M_0)\bigg].
\end{aligned}
$$

At finite regulator, the difference between the full physical and collective Ward insertions is

<a id="np-b0509e2b7b2fc17a-1"></a>

$$
\int_\Gamma d\mu_0\,
\left[\Phi(Q_{\rm physical}\mathcal O)-Q_L\Phi(\mathcal O)\right].
$$

Here $\mathcal O$ is the source operator, and $\Phi$ evaluates its renormalized physical insertion.

The operator $Q_{\rm physical}$ is the full supersymmetry variation; $Q_L$ is its completed action on the retained collective coordinates.

Their difference measures field variations transverse to that family and the correction of the local insertion dictionary.

Matching means evaluating this difference together with the action, measure, boundary, and coincident-contact contributions in one prescription.

The physical scalar coefficient of the pairing is the coefficient of the resulting local operator identity; its evaluation remains open.

The cycle $\Gamma$ retains both size faces and the spatial face.

The insertion also retains derivatives with respect to the external even potentials.

Combining this residual with the boundary terms and common composite subtraction determines the physical scalar coefficient of the original $C_2$ pairing.

The value of this coefficient and its strong-coupling/weak-coupling instanton matching remain open problems.

<a id="np-a894ecb45e143187-1"></a>

### References for Sections 9–13

- Budzik–Gaiotto–Kulp–Williams–Wu–Yu, *Semi-Chiral Operators in 4d N=1 Gauge Theories*, especially the instanton-pairing proposal in §8: [arXiv:2306.01039](https://arxiv.org/abs/2306.01039).
- Witten, *Chiral Ring of Sp(N) and SO(N) Supersymmetric Gauge Theory in Four Dimensions*, §5, discussion of coincident and separated products: [hep-th/0302194](https://arxiv.org/abs/hep-th/0302194).
- Casher–Shamir, comparison of point splitting, Wilson-line variation, and finite-size boundaries: [hep-th/9908074](https://arxiv.org/abs/hep-th/9908074).
