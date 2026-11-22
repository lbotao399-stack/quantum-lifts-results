<a id="np-section-08"></a>

## 10. One-instanton profiles and the regulated Ward identity

The following physical calculation uses the source convention of `physical_profiles`.

Let $x_0$ and $\rho>0$ be the instanton position and size, and set $y=x-x_0$.

A one-instanton configuration is a finite-action Euclidean gauge field in the sector of topological charge one.

Its collective coordinates describe continuous changes of the solution: position, size, color orientation, and fermion zero modes.

A zero mode solves the homogeneous linearized field equation.

The odd coordinates anticommute and square to zero; their Berezin integral extracts the coefficient of their ordered product.

The four odd collective coordinates are $\xi_1,\xi_2,\eta^{\dot1},\eta^{\dot2}$, in the stated Berezin integration order.

The Berezin integral of their ordered top monomial is one.

$$
\chi_\alpha=\xi_\alpha+y_{\alpha\dot\beta}\eta^{\dot\beta},
\qquad
f(y)=\frac{\rho^2}{(y^2+\rho^2)^2},
$$

$$
C_0(y)=3f(y)^2\chi_1\chi_2,
\qquad
x_2(y)=4if(y)^3\chi_1\chi_2\eta^{\dot2}.
$$

These are the evaluations of $C_0,x_2$ in the fixed physical dictionary.

The separated $C_0C_0$ Grassmann coefficient is $9$ in this convention.

The earlier displayed convention assigns it the coefficient $36$.

All physical-profile values below use the first convention, with the collective-measure coefficient $\mathcal N$ explicit.

<a id="np-bd7195b0c005f775-1"></a>

### Coincident and separated points

At one point, the first-jet gaugino profiles depend on the two odd combinations $\chi_1,\chi_2$, so

$$
C_0(y)^2=0.
$$

At two points with separation $r=x_1-x_2\ne0$, there are four independent odd combinations.

The separated-point integral evaluates both operator profiles on the same instanton and integrates over its possible centers and sizes.

The insertion points remain distinct while this integral is performed.

The squared lengths in the denominators are Euclidean; the factor $r^2$ comes from integrating the four fermion collective coordinates.

After removing the factor $\mathcal N\Lambda^6$, the position and size integral is

$$
\int d^4x_0\int_0^\infty d\rho\,
\frac{\rho^9r^2}{((x_1-x_0)^2+\rho^2)^4((x_2-x_0)^2+\rho^2)^4}
=\frac{\pi^2}{90}.
$$

> **Proof.**
>
> The denominator vectors are the displacements of the two insertion points from the instanton center.
>
> Introduce a Feynman parameter $\alpha$:
>
> $$
> \int\frac{d^4x_0}{(y_1^2+\rho^2)^4(y_2^2+\rho^2)^4}
> =\frac{10\pi^2}{3}\int_0^1
> \frac{\alpha^3(1-\alpha)^3\,d\alpha}
> {(\rho^2+\alpha(1-\alpha)r^2)^6}.
> $$
>
> $$
> \int_0^\infty\frac{\rho^9\,d\rho}
> {(\rho^2+\alpha(1-\alpha)r^2)^6}
> =\frac1{10\alpha(1-\alpha)r^2}.
> $$
>
> The result is
>
> $$
> \frac{\pi^2}{3}\int_0^1\alpha^2(1-\alpha)^2d\alpha
> =\frac{\pi^2}{90}.
> $$
>
> Multiplication by the Grassmann coefficient $9$ gives
>
> $$
> \langle C_0(x_1)C_0(x_2)\rangle_{1\text{-inst}}
> =\frac{\mathcal N\pi^2}{10}\Lambda^6.
> $$
>
> The separation cancels exactly.
>
> $\square$

This integrated correlator supplies a nonzero instanton contribution.

Its identification with the exact chain identity of Section 9 requires local-composite and Ward matching in the same subtraction prescription.

<a id="np-6d6750afea048494-1"></a>

### The collective vector, measure, and two size faces

The source supercharge satisfies $Q(gV)_{\alpha\dot1}=c_\alpha$.

Its action on the instanton collective coordinates is

$$
(Qx_0^1,Qx_0^2,Qx_0^3,Qx_0^4)
=\frac18(\xi_2,-i\xi_2,\xi_1,i\xi_1),
$$

$$
Q\rho=\frac{i\rho}{8}\eta^{\dot1},
\quad Q\xi_1=Q\xi_2=Q\eta^{\dot1}=0,
\quad Q\eta^{\dot2}=\frac i4\eta^{\dot1}\eta^{\dot2}.
$$

The orientation vector is $(\eta^{\dot2}u_1-i\eta^{\dot2}u_2+\eta^{\dot1}u_3)/4$.

The orientation generators $u_i$ obey $[u_i,u_j]=-\epsilon^k{}_{ij}u_k$ and preserve the normalized Haar measure $dU$.

The collective measure is

$$
d\mu_0=\mathcal N\rho\,d\rho\,d^4x_0\,dU\,d^2\xi\,d^2\eta.
$$

The supersymmetry transformation induces an odd vector field on the collective coordinates.

Its divergence measures the infinitesimal change of the integration density.

Zero divergence converts the integrated variation of an insertion into its outward flux through the regulated integration boundary.

This equality is the collective-coordinate Ward identity.

Both size faces, the spatial face, and the variation of any regulator enter the same identity.

> **Proof: nilpotency and divergence.**
>
> For example, the first orientation component satisfies
>
> $$
> Q\left(\frac{\eta^{\dot2}}4\right)
> -\left(-\frac{i\eta^{\dot2}}4\right)\frac{\eta^{\dot1}}4=0.
> $$
>
> The remaining orientation components, position, size, and odd coordinates also have vanishing $Q^2$.
>
> With left Grassmann derivatives, the measure divergence is
>
> $$
> \frac1\rho\partial_\rho(\rho Q\rho)
> +\partial_{\eta^{\dot2}}(Q\eta^{\dot2})
> =\frac i4\eta^{\dot1}-\frac i4\eta^{\dot1}=0.
> $$
>
> On the smooth profiles, $QC_0=Qx_2=0$.
>
> The flux through a constant-size face is
>
> $$
> [\rho(Q\rho)x_2]_{\xi_1\xi_2\eta^{\dot1}\eta^{\dot2}}
> =-\frac{\rho^2}{2}f^3,
> $$
>
> $$
> \mathcal N\int d^4x_0\,d^2\xi\,d^2\eta\,
> \rho(Q\rho)x_2=-\frac{\mathcal N\pi^2}{40}.
> $$
>
> On $\varepsilon\le\rho\le R$, the upper and lower faces cancel exactly:
>
> $$
> \langle Qx_2\rangle_{\varepsilon\le\rho\le R}
> =-\frac{\mathcal N\pi^2}{40}
> -\left(-\frac{\mathcal N\pi^2}{40}\right)=0.
> $$
>
> $\square$

The $Q$ variation of the explicit regulator $e^{-m^2\rho^2}$ with $m>0$ is

$$
Qe^{-m^2\rho^2}=-\frac{im^2\rho^2}{4}\eta^{\dot1}e^{-m^2\rho^2}.
$$

The integrated insertion of this variation is

$$
\int_{\varepsilon\le\rho\le R}d\mu_0\,
(Qe^{-m^2\rho^2})x_2
=\frac{\mathcal N\pi^2}{40}
(e^{-m^2\varepsilon^2}-e^{-m^2R^2}),
$$

This equals the regulated boundary flux.

Write the renormalized local representative as $[C_0^2]_{\rm ren}$.

A local product has coincident insertion points.

Its renormalized representative is defined by subtracting singular short-distance terms, with possible finite local counterterms.

The same representative must be used in the quantum-product equation and the supercharge equation.

Their two displayed scalar coefficients then change together under a finite subtraction.

In a common subtraction prescription,

$$
C_0\star C_0=[C_0^2]_{\rm ren}+\nu_{\rm prod}\Lambda^6,
\qquad
Qx_2=[C_0^2]_{\rm ren}+\nu_Q\Lambda^6.
$$

> **Proof: the scheme-invariant combination.**
>
> The finite redefinition $[C_0^2]_{\rm ren}\mapsto[C_0^2]_{\rm ren}+s\Lambda^6$ subtracts $s$ from both coefficients.
>
> Thus $\nu_{\rm prod}-\nu_Q$ is invariant.
>
> In a state satisfying the full Ward identity,
>
> $$
> \langle[C_0^2]_{\rm ren}\rangle=-\nu_Q\Lambda^6,
> \qquad
> \langle C_0\star C_0\rangle
> =(\nu_{\rm prod}-\nu_Q)\Lambda^6.
> $$
>
> The invariant difference $\nu_{\rm prod}-\nu_Q$ is the coefficient selected by this local Ward comparison.
>
> Identifying it with a separated correlator also uses the specified point-splitting limit, contact terms, and normalization of the operator dictionary.
>
> $\square$

<a id="np-f863edf6248e0361-1"></a>

### The exact straight-Wilson-line contribution

Point splitting defines a composite operator by separating its factors, transporting their color indices to a common point, and then prescribing the coincident limit.

The chosen path and the subtraction of singular terms are part of this definition.

At fixed separation $r>0$, move the second ghost jet of the primitive to $r e_4$ and transport it back to the origin with the adjoint transporter $U(0,r)$.

The Wilson line $U(0,r)$ is parallel transport in the adjoint color representation along the straight segment.

It expresses the shifted jet in the color frame at the origin, so the color contraction in the split operator is gauge invariant.

Supersymmetry varies the connection along this segment; inserting that variation in the transporter gives the Wilson contribution computed below.

Here $\operatorname{tr}U$ is the trace of the adjoint transport matrix.

Here $e_4$ is the fourth Euclidean coordinate vector.

The source color labels $i,j,k,l=1,2,3$ use the anti-Hermitian basis $t_i=-i\sigma_i/2$ and metric $\delta_{ij}$.

$$
\begin{aligned}
(x_2)_{\rm split}=-\frac16\epsilon_{ijk}\bigl[&
c_{01}^i(0)c_{01}^j(0)U^k{}_l(0,r)c_{20}^l(r)\\
&-2c_{10}^i(0)c_{01}^j(0)U^k{}_l(0,r)c_{11}^l(r)\\
&+c_{10}^i(0)c_{10}^j(0)U^k{}_l(0,r)c_{02}^l(r)\bigr].
\end{aligned}
$$

The source jet notation is $c_{rt}=\partial_1^r\partial_2^tc$.

The Grassmann coefficient of the second-jet variation is $-4\operatorname{tr}U/3$; that of the separated product is $9$.

The transporter contributes the variation

$$
QU(0,r)=\int_0^1ds\,
U(0,sr e_4)\operatorname{ad}(rQ(gV)_4(sr e_4))U(sr e_4,r).
$$

The two size-boundary limits are

$$
\lim_{\rho\to0}\mathrm{flux}=0,
\qquad
\lim_{\rho\to\infty}\mathrm{flux}=-\frac{\mathcal N\pi^2}{40}.
$$

Define the auxiliary average $\mathbb E$ with a positive scalar measure, retaining the field normalization:

$$
\mathbb E[F]=\frac{90}{\pi^2}
\int_0^\infty d\rho\int_{\mathbb R^4}d^4y\,
\frac{\rho^9r^2F}
{(y^2+\rho^2)^4((y+r e_4)^2+\rho^2)^4}.
$$

The exact Wilson insertion ratio is

$$
\frac{\text{Wilson insertion}}{\langle C_0(0)C_0(r)\rangle}
=-\frac14+\frac4{27}\mathbb E[\operatorname{tr}U].
$$

Write $y=(\boldsymbol y,v)$, $w=\sqrt{|\boldsymbol y|^2+\rho^2}$, and $t=|\boldsymbol y|/w$.

Then

$$
\operatorname{tr}U=1+2\cos(2t\Delta),
\qquad
\Delta=\arctan\frac{v+r}{w}-\arctan\frac vw.
$$

The variables $t$ and $\Delta$ are independent under this average.

Their exact moments are

$$
\mathbb E[t^2]=\frac3{13},
\qquad
\mathbb E[\cos\Delta]=\frac{36}{49},
\qquad
\mathbb E[\Delta^2]=\frac{100\pi^2}{693}-\frac{91}{108}.
$$

> **Proof: a strictly positive bound.**
>
> The inequality $\cos(2t\Delta)\ge1-2t^2\Delta^2$ gives
>
> $$
> \mathbb E[\operatorname{tr}U]
> \ge3-4\frac3{13}\left(\frac{100\pi^2}{693}-\frac{91}{108}\right)
> =\frac{34}{9}-\frac{400\pi^2}{3003}.
> $$
>
> Consequently,
>
> $$
> \frac{\text{Wilson insertion}}{\langle C_0(0)C_0(r)\rangle}
> \ge\frac{301}{972}-\frac{1600\pi^2}{81081}
> >\frac{5470949}{47675628}>0.
> $$
>
> The final rational bound uses $\pi<22/7$.
>
> $\square$

In terms of the standard sine and cosine integrals

$$
\operatorname{Si}(x)=\int_0^x\frac{\sin u}{u}\,du,
\qquad
\operatorname{Ci}(x)=-\int_x^\infty\frac{\cos u}{u}\,du
\quad(x>0)
$$

the Wilson average has the finite closed expression

$$
\begin{aligned}
\mathbb E[\operatorname{tr}U]={}&
6480[\log2+\operatorname{Ci}(4\pi)-\operatorname{Ci}(8\pi)]\\
&+\frac{10935}{8}[\log3+\operatorname{Ci}(2\pi)-\operatorname{Ci}(6\pi)]\\
&+\frac{2916}{\pi}[\operatorname{Si}(4\pi)+\operatorname{Si}(6\pi)
-\operatorname{Si}(2\pi)-\operatorname{Si}(8\pi)]
-\frac{23963}{4}.
\end{aligned}
$$

This is the exact integrated contribution for the specified straight-line prescription.

The full local BV differential coefficient is determined by the common Ward matching that includes this term.
