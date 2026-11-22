<a id="np-section-09"></a>

## 11. Background response and compatible collective coordinates

<a id="np-cfc7bbb8238456aa-1"></a>

### Opposite-curvature connection

Continue with the physical source coupling $g$.

The background $b_H=1$ corresponds to the physical curvature $f_{++}^3=2ig$.

For the connection $gV$, the fixed opposite-curvature coefficients are

$$
q^3{}_1=\frac{g^2}{2},
\qquad q^3{}_2=-\frac{ig^2}{2},
$$

Every remaining component vanishes.

Set $s=y^2$.

The adjoint harmonic matrix is $H=Y/(s+\rho^2)$, with $Y=sR$ and $R$ the regular-to-singular gauge rotation.

The source radial function $h$ satisfies

$$
s h''+\frac{2(s+2\rho^2)}{s+\rho^2}h'=\frac14,
\qquad
\frac{d}{ds}\left(\frac{s^4h'}{(s+\rho^2)^2}\right)
=\frac{s^3}{4(s+\rho^2)^2}.
$$

Its regular solution is $h(s)=\rho^2F(s/\rho^2)+h_0$, where $h_0$ is the homogeneous constant and

$$
F(t)=\frac t8-\frac{\log(1+t)}4
+\frac{5t+2}{8t^2}
-\frac{(3t^2+3t+1)\log(1+t)}{4t^3}
+\frac{11}{24}.
$$

$$
F(0)=0,
\quad F'(0)=\frac1{16},
\quad F(1)=\frac{35}{24}-2\log2,
\quad F'(1)=-2+3\log2,
\quad\lim_{t\to\infty}\frac{F(t)}t=\frac18.
$$

In the supplied Euclidean metric, $a,b$ are abstract cotangent slots, $i,j$ are color components, and $k$ labels opposite duality.

The exact connection is

$$
\delta(gV)^i_a=-2\bar\eta^k{}_{ab}D^b(H^i{}_j h q^j{}_k).
$$

Here $\bar\eta^k{}_{ab}$ is the source anti-self-dual 't Hooft tensor, and $D$ is the instanton covariant derivative.

Substitution verifies the background-gauge equation and the prescribed opposite curvature.

<a id="np-733664d1f7c344fe-1"></a>

### The regular Green response and its Grassmann-rank witness

A background changes the instanton gauge connection and hence the fermion equation.

The linear response is the first-order change of its fermion profile.

A Green response solves the resulting inhomogeneous equation with the stated center and infinity conditions.

A homogeneous solution is a response to zero forcing; zero modes with finite integrated squared norm are removed by orthogonal projection when fixing collective coordinates.

This fixes which part is field deformation and which part is motion within the instanton family.

The fermion response satisfies the linearized Weyl equation

$$
\bar\sigma^aD_a\delta c
=-\bar\sigma^a[\delta(gV)_a,c_{\rm inst}],
$$

The matrices $\sigma^a,\bar\sigma^a$ are the source Euclidean spinor bridge.

The color-angular polynomials decompose into three harmonic channels of degrees $l=1,3,2$, with color-angular eigenvalues $c=2,5,4$, respectively.

In the following table, $c$ labels an eigenvalue; in field expressions it denotes the ghost.

The radial equation is

$$
4sk''+4(l+2)k'
+\left[\frac{4c}{s+\rho^2}-\frac{8s}{(s+\rho^2)^2}\right]k=f.
$$

<a id="np-a6cb89ed04a182b4-1"></a>

| Channel | $(l,c)$ | regular homogeneous function $u$ | forcing $f$ |
|---|---|---|---|
| degree one | $(1,2)$ | $(s+3\rho^2)/(s+\rho^2)$ | $g^2\rho^4h/(s+\rho^2)^4+g^2\rho^2s h'/[3(s+\rho^2)^3]$ |
| degree three | $(3,5)$ | $1/(s+\rho^2)$ | $g^2\rho^2h'/(s+\rho^2)^3$ |
| $\eta^{\dot1}$ | $(2,4)$ | $1/(s+\rho^2)$ | $g^2\rho^4h/(s+\rho^2)^4$ |

$$
k(s)=-u(s)\int_s^\infty\frac{dt}{t^{l+2}u(t)^2}
\int_0^t\frac{dz}{4}z^{l+1}u(z)f(z).
$$

The inner integral starts at zero and enforces regularity at the instanton center.

The outer integral ends at infinity and selects the stated homogeneous component.

The table supplies $u$ and the forcing $f$ for each angular channel, so the formula determines its radial response $k$.

> **Proof.**
>
> Using the homogeneous equation for $u$ gives
>
> $$
> \frac{d}{ds}\left[s^{l+2}u^2\frac{d}{ds}\left(\frac ku\right)\right]
> =\frac14s^{l+1}uf.
> $$
>
> Two integrations give the displayed Green formula.
>
> At the origin, the inner integral vanishes to order at least $l+2$, and the outer integrand is finite.
>
> The outer integral converges at infinity.
>
> Its chosen endpoint fixes the coefficient of the additional growing auxiliary solution to zero.
>
> $\square$

Denote the three Green coefficients by $k_1,k_3,k_{\eta^{\dot1}}$.

Their exact limits are

$$
\begin{aligned}
\lim_{s\to\infty}sk_1&=-\frac{g^2\rho^2}{96},&
\lim_{s\to\infty}s^2k_1'&=\frac{g^2\rho^2}{96},&
\lim_{s\to\infty}s^3k_1''&=-\frac{g^2\rho^2}{48},\\
\lim_{s\to\infty}s^2k_3&=-\frac{g^2\rho^2}{32},&
\lim_{s\to\infty}s^3k_3'&=\frac{g^2\rho^2}{16},&
\lim_{s\to\infty}s^4k_3''&=-\frac{3g^2\rho^2}{16},\\
\lim_{s\to\infty}\frac{s^2k_{\eta^{\dot1}}}{\log(s/\rho^2)}&=-\frac{g^2\rho^4}{32},&
\lim_{s\to\infty}\frac{s^3k_{\eta^{\dot1}}'}{\log(s/\rho^2)}&=\frac{g^2\rho^4}{16},&
\lim_{s\to\infty}\frac{s^4k_{\eta^{\dot1}}''}{\log(s/\rho^2)}&=-\frac{3g^2\rho^4}{16}.
\end{aligned}
$$

Here $E,F,H$ denote source Chevalley color components, and $e_1,e_2$ are external fermion parameters.

At $y=(0,0,0,R)$, the six-row jet determinant of the full regular response satisfies

<a id="np-d147d35f585789c5-1"></a>

$$
\lim_{R\to\infty}R^{14}
\left.\frac{\partial}{\partial b_H}
\det\left[
\frac{\partial(c_{F,10},c_{E,01},c_{E,02},c_{F,20},c_{H,10},c_{H,01})}
{\partial(\xi_1,\xi_2,\eta^{\dot1},\eta^{\dot2},e_1,e_2)}
\right]\right|_{b_H=0}
=g^2\rho^8\ne0.
$$

This first-order determinant is therefore nonzero for sufficiently large finite $R$.

The corrected ghost jets span six Grassmann directions.

The corresponding radial limits of $h_0$ and normalizable homogeneous zero modes vanish, leaving this coefficient unchanged.

The background response thereby removes the original five-dimensional-span obstruction.

<a id="np-835db4476c0bd4be-1"></a>

### External fermions, orthogonal projection, and even boundary variables

The external source uses $c_{\rm ext}=2i\lambda_{\rm ext}$.

The squared norms of the supersymmetric and superconformal zero modes in a centered ball of radius $L$ are, respectively,

$$
\frac{\pi^2L^4(L^2+3\rho^2)}{2(L^2+\rho^2)^3},
\qquad
\frac{\pi^2\rho^2L^6}{(L^2+\rho^2)^3}.
$$

All cross inner products vanish.

The four nonzero whole-space projections of the mixed response are

$$
\xi_2e_1\longmapsto-\frac12\eta^{\dot2},
\quad
\xi_1e_2\longmapsto\frac12\eta^{\dot2},
\quad
\eta^{\dot1}e_\alpha\longmapsto\frac{\rho^2}{2}\xi_\alpha.
$$

The first two finite-ball coefficients have magnitude $(L^2+2\rho^2)/(2(L^2+\rho^2))$; the last two have magnitude $\rho^2L^4/[2(L^2+\rho^2)(L^2+3\rho^2)]$.

All $32$ overlaps converge absolutely.

The squared norms of the four projected $\xi e$ responses are

$$
\frac{7\pi^2\rho^2}{30},\quad
\frac{\pi^2\rho^2}{6},\quad
\frac{\pi^2\rho^2}{6},\quad
\frac{7\pi^2\rho^2}{30}.
$$

For each $\eta^{\dot\alpha}e_\beta$ response,

$$
\lim_{L\to\infty}
\frac{\|\text{response}\|_{|y|\le L}^{2}}{\log(L/\rho)}
=2\pi^2\rho^4.
$$

The finite overlaps fix orthogonality; the response norm has the displayed logarithmic divergence.

The Cartan external connection response decomposes into a harmonic gauge transformation and a translation zero mode.

Its translation coefficients are

$$
(Qx_0^1,Qx_0^2,Qx_0^3,Qx_0^4)_{\rm ext}
=\frac{\rho^2}{4}(-ie_2,-e_2,-ie_1,e_1).
$$

A general body-frame color $n^i$ satisfies $n^in_i=1$ in the supplied color metric, and transforms as $Qn^i=\epsilon^i{}_{jk}\omega^jn^k$.

The additional fermion term in the complete first-order collective vector is

$$
(Q\xi_\alpha)_{\rm ext}
=-\frac{\rho^2}{2}
[n^3\eta^{\dot1}+(n^1-in^2)\eta^{\dot2}]e_\alpha.
$$

These additional first-order terms integrate to zero in the partner component with four instanton fermions and two external fermions.

The $Q\xi$ term vanishes through repeated $\eta$ factors; the translation term is a position derivative whose boundary integral vanishes.

Second-order closure includes the even boundary potentials:

$$
Q(gV)_{1\dot1,\infty}=2ie_1,
\qquad Q(gV)_{2\dot1,\infty}=2ie_2,
\qquad Qe_1=Qe_2=0.
$$

The orientation coefficient uses the full combination

$$
\eta^{\dot2}
+\rho^2[(gV)_{1\dot1,\infty}e_2-(gV)_{2\dot1,\infty}e_1],
$$

and

$$
Q\eta^{\dot2}
=\frac i4\eta^{\dot1}\eta^{\dot2}-4i\rho^2e_1e_2.
$$

> **Proof.**
>
> The boundary variables give
>
> $$
> Q[(gV)_{1\dot1,\infty}e_2-(gV)_{2\dot1,\infty}e_1]
> =4ie_1e_2.
> $$
>
> The background derivative contributes $\rho^2(i,1,0)e_1e_2$, canceling the orientation square $\rho^2(-i,-1,0)e_1e_2$ on the zero-potential slice.
>
> Evaluating $Q^2$ on position, size, the four instanton fermions, body-frame color, boundary potentials, and the framing Maurer–Cartan relation gives zero in every case.
>
> $\square$

On a gauge-orbit representative, the boundary-potential derivative acts on covariant jets by the same gauge transformation.

The gauge invariance of the partner polynomials makes this insertion vanish pointwise.

The difference for the centered-ball orthogonal representative is a translation whose integrated boundary term also vanishes.

These calculations establish collective-coordinate closure and the vanishing of the specified gauge contributions.
