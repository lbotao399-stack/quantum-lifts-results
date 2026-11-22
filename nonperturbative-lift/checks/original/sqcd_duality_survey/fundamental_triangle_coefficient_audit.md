# Fundamental triangle coefficient and polynomial Higgs completion

Authoritative input: current /Users/libotao/np_lifts_engine/work/allloop_repo/main.tex, read without edits. Anchors: eq-B-15 (line 1246), eq-W-fibre-explicit (1303), eq-D-lambda (1297), eq-B-10 (1423), eq-A-20 (1444), eq-1-1 (1518), eq-1-7 (1810), eq-1-25 (2201), and eq-2-26 (4423). Equation labels are authoritative.

## 1. Direct fundamental triangle contraction

The source uses even matter scalars $\gamma_i$, odd dual matter letters $\beta^i$, and odd ghosts $\mathsf c^A$. Its gauge-matter interaction and matter pairing are

$$
\mathcal I_{\mathrm{matter}}
=-ig\,\mathsf c^A(t_A)^j{}_i\gamma_j\beta^i,
$$

$$
\Delta_{vw}^{\mathrm{matter}}
=-\frac{\partial_L}{\partial\beta_v^i}
 \frac{\partial}{\partial\gamma_{i,w}}
-\frac{\partial}{\partial\gamma_{i,v}}
 \frac{\partial_L}{\partial\beta_w^i}.
$$

Here $i,j$ are component labels in the matter representation $R$, and $A$ is a Lie-algebra component label. $\partial_L$ is the left derivative of an odd variable. The observable is $\gamma_i\beta^i$.

There are two oriented matter cycles for ordered triangle edges $(12,13,23)$:

| Cycle | Edge 12 | Edge 13 | Edge 23 | Derivative sign including three pairing minus signs |
|---|---|---|---|---|
| First | $\beta_1,\gamma_2$ | $\gamma_1,\beta_3$ | $\beta_2,\gamma_3$ | $+1$ |
| Second | $\gamma_1,\beta_2$ | $\beta_1,\gamma_3$ | $\gamma_2,\beta_3$ | $+1$ |

The vertex numbers in this table label fields at vertices, not color slots. Apply edges from right to left, and the two derivatives in each pairing from right to left, as prescribed after eq-W-fibre-explicit. The initial odd order is
$\mathsf c_1\beta_1\mathsf c_2\beta_2\beta_3$.
Both cycles leave $\mathsf c_1^D\mathsf c_2^E$ with positive derivative sign. The two vertex coefficients give $(-ig)^2=-g^2$. The fiber contraction is

$$
-g^2\left[
\operatorname{Tr}_R(t_Dt_E)+
\operatorname{Tr}_R(t_Et_D)\right]
\mathsf c_1^D\mathsf c_2^E.
$$

Only the numerical generator matrices were cyclically permuted. No odd ghost was moved past another odd ghost in this trace step.

The kernel normalization, simplex integral and momentum conversion are

$$
\{\mathcal I,\mathcal I,\mathcal O\}_{\triangle}
=-\frac1{8\pi^2}\,
\mathcal I_\triangle\ \text{acting on the fiber contraction},
$$

$$
\mathcal I_\triangle(\lambda;0)
=(\lambda_{1,+\dot\alpha}\lambda_{2,+}{}^{\dot\alpha})
\int_0^1dr_1\int_0^{1-r_1}dr_2
=\frac12\lambda_{1,+\dot\alpha}\lambda_{2,+}{}^{\dot\alpha},
$$

$$
\lambda_{v,+\dot\alpha}
\longmapsto-\frac12\partial_{v,+\dot\alpha},
\qquad
Q_1\mathcal O=\frac1{2!}\{\mathcal I,\mathcal I,\mathcal O\}_0.
$$

Thus the complete multiplier is

$$
\frac12
\left(-\frac1{8\pi^2}\right)
\frac12\frac14(-2g^2)
=\frac{g^2}{64\pi^2}.
$$

The direct result for the actual representation is

$$
\boxed{
Q_1(\gamma_i\beta^i)
=
\frac{g^2}{64\pi^2}
\operatorname{Tr}_R(t_Dt_E)
(\partial_{+\dot\alpha}\mathsf c^D)
(\partial_+{}^{\dot\alpha}\mathsf c^E).
}
$$

The dotted contraction uses the manuscript's spinor metric. The contracted odd bilinear is symmetric in $D,E$:

$$
(\partial_{+\dot\alpha}\mathsf c^D)
(\partial_+{}^{\dot\alpha}\mathsf c^E)
=
(\partial_{+\dot\alpha}\mathsf c^E)
(\partial_+{}^{\dot\alpha}\mathsf c^D).
$$

The odd-field interchange and the spinor-metric interchange each supply a minus sign. No additional factor two remains. The adjoint equation eq-2-26 is only a check: with $(t_D^{\mathrm{adj}})^B{}_C=iC^B{}_{DC}$, its contracted color tensor is $\operatorname{Tr}_{\mathrm{adj}}(t_Dt_E)$. No fundamental Dynkin index or Pauli-matrix normalization was inserted.

## 2. Lock the antichiral matrix action

For a column of the source's antichiral scalars, define $\mathsf c_\gamma$ by its actual action,

$$
(\mathsf c_\gamma)_i{}^j
=\mathsf c^A(t_A)^j{}_i,
\qquad
Q_0\gamma=-ig\,\mathsf c_\gamma\gamma.
$$

These are the transpose matrices $t_A^{\mathsf T}$. Since the original matrices obey $[t_A,t_B]=iC^C{}_{AB}t_C$,

$$
[t_A^{\mathsf T},t_B^{\mathsf T}]
=-iC^C{}_{AB}t_C^{\mathsf T},
\qquad
Q_0\mathsf c_\gamma=-ig\,\mathsf c_\gamma^2.
$$

The quadratic trace is unchanged:

$$
\operatorname{Tr}(t_D^{\mathsf T}t_E^{\mathsf T})
=\operatorname{Tr}(t_Et_D)
=\operatorname{Tr}(t_Dt_E).
$$

The candidate $A$ is the Euler odd cochain $M\beta_M-q_{\mathrm{light}}\beta_q$, with all three light meson rows and all four meson columns included. The light quarks give three fundamental copies; the mesons are singlets. Hence

$$
\boxed{
Q_1A
=-\frac{3g^2}{64\pi^2}
\operatorname{Tr}_{\mathrm{fund}}
\left[
(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)
\right].
}
$$

The source convention $Q_0\beta=-\partial\bar W/\partial\gamma$ gives $Q_0A=0$ by the opposite Euler weights of the matter and meson factors in the parent superpotential.

If one instead uses the original matrix $\mathsf c^At_A$, its matrix BRST law has the opposite displayed sign. An epsilon-converted Higgs column also has the opposite displayed action. The correction must be transformed consistently.

## 3. Graded trace and the correction sign

Let $v$ be the heavy antichiral fundamental column, $t$ its heavy antifundamental row, and $h=tv$:

$$
Q_0v=-ig\,\mathsf c_\gamma v,
\qquad
Q_0t=ig\,t\mathsf c_\gamma,
\qquad Q_0h=0.
$$

In a color basis the GL(2) Higgs frame is

$$
G=
\begin{pmatrix}
v_1&-t_2\\
v_2&t_1
\end{pmatrix},
\qquad
\det G=h,
\qquad
Q_0G=-ig\,\mathsf c_\gamma G.
$$

The last equality uses $\operatorname{Tr}\mathsf c_\gamma=0$. No determinant-one normalization is needed. For $h\ne0$,

$$
K=
\operatorname{Tr}_{\mathrm{fund}}
\left[
(\partial_{+\dot\alpha}G)G^{-1}
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)
\right].
$$

The frame and its derivative are even; $\mathsf c_\gamma$ and its derivative are odd. Direct differentiation gives

$$
\begin{aligned}
Q_0[(\partial G)G^{-1}]
&=-ig\,\partial\mathsf c_\gamma
-ig\,\mathsf c_\gamma(\partial G)G^{-1}
+ig\,(\partial G)G^{-1}\mathsf c_\gamma,\\
Q_0(\partial\mathsf c_\gamma)
&=-ig\left[
(\partial\mathsf c_\gamma)\mathsf c_\gamma
+\mathsf c_\gamma(\partial\mathsf c_\gamma)\right].
\end{aligned}
$$

Here each derivative has the dotted slot shown in $K$. The cubic terms cancel using

$$
\operatorname{Tr}\left[
\mathsf c_\gamma(\partial G)G^{-1}
(\partial\mathsf c_\gamma)\right]
=
-\operatorname{Tr}\left[
(\partial G)G^{-1}
(\partial\mathsf c_\gamma)\mathsf c_\gamma\right].
$$

Therefore

$$
\boxed{
Q_0K
=-ig\,\operatorname{Tr}_{\mathrm{fund}}
\left[
(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right],
}
$$

and the numerical cancellation is

$$
-\frac{3g^2}{64\pi^2}
+\frac{3ig}{64\pi^2}(-ig)=0.
$$

## 4. Polynomial completion and all-order contraction check

Since $G^{-1}=\operatorname{adj}(G)/h$,

$$
hK=
\operatorname{Tr}_{\mathrm{fund}}
\left[
(\partial_{+\dot\alpha}G)\operatorname{adj}(G)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]
$$

is polynomial. Explicit multiplication gives

$$
\boxed{
hK=
t(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)
(\partial_{+\dot\alpha}v)
-
(\partial_{+\dot\alpha}t)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)v.
}
$$

Thus no inverse, square root, field equation, or extra normalization occurs. The equality $Q_0(hK)=-ig\,h\operatorname{Tr}(\partial\mathsf c_\gamma\partial\mathsf c_\gamma)$ extends across $h=0$ as a polynomial identity.

For $Q_1(hA)$, both triangle interaction vertices must be gauge-matter vertices: the root has no $\mathsf b$ or $\mathsf c$ for a gauge-vertex spoke. If its $\beta$ has light-quark flavor, the internal matter line forces the other root contraction to have that same light flavor. The two heavy scalars in $h$ are untouched. If the root $\beta$ is a meson, the necessary gauge vertex vanishes. Hence

$$
Q_1(hA)=hQ_1A.
$$

This is a contraction result, not a Leibniz assumption for $Q_1$.

Each edge consumes a $\beta$ or $\mathsf b$ terminal. An $n$-loop graph has $n+1$ interaction vertices and $2n+1$ edges. Each interaction supplies at most one terminal. The root $hA$ supplies one, so for $n\ge2$ the deficit is

$$
(n+1)+1-(2n+1)=1-n<0.
$$

The root $hK$ supplies none, giving

$$
(n+1)-(2n+1)=-n<0
\qquad(n\ge1).
$$

Consequently,

$$
Q_n(hA)=0\quad(n\ge2),
\qquad
Q_n(hK)=0\quad(n\ge1).
$$

With the manuscript's $\mathbf Q=\sum_{\substack{n\ge0\\\mathrm{loop\ orders}}}\hbar^nQ_n$,

$$
\boxed{
\mathbf Q\left[
hA+\frac{3i\hbar g}{64\pi^2}hK
\right]=0
}
$$

to all perturbative orders.

This audit verifies closure and its coefficient in the manuscript's twisted complex. It does not alone prove nonexactness of the odd cocycle, compute its instanton differential, or match it to a renormalized physical composite current.
