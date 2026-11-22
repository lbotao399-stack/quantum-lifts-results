# Massless secondary class in the perturbative Ward complex

Status: the proposed secondary class is closed and nonexact over formal series with nonnegative powers of $\hbar$. The associated polynomial incoming class is $3\hbar g^2/(64\pi^2)$ times this class. Its nonexactness, the absence of $\hbar$-torsion, and its nonperturbative differential are not established by this audit.

Source: `/Users/libotao/np_lifts_engine/work/allloop_repo/main.tex`, especially `eq-B-15`, `eq-W-fibre-explicit`, `eq-antifield-loss`, and the spinor convention $\epsilon^{\dot+\dot-}=1=-\epsilon_{\dot+\dot-}$ near line 799. The coefficient and polynomial Higgs completion are checked in `fundamental_triangle_coefficient_audit.md`; the source-compatible antichiral variables and Euler cochain are specified in `microscopic_primitive_audit.md`. No manuscript was modified.

## 1. Operators and closure

Use the existing light-flavor Euler cochain $A$, heavy antichiral scalar contraction $\bar h=tv$, ghost matrix $\mathsf c_\gamma=\mathsf c^At_A^{\mathsf T}$, and Higgs-frame expression $K$ from the two preceding audits. The zero jet $\beta^{M^4{}_{\tilde4}}$ is the odd dual letter of the gauge-singlet heavy meson. Its equation is

$$
Q_0\beta^{M^4{}_{\tilde4}}
=-\frac{\bar h}{\bar\mu}-\bar m.
$$

Here $\bar\mu\ne0$ is the source superpotential normalization, and $g\ne0$ is held fixed. The expression $\bar hK$ is the polynomial

$$
\bar hK
=t_a(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b
(\partial_{+\dot\alpha}v)^b
-(\partial_{+\dot\alpha}t)_a
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)^a{}_b v^b.
$$

The indices $a,b$ are abstract color slots, while $\dot\alpha$ labels the source's two dotted spinor components; raising uses its stated spinor tensor. In particular, this expression is defined at $\bar h=0$ without inverting any field.

The proposed secondary operator is

$$
C=\bar\mu\beta^{M^4{}_{\tilde4}}
\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]
+\frac{i}{g}\bar hK.
$$

The quadratic ghost-jet trace is closed globally, not just after Cartan evaluation. Indeed $Q_0\mathsf c_\gamma=-ig\mathsf c_\gamma^2$, and the graded Leibniz rule gives

$$
\begin{aligned}
&Q_0\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]\\
&=-ig\operatorname{Tr}_{\mathrm{fund}}\left[
\mathsf c_\gamma(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)
-(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\mathsf c_\gamma
\right]=0.
\end{aligned}
$$

The two omitted intermediate terms cancel directly. In the displayed pair, moving the odd ghost past two odd derivatives gives sign $+1$.

Using $Q_0(\bar hK)=-ig\bar h\operatorname{Tr}_{\mathrm{fund}}[(\partial\mathsf c_\gamma)(\partial\mathsf c_\gamma)]$ with the same dotted contraction,

$$
\begin{aligned}
Q_0C
&=-(\bar h+\bar\mu\bar m)
\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]\\
&\quad+\bar h\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]\\
&=-\bar\mu\bar m\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right].
\end{aligned}
$$

For positive-loop operations, each graph edge consumes one $\beta$ or $\mathsf b$ terminal. At loop order $n$, a root with $r$ terminals has $r-n-\#(\bar W\text{ vertices})$ output terminals. The root $\bar hK$ and the ghost trace have none, so all their positive-loop images vanish. The other term in $C$ has one terminal, the singlet $\beta^{M^4{}_{\tilde4}}$. For $n\ge2$ counting excludes every graph. At $n=1$, a superpotential vertex is excluded by counting. Gauge-only graphs cannot consume this singlet terminal; treating it as a spectator leaves a root with no terminals, again excluded. Consequently

$$
\mathbf Q_{\mathrm{pert}}C
=-\bar\mu\bar m\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right],
\qquad
\left.\mathbf Q_{\mathrm{pert}}C\right|_{m=0}=0.
$$

The operator is relative to constant gauge transformations: the trace and $\bar hK$ are invariant under constant conjugation, the meson beta is a singlet, and contraction with a constant ghost shift annihilates every differentiated ghost. The use of $1/g$ is only a coefficient, not a Higgs-field localization. If regular dependence on $g$ at zero is required, $gC$ is polynomial in $g$ and obeys the multiplied identities.

## 2. A detector that includes every jet

Work at $m=0$ in the polynomial jet algebra. Evaluate all scalar jets at zero; all beta jets at zero except the zero jet $\beta^{M^4{}_{\tilde4}}$; and all $\mathsf b$ jets and the undifferentiated ghost at zero. Retain only

$$
\partial_{+\dot+}\mathsf c_\gamma=t_3^{\mathsf T}\xi_{\dot+},
\qquad
\partial_{+\dot-}\mathsf c_\gamma=t_3^{\mathsf T}\xi_{\dot-}.
$$

Here $t_3$ is a nonzero Cartan generator of the source $SU(2)$ fundamental representation, with its original normalization, and $\xi_{\dot+},\xi_{\dot-}$ are independent odd evaluation variables. All other ghost jets vanish. The retained beta is a third independent odd variable. The chosen Cartan generator has $\operatorname{Tr}_{\mathrm{fund}}(t_3^2)\ne0$; no numerical trace normalization is imposed.

This evaluation kills the $Q_0$-image of every jet generator:

1. Every differentiated $Q_0\gamma=-ig\mathsf c\gamma$ contains a scalar jet.
2. At zero mass every differentiated superpotential derivative in $Q_0\beta$ is quadratic in scalar jets. The gauge term of every charged beta vanishes because those beta jets vanish. The only retained beta belongs to a trivial representation, so no generator matrix can mix it into a charged beta equation.
3. Every differentiated $Q_0\mathsf b$ contains either an $\mathsf b$ jet or a scalar jet in its matter term.
4. For every pair of nonnegative jet orders $r,s$, the differentiated ghost equation is

$$
Q_0\bigl(\partial_{+\dot+}^{r}\partial_{+\dot-}^{s}\mathsf c^A\bigr)
=-\frac g2 C^A{}_{BC}
\partial_{+\dot+}^{r}\partial_{+\dot-}^{s}
(\mathsf c^B\mathsf c^C).
$$

The labels $A,B,C$ here are Lie-algebra components and $C^A{}_{BC}$ are the source structure constants. Any potentially surviving term after the full derivative expansion has both ghost factors in Cartan color $3$, hence is multiplied by $C^A{}_{33}=0$. This also handles second jets whose image contains two retained first jets, and arbitrary higher jets.

Because $Q_0$ is a derivation, the evaluation annihilates every $Q_0$-boundary in the full polynomial jet algebra. It therefore annihilates boundaries in the relative subcomplex as well; the evaluation need not itself preserve the full gauge group.

The polynomial $\bar hK$ evaluates to zero. With the source raising convention, the ordered surviving ghost contraction is

$$
\begin{aligned}
&\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right]\\
&\longmapsto
\operatorname{Tr}_{\mathrm{fund}}(t_3^2)
\bigl(\xi_{\dot+}\xi_{\dot-}-\xi_{\dot-}\xi_{\dot+}\bigr)
=2\operatorname{Tr}_{\mathrm{fund}}(t_3^2)
\xi_{\dot+}\xi_{\dot-}.
\end{aligned}
$$

Hence

$$
C\longmapsto
2\bar\mu\operatorname{Tr}_{\mathrm{fund}}(t_3^2)
\beta^{M^4{}_{\tilde4}}\xi_{\dot+}\xi_{\dot-}\ne0.
$$

Thus $C$ is not $Q_0$-exact. Over $\mathbb C[[\hbar]]$, meaning formal series with no negative powers of $\hbar$ and with the stated nonzero parameters fixed, any equation $C=\mathbf Q_{\mathrm{pert}}U$ would imply $C=Q_0U_0$ at order $\hbar^0$, where $U_0$ is the constant coefficient of $U$. This contradiction proves that $C$ is not $\mathbf Q_{\mathrm{pert}}$-exact.

## 3. The heavy/light cross term and the original incoming operator

The root $\beta^{M^4{}_{\tilde4}}A$ has two terminals. Its one-loop operation cannot be evaluated by assuming a quantum Leibniz rule. The possible contractions are instead as follows.

With gauge-only vertices, the singlet beta cannot contract and remains a spectator. The sign from the odd operation passing this odd spectator gives

$$
Q_1\bigl(\beta^{M^4{}_{\tilde4}}A\bigr)
=-\beta^{M^4{}_{\tilde4}}Q_1A
=\frac{3g^2}{64\pi^2}\beta^{M^4{}_{\tilde4}}
\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right].
$$

A one-loop graph containing a superpotential vertex has the triangle with one superpotential and one gauge vertex. To contract the heavy meson beta, the superpotential vertex must be the monomial $\overline M^4{}_{\tilde4}\,\overline q_4\,\overline{\widetilde q}^{\tilde4}/\bar\mu$. Its edge to the gauge vertex consumes that vertex's matter beta and fixes its flavor to $4$. The remaining edge from the gauge vertex to the root would then require a heavy charged beta in $A$. The charged betas of $A$ have only flavors $1,2,3$. Its meson betas are singlets and cannot contract with that gauge vertex. This exhausts the cross contractions. A second superpotential vertex has insufficient terminals and vanishes independently.

At $n=2$, counting allows only gauge vertices, which again cannot consume the heavy singlet beta. Removing that spectator leaves only one effective root terminal, so the graph vanishes. At $n\ge3$, the original root terminal count already excludes every graph. Thus the displayed one-loop formula is the complete positive-loop correction.

The classical graded Leibniz rule, $Q_0A=0$, and the checked spectator sign yield

$$
\begin{aligned}
\mathbf Q_{\mathrm{pert}}
\bigl(\bar\mu\beta^{M^4{}_{\tilde4}}A\bigr)
&=-(\bar h+\bar\mu\bar m)A\\
&\quad+\frac{3\hbar g^2}{64\pi^2}\bar\mu\beta^{M^4{}_{\tilde4}}
\operatorname{Tr}_{\mathrm{fund}}
\left[(\partial_{+\dot\alpha}\mathsf c_\gamma)
(\partial_+{}^{\dot\alpha}\mathsf c_\gamma)\right].
\end{aligned}
$$

Consequently the exact polynomial identity, valid also at nonzero mass, is

$$
-\mathbf Q_{\mathrm{pert}}
\bigl(\bar\mu\beta^{M^4{}_{\tilde4}}A\bigr)
+\frac{3\hbar g^2}{64\pi^2}C
=(\bar h+\bar\mu\bar m)A
+\frac{3i\hbar g}{64\pi^2}\bar hK.
$$

At $m=0$, writing brackets for perturbative cohomology classes, this gives

$$
\left[\bar hA+\frac{3i\hbar g}{64\pi^2}\bar hK\right]
=\frac{3\hbar g^2}{64\pi^2}[C].
$$

The proof $[C]\ne0$ does not prove $\hbar[C]\ne0$: a nonzero class can be annihilated by a positive power of $\hbar$, which is the meaning of $\hbar$-torsion here. The detector above annihilates $Q_0$-boundaries, not all positive-loop images, and cannot by itself rule this out. Neither inverting $\hbar$ nor specializing $\hbar=1$ is licensed by the formal leading-coefficient argument. The nonperturbative image of either operator still requires the local instanton/source matching identified in `microscopic_primitive_audit.md`.
