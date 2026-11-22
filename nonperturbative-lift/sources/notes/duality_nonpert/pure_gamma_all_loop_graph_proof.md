# Pure-gamma closure and zero-jet nonexactness

Checked 2026-11-13 against Budzik–Kulp, arXiv:2512.07771v3, §§2.3, 3: Bochner–Martinelli kernel eq. (2.31), shifted kernel eq. (2.36), master integral eq. (2.43), action eq. (3.1), and classical differential eqs. (3.4)–(3.5). The earlier equation references 56/62 are not the current v3 numbering. [Primary source](https://arxiv.org/pdf/2512.07771).

Current local source was also checked at `/Users/libotao/np_lifts_engine/work/allloop_repo/main.tex`: interaction `eq-B-10`, cut-edge kernels `eq-C-heat-spatial` / `eq-C-explicit-propagator`, antifield count `eq-antifield-loss`, finite-loop corollary `b-00235`, and derivative lemma `b-00094`. No project source was edited.

## 1. The proposed cubic half-edge count is correct

Consider one source containing only gamma letters. Let $r$ count contracted source letters; uncontracted source gamma spectators are initially omitted. Let $g,m,w$ count bulk vertices of types $bcc$, $\beta c\gamma$, and $\gamma^3$, respectively, and $V_b=g+m+w$. For a contributing connected 2-Laman anomaly graph, $E$ counts all edges, including the single cut edge, and $n$ is the loop number:

$$
E=2(V_b+1)-3=2V_b-1,
\qquad n=E-(V_b+1)+1=V_b-1.
$$

Each source edge consumes one bulk beta. If $o$ is the number of uncontracted bulk half-edges, cubic valence gives

$$
o=3V_b-(2E-r)=r-V_b+2.
$$

Write $b_{\mathrm{out}},c_{\mathrm{out}},\beta_{\mathrm{out}},\gamma_{\mathrm{out}}$ for nonnegative output-letter counts. The two propagator species give

$$
c_{\mathrm{out}}-b_{\mathrm{out}}=g+m,
\qquad
\gamma_{\mathrm{out}}-\beta_{\mathrm{out}}=r+3w.
$$

Therefore

$$
r-V_b+2=o
\ge g+m+r+3w=r+V_b+2w,
$$

$$
V_b+w\le1,
\qquad n+w\le0.
$$

Every $n\ge1$ contribution is impossible. Adding gamma spectators adds the same number to $o$ and to $\gamma_{\mathrm{out}}-\beta_{\mathrm{out}}$, and cancels from the inequality. No assumption about fan topology, planarity, color contraction, root degree, or Wick sign entered.

## 2. A shorter count removes the cubic restriction on W

Every $b$–$c$ or $\beta$–$\gamma$ edge consumes exactly one $b$ or $\beta$. Gauge/matter vertices supply one such letter each; a $W(\gamma)$ vertex supplies none, at any polynomial degree. The source supplies none. Consequently

$$
\#(b,\beta\text{ in output})
=(g+m)-E
=V_b-w-(2V_b-1)
=-n-w.
$$

Nonnegativity gives the same exclusion. In the source higher-bracket scheme,

$$
Q_nO(\gamma)=0\qquad(n\ge1).
$$

For gauge-invariant undifferentiated $O(\gamma)$, $Q_0\gamma=c\gamma$ gives $Q_0O=0$. Thus $Q_{\mathrm{pert}}O=0$ without adding any perturbative completion.

The cut edge carries the same field pairing as an uncut edge. An internal QME-defect cut is an edge joining two original interaction vertices, not an extra independent interaction type. Normal ordering removes self-contractions; allowing a self-contraction formally would still consume one $b$ or $\beta$ and would not evade this count. Vacuum components are absent from the connected operator anomaly. Gauge-anomaly cancellation is needed separately for $Q_{\mathrm{pert}}^2=0$ and hence for speaking of its cohomology.

## 3. Independent derivative-degree calculation

Use the primary source normalization of the Bochner–Martinelli form:

$$
P(z,\bar z)=\frac{\bar z^2d\bar z^1-\bar z^1d\bar z^2}{|z|^4}.
$$

Here $z^1,z^2$ are the two holomorphic coordinates; $t>0$ below is a uniform dilation parameter, distinct from edge Schwinger times. Under pullback by $z\mapsto tz$, including the differential forms,

$$
P(tz,t\bar z)=t^{-2}P(z,\bar z),
\qquad
(\bar\partial P)(tz,t\bar z)=t^{-2}\bar\partial P(z,\bar z),
\qquad
d^2(tz)=t^2d^2z.
$$

The regulated identity holds with each Schwinger time and cutoff rescaled by $t^2$. The cut heat kernel has the same weight $t^{-2}$; the existing local-source normalization factors are unchanged. The finite local UV anomaly has this homogeneous weight. There is no additional counterterm vertex in the cited UV-finite graph scheme.

For a graph $\Gamma$ with $V$ vertices, $E$ edges, and a fixed source vertex, the master kernel $I_\Gamma$ uses output momenta $\lambda$ and holomorphic input shifts $w$. Rescaling the $V-1$ integration variables gives

$$
I_\Gamma[\lambda/t;tw]
=t^{2(V-1)-2E}I_\Gamma[\lambda;w]
=t^{-2n}I_\Gamma[\lambda;w],
\qquad n=E-V+1.
$$

For a Taylor monomial with total shift degree $p$ and output-momentum degree $q$, this implies

$$
t^{p-q}=t^{-2n},\qquad q=p+2n.
$$

Input derivatives on contracted letters are the Taylor degree $p$; input derivatives on spectators remain on those spectators. Hence, for derivative-free interactions,

$$
\text{total output holomorphic derivatives}
=\text{total input holomorphic derivatives}+2n.
$$

This total-degree statement suffices; it does not require claiming a separate count for each coordinate. It applies to every 2-Laman topology. Holomorphic locality permits the kernel Taylor coefficients to be read as polynomial differential operators; no inverse derivatives are permitted in the local jet complex.

## 4. Projection excludes all perturbative primitives of a nonzero scalar class

Let $\pi_0$ denote setting every positive holomorphic jet to zero, while retaining all zero-jet fields. For derivative-free interactions $Q_0$ preserves total derivative degree, and the previous calculation yields

$$
\pi_0Q_0=Q_0\pi_0,
\qquad \pi_0Q_n=0\quad(n\ge1),
\qquad
\pi_0Q_{\mathrm{pert}}=Q_0\pi_0.
$$

Suppose the gauge-invariant scalar $O(\gamma)$ represents a nonzero class of the full zero-jet $Q_0$ complex. If a local formal or finite-Laurent perturbative primitive $U$ existed with

$$
Q_{\mathrm{pert}}U=\hbar^rO,
$$

then

$$
Q_0\pi_0U=\hbar^rO.
$$

The coefficient of $\hbar^r$ says that $O$ is $Q_0$-exact in the zero-jet complex, a contradiction. This proves nonexactness of $O$ and absence of $\hbar$-torsion for this class, also after allowing finitely many negative powers of $\hbar$.

The needed zero-jet scalar-boundary test is direct. The following displayed equality uses **only Budzik–Kulp v3 eq. (3.5)**, whose convention is $Q_0\beta_I=\beta_Ic+\partial W/\partial\gamma^I$. Here $\gamma^I$ and $\beta_I$ have dual matter-component indices, and $\partial_L/\partial\beta_I$ is the left derivative for an odd beta letter. This gives

$$
\left.Q_0U\right|_{b=c=\beta=0}
=\frac{\partial W}{\partial\gamma^I}
\left.\frac{\partial_LU}{\partial\beta_I}\right|_{b=c=\beta=0}.
$$

The original local manuscript instead uses a negative derivative of $\bar W$ in `eq-CY-completion-images`. Its scalar-extraction equality therefore carries that minus sign and its own antichiral potential. No identification of its $\bar W$ normalization with the primary source's $W$ is made here. The boundary-ideal membership statement is unchanged by that sign.

Every term from $Q_0c$, $Q_0\gamma$ or $Q_0b$ retains a $c$ or $\beta$ and vanishes in this extraction. Thus any pure-gamma scalar boundary lies in the source's F-term ideal, even when its primitive contains $b,c,\beta$. Invariant coefficients can be recovered by Reynolds averaging as in the separate $SU(5)$ proof. For a traceless adjoint $X$ with the stated F-term $F=X^2-(\operatorname{Tr}X^2/5)I$, the exact invariant-ideal certificate therefore excludes a zero-jet primitive for $(\operatorname{Tr}X^2)^2$. Additional fundamental gamma letters with no superpotential coupling cannot change this conclusion: set those letters to zero in the extracted identity.

No ring injection from the entire classical F-term quotient is asserted. A classical primitive of $(\operatorname{Tr}X^2)^3$ can acquire positive-jet Konishi terms under $Q_{\mathrm{pert}}$; the classical cubic relation, a quantum-product relation, and nonperturbative lifting remain distinct questions.

## Audit result

The proposed half-edge argument is valid. The stronger $\#(b,\beta)_{\mathrm{out}}=-n-w$ count proves closure without restricting $W$ to cubic degree. The independent scaling proof validates zero-jet projection on arbitrary derivative-input primitives. The conclusions require the specified local holomorphic-jet graph scheme and an anomaly-free differential; they do not cover nonlocal inverse derivatives, arbitrary extra BV vertices, or instanton corrections.
