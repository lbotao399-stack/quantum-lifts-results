<a id="np-section-12"></a>

## 14. Massive matter, the Konishi identity, and the vacuum algebra

For a chiral multiplet, the twisted matter fields are even $\gamma$ and odd $\beta$ in dual representations.

In the source convention, the tree differential contains

$$
Q_0\beta=\frac{\partial W}{\partial\gamma}+c\cdot\beta,
\qquad
Q_0\gamma=c\cdot\gamma.
$$

For a massive adjoint of $SU(2)$, the superpotential is $W=mu$, where $u$ is the quadratic invariant in the source metric:

$$
u=\frac12\gamma^i\gamma_i.
$$

For one fundamental flavor, the source uses two doublets $q^a,t^b$.

Here $a,b$ are abstract color slots and $\epsilon_{ab}$ is the invariant alternating form.

The meson and mass term are

$$
M=\epsilon_{ab}q^at^b,
\qquad W=mM.
$$

Assume $m\ne0$ below.

<a id="np-27b3f0c412c228c5-1"></a>

### Exact bilinear identity

The adjoint Dynkin index is $T=2$; the doublet pair has total index $T=1$.

In the source loop convention, the one-loop Konishi term is $\tfrac12\hbar T C_0$.

The expression $\beta_i\gamma^i$ is the canonical dual pairing; for the fundamentals, the two pairings are $\beta_q q$ and $\beta_t t$.

$$
Q_{\rm pert}(\beta_i\gamma^i)=2mu+\hbar C_0,
$$

$$
Q_{\rm pert}(\beta_q q+\beta_t t)=2mM+\frac{\hbar}{2}C_0.
$$

> **Proof.**
>
> The $Q_0$ gauge actions cancel in the canonical pairing.
>
> Euler's identity for the quadratic $W$ gives $2W$.
>
> The one-loop color trace gives $T/2$.
>
> The operator $Q_r$ acts on $r+1$ letters simultaneously, so its action on a bilinear vanishes for $r\ge2$.
>
> Hence, in cohomology,
>
> $$
> [u]=-\frac{\hbar}{2m}[C_0],
> \qquad [M]=-\frac{\hbar}{4m}[C_0].
> $$
>
> $\square$

<a id="np-7f0897b11f690436-1"></a>

### Holomorphic scale matching and two vacua

Denote the low-energy pure $SU(2)$ scale by $\Lambda_L$.

The holomorphic thresholds give

$$
\Lambda_L^6=m^2\Lambda^4
\quad\text{(massive adjoint)},
\qquad
\Lambda_L^6=m\Lambda^5
\quad\text{(one massive flavor)}.
$$

Assuming the exact low-energy chiral relation $[C_0]^2=\nu\Lambda_L^6$ gives

$$
[u]^2=\frac{\hbar^2\nu}{4m^2}\Lambda_L^6
=\frac{\hbar^2\nu}{4}\Lambda^4,
$$

$$
[M]^2=\frac{\hbar^2\nu}{16m^2}\Lambda_L^6
=\frac{\hbar^2\nu}{16}\frac{\Lambda^5}{m}.
$$

The first equation selects two vacua of the massive deformation of the $\mathcal N=2$ Coulomb branch.

Its source coefficient remains explicit for matching to the Seiberg–Witten normalization.

The second equation is the stationary equation of

$$
W_{\rm eff}=mM+\frac{\hbar^2\nu}{16}\frac{\Lambda^5}{M},
\qquad
\frac{\partial W_{\rm eff}}{\partial M}
=m-\frac{\hbar^2\nu}{16}\frac{\Lambda^5}{M^2}=0.
$$

The ADS coefficient is expressed in the same source convention for $C_0$.

Determining $\nu$ microscopically is the physical matching problem.

The source finite-sector page comparison includes $31$ sectors for the $SU(2)$ massive adjoint, $10$ for the massive fundamental, and $11$ for the $SU(3)$ massive adjoint.

The compared $E_1$ dimensions, $d_1,d_2$ ranks, and survivor counts agree with the pure theory, including the sectors of $C_0,C_2,C_4$ and the first secondary class.

These are finite-field census results; the preceding rational certificates establish the sectors specified there over the rational numbers.

<a id="np-2f6c6805ba565142-1"></a>

### The massless limit and finite vacuum algebra

At adjoint mass $m=0$, the bilinear identity becomes

$$
Q_0(\beta_i\gamma^i)=0,
\qquad
Q_1(\beta_i\gamma^i)=C_0.
$$

Thus $[C_0]$ becomes a boundary on the first loop page.

> **Calculation: the coefficient ring at zero mass.**
>
> At $m=0$, the exact bilinear identity and the vanishing of higher-loop actions on two letters give
>
> $$
> Q_{\rm pert}(\beta_i\gamma^i)=\hbar C_0,
> \qquad \hbar[C_0]=0.
> $$
>
> After adjoining $\hbar^{-1}$, the same equation has the primitive
>
> $$
> Q_{\rm pert}(\hbar^{-1}\beta_i\gamma^i)=C_0.
> $$
>
> The first-loop page records $d_1[\beta_i\gamma^i]=[C_0]$.
>
> In the nonnegative-power formal ring the exact equation is $\hbar[C_0]=0$; deciding whether $[C_0]$ itself vanishes there is the corresponding formal-cohomology question.
>
> $\square$

Within the computed range, $u,u^2,u^3$ survive.

Among the low-charge sectors with every specified incoming map computed, $37$ have survivors, compared with $3$ in the pure theory.

The massive identities apply on the open locus $m\ne0$.

At $m=0$, the Coulomb phase has its own complex and cohomology.

At fixed $\nu\Lambda^{3N}\ne0$, the vacuum decomposition of the gaugino subalgebra can be computed directly.

Let $s_k=\nu^{1/N}\Lambda^3e^{2\pi ik/N}$ be the $N$ distinct roots, with $k$ labeling vacua.

The corresponding idempotents are

$$
e_k=\frac1N\sum_{\substack{j=0,\ldots,N-1\\\mathrm{finite\ geometric\ series}}}
\left(\frac{C_0}{s_k}\right)^j.
$$

> **Proof.**
>
> At the root $s_l$, the geometric series gives
>
> $$
> e_k(s_l)=\frac1N
> \sum_{\substack{j=0,\ldots,N-1\\\mathrm{finite\ geometric\ series}}}
> e^{2\pi i(l-k)j/N}=\delta_{kl}.
> $$
>
> Since $C_0^N-\nu\Lambda^{3N}$ has distinct roots, evaluation is an algebra isomorphism.
>
> Thus
>
> $$
> \mathbb C[C_0]/(C_0^N-\nu\Lambda^{3N})
> =\bigoplus_{\substack{k=0,\ldots,N-1\\\mathrm{vacua}}}\mathbb C e_k,
> \qquad
> C_0e_k=s_ke_k,
> \qquad e_ke_l=\delta_{kl}e_k.
> $$
>
> $\square$

This is the gaugino vacuum algebra at fixed nonzero scale.

The NP pairing and microscopic matching provide the proposed connection from the full semichiral complex to the infrared cohomology of each vacuum.
