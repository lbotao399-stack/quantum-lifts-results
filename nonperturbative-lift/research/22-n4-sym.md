<a id="np-section-20"></a>

## 22. $\mathcal N=4$ SYM: known instanton effects and proved protection

<a id="np-fc53eb8790b564e6-1"></a>

### Exact vacuum detector for the pure-scalar $SU(N)$ sector

The pure-scalar sector consists of gauge-invariant polynomials in the three complex scalars, with no fermions or spacetime derivatives.

The quotient by their F-relations retains polynomial operator relations, including possible nilpotents.

Its reducedness means that every nonzero class has some nonzero value on an ordinary vacuum.

This property lets exact supersymmetric vacua test whether a scalar operator can become a $Q$-variation.

The three traceless complex adjoint scalars $X,Y,Z$ obey the F-relations

$$
[X,Y]=[Y,Z]=[Z,X]=0.
$$

The commuting-matrix scheme here means imposing each matrix entry of the three commutators as a polynomial relation.

The invariant quotient keeps the classes unchanged by simultaneous conjugation.

Multisymmetric functions are polynomials in the diagonal triples that remain unchanged when all triples are permuted together.

The theorem identifies these two explicitly described algebras, including any nilpotent classes of the original quotient.

In characteristic zero, Vaccarino's invariant-theory theorem identifies the full invariant quotient of the commuting-matrix scheme with multisymmetric functions.

Let $x_i,y_i,z_i$ be the simultaneous diagonal entries, with eigenvalue label $i$.

The group $S_N$ simultaneously permutes the $N$ triples.

Separating the three trace coordinates gives

$$
\left(\frac{\mathbb C[X,Y,Z]}{([X,Y],[Y,Z],[Z,X])}\right)^{SL(N,\mathbb C)}
\cong
\left(
\frac{\mathbb C[x_1,\ldots,x_N,y_1,\ldots,y_N,z_1,\ldots,z_N]}
{\left(
\sum_{\substack{i=1,\ldots,N\\\mathrm{eigenvalues}}}x_i,
\sum_{\substack{i=1,\ldots,N\\\mathrm{eigenvalues}}}y_i,
\sum_{\substack{i=1,\ldots,N\\\mathrm{eigenvalues}}}z_i
\right)}
\right)^{S_N}.
$$

The right-hand side is reduced, so every nonzero scalar class has a nonzero value on some diagonal vacuum.

[Vaccarino, Theorem 3](https://arxiv.org/abs/math/0602660) The exact $\mathcal N=4$ vacuum geometry preserves these vacua.

[Argyres–Bourget–Martone, §2.3](https://arxiv.org/abs/1904.10969)

> **Proof: exclusion of instantonic exactness.**
>
> Let $O$ represent the scalar class and be closed under the exact supercharge, $v$ be a diagonal vacuum where its value is nonzero, and $V$ be a candidate local primitive.
>
> $$
> O=[Q,V\}
> \quad\Longrightarrow\quad
> \langle O\rangle_v=\langle[Q,V\}\rangle_v=0,
> $$
>
> This contradicts $\langle O\rangle_v\ne0$.
>
> Every nonzero class of the pure-scalar quotient therefore remains nonzero in exact cohomology.
>
> $\square$

This argument covers the surviving scalar half-, quarter-, and one-eighth-BPS sectors.

The general fermion/derivative one-sixteenth-BPS complex requires a separate calculation.

<a id="n4-descendant-detector"></a>

**Proposition 22.1.** Consider the half-BPS scalar primaries built from one traceless complex scalar $Z$ in unitary $SU(N)$ theory. At scalar degree $p\ge2$ and total holomorphic translation order $k\ge0$, their *total translation descendants* span an exact cohomology subspace of dimension

$$
(k+1)[t^p]\prod_{\substack{j=2,\ldots,N\\\mathrm{invariant\ degrees}}}
\frac{1}{1-t^j}.
$$

Here $[t^p]$ extracts the coefficient of $t^p$. A total translation acts on the entire gauge-invariant operator.

> **Proof.** Symmetric polynomials in the eigenvalues of one traceless matrix are freely generated in degrees $2,\ldots,N$. The preceding vacuum detector and half-BPS shortening supply the independent scalar primaries at dimension $p$. In the radial-quantization convention of Chang and Yin, $E$ is dilation energy, $J$ is the left rotation Cartan charge and $H_1,H_2,H_3$ are the three internal Cartan charges. Their algebra gives
>
> $$
> 2\{Q,Q^\dagger\}=E-2J-H_1-H_2-H_3 .
> $$
>
> Each of the two holomorphic translations raises $E$ by $1$ and $J$ by $1/2$, with unchanged internal charges. Every translated state still saturates the bound. Positivity therefore makes it annihilated by both $Q$ and $Q^\dagger$. A harmonic state equal to $Q$ applied to another state has zero norm.
>
> Let $G_{IJ}$ be the positive definite two-point matrix of the independent degree-$p$ primaries, labeled by $I,J$ in this equation. Their separated two-point function is $G_{IJ}/(x^2)^p$ with $x^2=x_a x^a$ in the Euclidean metric. Take independent complex translation vectors $u^a,v^a$ in the holomorphic plane. The complexified metric pairs any two such vectors to zero. With $\partial^{(y)}_a$ the flat derivative in the second insertion point, for $0\le r\le k$ the translated correlator is exactly
>
> $$
> \left.
> (u^a\partial^{(y)}_a)^r(v^b\partial^{(y)}_b)^{k-r}
> \frac{G_{IJ}}{((x-y)^2)^p}\right|_{y=0}
> =
> \frac{2^k(p+k-1)!}{(p-1)!}
> \frac{G_{IJ}(x_a u^a)^r(x_b v^b)^{k-r}}{(x^2)^{p+k}} .
> $$
>
> The $k+1$ displayed monomials are linearly independent as polynomials in $x$. Invertibility of $G$ therefore rules out any zero state in this descendant space. Harmonicity excludes exactness. For $SU(2)$ and $p=2$, $\operatorname{Tr}Z^2$ alone gives exactly $k+1$ independent protected descendants at every order $k$. $\square$

This uses the superconformal representation and the exact form of a scalar two-point function, allowing the positive matrix $G$ to depend on the coupling. It identifies a derivative subspace that can be removed from a search for a nonzero instanton pairing. Spinful primaries and derivatives distributed among separate factors beyond this subspace require their own tests. [Chang–Yin, §2](https://arxiv.org/abs/1305.6314) · [Córdova–Dumitrescu–Intriligator](https://arxiv.org/abs/1612.00809)

Representation theory establishes absolute protection for a narrower range.

In Córdova–Dumitrescu–Intriligator notation, $R_1,R_2,R_3$ are $SU(4)_R$ Dynkin labels, $[j;\bar j]$ are Lorentz labels, and $B_1$ denotes BPS shortening:

$$
B_1\overline B_1[0;0]_{R_1+R_2+R_3}^{(R_1,R_2,R_3)},
\qquad R_1=R_3\in\{0,1\}
$$

These are the absolutely protected scalar families.

For example, a long scalar multiplet reaching $\Delta=2$ admits the following recombination.

Here $L$ denotes the long type, and $A_2$ the scalar semishort type:

$$
\begin{aligned}
L\overline L[0;0]_{\Delta}^{(0,0,0)}\longrightarrow{}&
A_2\overline A_2[0;0]_2^{(0,0,0)}
\oplus B_1\overline A_2[0;0]_3^{(2,0,0)}\\
&\oplus A_2\overline B_1[0;0]_3^{(0,0,2)}
\oplus B_1\overline B_1[0;0]_4^{(2,0,2)}.
\end{aligned}
$$

The vacuum detector supplies additional dynamical protection for the scalar quotient.

[Córdova–Dumitrescu–Intriligator, §2.2.4](https://arxiv.org/abs/1612.00809)

Known nonperturbative correlator theorems cover half-BPS three-point and extremal correlators, and three-point functions with one quarter-BPS and two half-BPS operators.

[Baggio–de Boer–Papadodimas](https://arxiv.org/abs/1203.1036) For external labels $[q_i,p_i,q_i]$, the extremal/next-to-extremal condition of the quarter-BPS four-point theorem is

$$
p_1+2q_1=p_2+2q_2+p_3+2q_3+p_4+2q_4-2k,
\qquad k\in\{0,1\}.
$$

The inputs are protected analytic multiplets.

These correlator statements address a separate question from classifying the full one-supercharge spectrum.

[Bissi–Fardelli–Puerta-Ramisa](https://arxiv.org/abs/2610.03606) Flatness of Higgs-branch primary bundles assumes the local absence of level crossing.

Comparison with the formal expansion at the weak-coupling cusp remains a separate step.

[Niarchos](https://arxiv.org/abs/1807.04296) For every nonnegative integer $k$,

$$
\log(g^{-k}e^{-8\pi^2/g^2})=-k\log g-8\pi^2/g^2\longrightarrow-\infty,
\qquad
\lim_{g\to0^+}g^{-k}e^{-8\pi^2/g^2}=0.
$$

Formal perturbative vanishing and finite-coupling vanishing therefore pose different analyticity questions.

<a id="np-87f770bad4f64cbb-1"></a>

### Perturbative nature of the known cohomology pairings

The two operators in the classical $SO(7)$–$USp(6)$ mismatch have Chang–Lin charges

$$
(J_1,J_2,q_1,q_2,q_3)
=(1/2,1/2,5/2,5/2,5/2),
\qquad(0,0,3,3,3).
$$

The Choi–Lee generalized Konishi correction gives

$$
Q_KO_f=O_m+Q_0\Lambda,
\qquad [Q_KO_f]=[O_m]\ne0.
$$

Here $O_f,O_m$ are their fortuitous and monotone classes, respectively, and $\Lambda$ is the source classical primitive.

This is a pairing at first quantum order.

[Chang–Lin](https://arxiv.org/abs/2510.24008); [Choi–Lee](https://arxiv.org/abs/2511.09519) In the subsequent source normalization, the exact detector coefficients are

$$
\rho_J(Q_1F_4)=-4608\kappa_Q\Omega,
\qquad \rho_J(\mathcal S^4)=384\Omega,
\qquad
[Q_1][F_4]=-12\kappa_Q[\mathcal S^4]\ne0.
$$

Here $F_4$ is the fortuitous representative used in that paper, $\rho_J$ is its root-direction restriction, $\Omega$ is the ordered product of the eight retained Grassmann variables, $\kappa_Q$ is the pair-action coefficient, and $\mathcal S$ is its gaugino bilinear.

[Kibe–van Leuven–Mathieson](https://arxiv.org/abs/2610.22332) Established protection of the $SU(3)$ brane-fused towers currently extends through first quantum order.

[Kim–Lee–Lee–Song](https://arxiv.org/abs/2610.20340)

The full one-loop supercharge construction also retains the chain-level relation

$$
Q_1^2+\{Q_0,Q_2\}=0.
$$

Thus $Q_1$ induces a differential on $H(Q_0)$.

Its chain-space square is controlled by the displayed $Q_2$ term.

[Budzik–Kulp](https://arxiv.org/abs/2512.07771)

<a id="np-516d83e0f53cc864-1"></a>

### Nonzero instanton coefficient in a long multiplet

A supersymmetry multiplet groups operators related by the supercharges.

A shortened multiplet obeys additional annihilation equations and can have a protected scaling dimension.

In a long multiplet, supersymmetry permits the dimension to vary with the coupling.

Let $\mathcal O$ be the Konishi scalar.

Its anomalous dimension $\gamma_K$ is specified by the response of its separated two-point correlation to a positive dilation $t$:

$$
\langle\mathcal O(tx)\mathcal O(0)\rangle=t^{-2(2+\gamma_K)}\langle\mathcal O(x)\mathcal O(0)\rangle,\qquad t>0,\quad x\ne0.
$$

The anomalous dimension $\gamma_K$ is the shift from the classical value two.

An instanton correction to this number changes the scaling of the correlation function.

Use the coupling $g$ and vacuum angle $\theta$ of Alday–Korchemsky:

$$
\tau=\frac\theta{2\pi}+\frac{4\pi i}{g^2},
\qquad e^{2\pi i\tau}=e^{-8\pi^2/g^2+i\theta}.
$$

The leading one-instanton and anti-instanton term of the source Konishi anomalous dimension is

$$
\left.\gamma_K\right|_{\mathrm{leading\ one\ instanton\ and\ anti\ instanton}}
=-\frac{27\kappa_N}{5(N^2-1)}
\left(\frac{g^2}{4\pi^2}\right)^2
(e^{2\pi i\tau}+e^{-2\pi i\bar\tau}),
$$

For the positive arguments used here, the Euler gamma function is defined by

$$
\Gamma(z)=\int_0^\infty t^{z-1}e^{-t}\,dt,\qquad z>0.
$$

$$
\kappa_N=\frac{2\Gamma(N-\tfrac12)}{\sqrt\pi\Gamma(N-1)}.
$$

Here $\Gamma$ is the Euler gamma function.

Direct evaluation gives $\kappa_2=1$ and the $SU(2)$ coefficient $-9/5$.

The source Konishi profile at fermion degree eight satisfies

$$
K^{(4)}=0,
\qquad
K^{(8)}(x)=-9\cdot2^{15}
\frac{\rho^6[\zeta(x)]^8}{[\rho^2+(x-x_0)^2]^6}.
$$

Here $x_0,\rho$ are the instanton center and size, and $[\zeta(x)]^8$ is the source ordered product of eight exact modes.

Two $K^{(8)}$ insertions saturate the sixteen exact modes.

The complete zero-mode calculation uses these higher-degree profiles.

Konishi belongs to a long multiplet.

[Alday–Korchemsky](https://arxiv.org/abs/1605.06346)

<a id="np-9487cf1fcea2eac8-1"></a>

### Positivity constraint on leading Fourier sectors

The supersymmetry Hamiltonian measures the failure of a state to remain annihilated by the chosen supercharge and its adjoint.

Its expectation value is nonnegative.

A Fourier sector is a coefficient of a phase $e^{ik\theta}$ in the dependence on the vacuum angle; the neutral sector has $k=0$.

The following assumption excludes a neutral term at the same leading exponential order.

Take a finite-dimensional space of perturbative BPS survivors at fixed charge, with positive radial inner product.

Let $H$ be its nonnegative supersymmetry Hamiltonian, and $A$ the proposed leading instanton Fourier matrix.

Assume the leading exponential order has vanishing neutral term:

$$
H_{\mathrm{lead}}(\theta)
=e^{-8\pi^2/g^2}(e^{i\theta}A+e^{-i\theta}A^\dagger).
$$

> **Proof.**
>
> Exact positivity requires a positive semidefinite leading coefficient for each fixed $\theta$.
>
> Its angular average vanishes:
>
> $$
> \int_0^{2\pi}H_{\mathrm{lead}}(\theta)\,d\theta=0,
> \qquad \operatorname{Tr}H_{\mathrm{lead}}(\theta)\ge0.
> $$
>
> The continuous nonnegative trace has zero integral and therefore vanishes identically.
>
> The positive semidefinite matrix vanishes identically, giving $A=0$.
>
> $\square$

If the transferred effective supercharge starts with $e^{2\pi i\tau}d$, where $d$ is the instanton coefficient, its leading Hamiltonian is

$$
\left.\{Q_{\mathrm{eff}},Q_{\mathrm{eff}}^\dagger\}\right|_{\mathrm{lead}}
=e^{-16\pi^2/g^2}\{d,d^\dagger\}\ge0.
$$

This neutral quadratic contribution satisfies positivity.

A vanishing linear instanton energy shift is therefore compatible with a nonzero instanton supercharge differential.

Finding a $\mathcal N=4$ survivor certified through the complete perturbative complex, together with a nonzero microscopic instanton image, remains open beyond the pure-scalar sector.
