# Beasley–Witten higher F-terms: source audit

Source: C. Beasley and E. Witten, *New Instanton Effects In Supersymmetric QCD*, [hep-th/0409149v1](https://arxiv.org/abs/hep-th/0409149), §§2.2–2.3, 3.2, 4.1–4.3. Page numbers below are printed page numbers; the PDF page number is one larger. Relevant formulas and bar placements were checked against rendered PDF pages.

## 1. SU(2), four flavors: what the instanton generates

The source writes the number of flavors as $n=N_f$; SU(2) then has $2n$ quark doublets. $N_c$ is the number of colors. In a one-instanton background, the gaugino and quark zero-mode counts are respectively $2N_c$ and $2N_f$. Yukawa interactions pair and lift $2(N_c-1)$ modes of each species (§4.1, pp.23–27):

$$
2N_c-2(N_c-1)=2,
\qquad
2N_f-2(N_c-1)=2(N_f-N_c+1).
$$

The source's $p$ counts fermion pairs in the chiral integrand. The two exact supersymmetry modes provide the chiral superspace measure $d^2\theta$:

$$
p=N_f-N_c+1,
\qquad
2p+2=2(N_f-N_c)+4.
$$

For SU(2), $N_f=4$:

$$
\underbrace{4}_{\text{gaugino}}+
\underbrace{8}_{\text{quark}}-
\underbrace{2}_{\text{lifted gaugino}}-
\underbrace{2}_{\text{lifted quark}}
=2+6=8,
\qquad p=3.
$$

Thus the chiral integrand has six fermions; the component interaction probed in Eqs.(4.5)–(4.6) has eight. With the source's holomorphic dynamical scale $\Lambda$, the one-instanton factor is

$$
\Lambda^{3N_c-N_f}=\Lambda^{6-4}=\Lambda^2.
$$

The two lifted gaugino modes are the superconformal modes. Two Yukawa insertions absorb them together with two quark modes. The six remaining quark modes are absorbed by six massless nonsinglet external fermions. The two exact supersymmetry modes acquire singlet-quark components through the coupled equations (4.11)–(4.13); two singlet external fermions absorb these modes. The correlator uses fields conjugate to those in the effective vertex because it includes external propagators (p.25, footnote 4).

The nonzero Higgs expectation value $v$ regulates the instanton size $\rho$. For the radial integral in Eq.(4.14), with $g$ the source's gauge coupling, the $n=4$ value is exactly

$$
\int_0^\infty d\rho\,\rho^{13}
\exp\!\left[-\frac{4\pi^2|v|^2}{g^2}\rho^2\right]
=\frac{6!}{2}
\left(\frac{g^2}{4\pi^2|v|^2}\right)^7
=360\left(\frac{g^2}{4\pi^2|v|^2}\right)^7.
$$

The SU(2) orientation integral is nonzero by the positivity argument after Eq.(4.15); its tensor structure follows from group invariance. The radial integral above does not fix the absolute interaction coefficient: the source leaves the collective-coordinate measure and zero-mode normalization unfixed.

## 2. The cohomology actually established

The source's $\mathcal M$ is the vacuum moduli space; $\Phi^i$ are its local chiral coordinates, and barred indices label antiholomorphic coordinate components. $\bar D_{\dot\alpha}$ is the antichiral superspace derivative. The coefficient tensor $\omega$ is antisymmetric within each of its two $p$-index sets and symmetric under exchange of the sets (Eq.(2.28)):

$$
\delta S=\int d^4x\,d^2\theta\,\mathcal O_\omega,
$$

$$
\mathcal O_\omega=
\omega_{\bar i_1\cdots\bar i_p\,\bar j_1\cdots\bar j_p}
\left(\bar D_{\dot\alpha_1}\bar\Phi^{\bar i_1}
      \bar D^{\dot\alpha_1}\bar\Phi^{\bar j_1}\right)
\cdots
\left(\bar D_{\dot\alpha_p}\bar\Phi^{\bar i_p}
      \bar D^{\dot\alpha_p}\bar\Phi^{\bar j_p}\right).
$$

Spinor indices are contracted with the source's spinor metric. Using the supplied Kähler metric to raise one antiholomorphic index set turns that set into holomorphic tangent slots. The source denotes antiholomorphic $p$-forms by $\overline\Omega^p_{\mathcal M}$, and $T\mathcal M$ denotes the holomorphic tangent bundle. Thus $\omega$ becomes a section of

$$
\overline\Omega^p_{\mathcal M}\otimes\wedge^p T\mathcal M,
\qquad \bar\partial\omega=0.
$$

Here $\bar\partial$ is the antiholomorphic differential. This is the **linearized, on-shell supersymmetry condition in the unperturbed IR sigma model** (§2.2, pp.12–13).

The quotient is not literally ordinary Dolbeault cohomology. Equation (2.34) symmetrizes the primitive contribution under exchange of the two index sets. A primitive arising from a three-quarter-superspace integral must also obey the additional supersymmetry condition stated on p.14. Therefore replacing the physical classification without qualification by $H^p(\mathcal M,\wedge^pT\mathcal M)$ loses a source condition.

For the meson matrix $M^{ij}$, whose flavor labels run from $1$ to $2n$, the classical constraint is

$$
M\wedge M=0,
\qquad \operatorname{rank}M\le2.
$$

Here the wedge is the source's antisymmetrized product of flavor tensors. Excising the singular origin gives the source's regular space $\widetilde{\mathcal M}=\mathcal M-\{0\}$, a $\mathbb C^*$ bundle over $B=\operatorname{Gr}(2,2n)$, the space of complex two-planes in $\mathbb C^{2n}$. Its homogeneous line bundle $\mathcal L$ records meson scaling. The higher F-term pulls back from a section of (§§3.1–3.2)

$$
\overline\Omega^{n-1}_B\otimes\wedge^{n-1}TB\otimes\mathcal L^{-n}.
$$

For $n=4$, this is

$$
B=\operatorname{Gr}(2,8),
\qquad
\overline\Omega^3_B\otimes\wedge^3TB\otimes\mathcal L^{-4}.
$$

The chirality and primitive tests on p.21 reduce to the unbroken SU(2) representation. The candidate $\bar\partial\omega$ contains $2p+1=7$ doublets; a candidate primitive contains $2p-1=5$ doublets. The central group element $-1$ acts as

$$
(-1)^7=-1,
\qquad (-1)^5=-1.
$$

Neither tensor product contains an invariant vector. Hence the invariant tensor is closed and has no flavor-invariant admissible primitive. For globally regular primitives, averaging over the compact flavor group would turn any primitive into an invariant one; the admissibility conditions are linear and flavor invariant. This extends the source's symmetry argument to exclusion of a global admissible primitive.

On a sufficiently small neighborhood of one Higgs vacuum the positive-degree class is trivial and the interaction is a D-term (§4.2, pp.30–31). The local primitives fail to define one global primitive. **This is nontrivial interaction data on the regular IR Higgs branch. It does not establish removal of an ordinary scalar chiral-ring element or a UV holomorphic-jet/BV class.** For $N_f>N_c$, the source explicitly says the higher F-terms generate no superpotential and do not deform the complex structure of the vacuum moduli space (Introduction).

## 3. What the Seiberg-dual calculation actually computes

Section 4.2 concerns **SU(2), $N_f=n=3$, six doublets**. Its dual gauge group is trivial. Its fields are the mesons, with the Wess–Zumino action (4.17)

$$
S=\frac1{\mu^2}\int d^4x\,d^4\theta\,\overline M M
+\int d^4x\,d^2\theta\,\Lambda^{-3}M\wedge M\wedge M
+\mathrm{c.c.}
$$

The scale $\mu$ allows the meson to retain engineering dimension two. The supplied Hermitian pairing is $\overline M M=\tfrac12\overline M_{ij}M^{ij}$; the conjugate meson belongs to the conjugate representation, with the pairing explicitly supplied. The wedge convention is the source's and has not been converted to another paper's Pfaffian convention. Its F equations give

$$
M\wedge M=0.
$$

Expand about $\langle M^{12}\rangle\ne0$. Integrating out massive meson components produces the six-point tree supergraph, Fig.4 and Eq.(4.28). The apparent zero-momentum pole cancels using $\bar D^2D^2=p^2$ on the relevant chiral field **in the source's superspace conventions**; here $p$ is momentum, unrelated to the earlier form degree.

The source defines $\lambda=\Lambda^{-3}$, Eq.(4.23). Extracting only the coupling factors from Eq.(4.28) gives the exact algebraic cancellation

$$
\frac{\bar\lambda}{\lambda\bar\lambda}\,
\lambda\,
\frac{\bar\lambda}{\lambda\bar\lambda}\,
\lambda\,
\frac{\bar\lambda}{\lambda\bar\lambda}
=\lambda^{-1}=\Lambda^3.
$$

The local D-term obtained in this way is the local form of the global higher F-term (4.18). This does not conflict with ordinary superpotential nonrenormalization: the generated local interaction is a D-term, and its globally nontrivial higher F-term class is not an ordinary superpotential.

For this theory,

$$
p=3-2+1=2,
\qquad 2p+2=6,
\qquad \Lambda^{6-3}=\Lambda^3.
$$

**This is not the magnetic dual calculation for $N_f=4$.** SU(2), $N_f=4$ is self-dual with magnetic SU(2), $N_f=4$. The relation $N_f=4\to3$ in the present paper is mass deformation, not that duality map.

## 4. Mass deformation and descent

With a holomorphic superpotential $W$, the source replaces the chirality differential by (2.37)

$$
\delta=\bar\partial+\iota_{dW},
\qquad
\delta^2=\bar\partial^2+
\{\bar\partial,\iota_{dW}\}+\iota_{dW}^{\,2}=0.
$$

Here $\iota_{dW}$ contracts a holomorphic tangent slot with the one-form $dW$; its anticommutator with $\bar\partial$ vanishes because $W$ is holomorphic. A mass deformation mixes form/polyvector degrees; Eq.(4.37) constructs an inhomogeneous $\delta$-closed completion. The source uses $W=m_{ij}M^{ij}$, with antisymmetric mass matrix $m_{ij}$, without inserting a new factor $1/2$.

The degree-$k$ term in the masses has the binomial coefficient (4.43)

$$
(-1)^k\binom{n-1}{k}.
$$

Integrating out the massive flavors gives the lower-flavor higher F-term. In particular,

$$
n:4\longrightarrow3\longrightarrow2\longrightarrow1,
\qquad
p:3\longrightarrow2\longrightarrow1\longrightarrow0,
$$

$$
8\text{-fermion}\longrightarrow6\text{-fermion}
\longrightarrow4\text{-fermion complex-structure deformation}
\longrightarrow\text{ADS superpotential}.
$$

At the final step the source writes $M=M^{12}$ for the sole light meson and obtains (4.32)

$$
W=\frac{\Lambda^5}{M}.
$$

Equations (4.44)–(4.45) include the heavy-mass Pfaffian, the source sign, and an additional factor two from contraction of the remaining meson indices. Ultraviolet/infrared scale matching then gives (4.32). These expressions should not be transferred between different mass, wedge/Pfaffian, or scale conventions without matching them.

## 5. Coefficients and the boundary of the result

- **p.26, after (4.7):** “we will not try to compute the absolute normalization of our interaction”; the zero modes are not normalized.
- **p.28, after (4.14):** powers of $g^2$ in the measure and Gaussian integration were not tracked. Explicit gauge-coupling dependence is assigned to external wavefunction normalization by holomorphy.
- **p.32, after (4.23):** numerical constants are ignored; the factor three in the expanded cubic superpotential is absorbed into the vacuum expectation value; $\mu$ is suppressed.

The source establishes a nonzero one-instanton higher F-term, its fermion count, holomorphic scale power, tensor structure, global cohomological nontriviality, and mass-deformation consistency. It does not supply an absolute coefficient in a separately normalized UV/BV composite-operator complex. Nor does it supply a chain map from that complex to the IR sigma-model complex, a nonperturbative differential on the UV complex, or a demonstration that a specified perturbative class becomes exact.

No manuscript files were changed by this audit.
