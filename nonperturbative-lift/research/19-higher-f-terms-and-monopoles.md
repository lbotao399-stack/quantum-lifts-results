<a id="np-section-17"></a>

## 19. Higher F-terms and monopole-extended rings

<a id="np-b5341a21dac18f63-1"></a>

### Instanton data beyond the ordinary superpotential

An ordinary superpotential contributes an integral over two anticommuting supersymmetry coordinates.

A higher F-term uses the same two-coordinate integral but inserts extra fermion factors.

A D-term integrates over all four supersymmetry coordinates.

These names specify the integration and field content of the interaction.

An instanton fermion zero mode is an anticommuting parameter left free by the linearized fermion equation in that instanton background.

Its integral vanishes unless an insertion supplies that parameter.

The two universal supersymmetry modes supply the chiral integration; the additional matter modes require the extra fermion factors below.

In the Beasley–Witten example with $SU(2)$ and $N_f=4$, one instanton retains two supersymmetry modes and six matter modes.

Here $p=N_f-N_c+1=3$ counts the fermion pairs in the chiral integrand:

$$
2N_c-2(N_c-1)=2,
\qquad
2N_f-2(N_c-1)=6,
\qquad 2p+2=8.
$$

This yields an eight-fermion component interaction, with instanton scale factor $\Lambda^2$.

For source Higgs expectation $v$ and instanton size $\rho$, the radial integral is

$$
\int_0^\infty d\rho\,\rho^{13}
\exp\left[-\frac{4\pi^2|v|^2}{g^2}\rho^2\right]
=360\left(\frac{g^2}{4\pi^2|v|^2}\right)^7.
$$

The absolute interaction coefficient remains unspecified in the source; determining it requires the zero-mode normalization and collective-coordinate measure together with this integral.

Let $\Phi^i$ be local chiral coordinates on the vacuum moduli space $\mathcal M$, and let $\omega$ be the source coefficient tensor:

The vacuum moduli space $\mathcal M$ is the continuous family of supersymmetric vacua after identifying gauge-equivalent configurations.

Each chiral field $\Phi^i$ packages one complex coordinate with its supersymmetry partners.

The tensor $\omega$ is the field-dependent interaction coefficient.

The operator $\bar D_{\dot\alpha}$ is the antichiral superspace derivative; its action on $\bar\Phi$ supplies a fermionic factor.

The barred flavor slots in the following expression contract these factors with $\omega$.

$$
\delta S=\int d^4x\,d^2\theta\,\mathcal O_\omega,
$$

$$
\mathcal O_\omega=
\omega_{\bar i_1\cdots\bar i_p\,\bar j_1\cdots\bar j_p}
\prod_{\substack{r=1,\ldots,p\\\mathrm{fermion\ pairs}}}
\left(\bar D_{\dot\alpha_r}\bar\Phi^{\bar i_r}
\bar D^{\dot\alpha_r}\bar\Phi^{\bar j_r}\right).
$$

The source Kähler metric converts one group of antiholomorphic covariant slots into holomorphic tangent slots.

Linearized chirality requires $\bar\partial\omega=0$.

On the regular Higgs branch of $SU(2)$ with $N_f=n$, the class descends from the Grassmannian $B=\operatorname{Gr}(2,2n)$ and has tensor type

$$
\overline\Omega_B^{n-1}\otimes\wedge^{n-1}TB\otimes\mathcal L^{-n}.
$$

Here $\overline\Omega_B^{n-1}$ denotes antiholomorphic forms, $TB$ is the holomorphic tangent bundle, and $\mathcal L$ records meson scaling.

Admissible primitives in the physical quotient also obey the source three-quarter-superspace condition and index-set symmetrization.

A primitive here is an expression whose supersymmetry derivative reproduces the interaction.

Admissibility requires the source index symmetries and regularity throughout the vacuum region under consideration.

Three-quarter superspace means integration over three of the four anticommuting supersymmetry coordinates.

A globally regular primitive of this allowed form would remove the higher F-term interaction class.

Local primitives on separate coordinate patches can fail to agree on their overlaps; this is the global obstruction tested below.

> **Proof: source symmetry test for global nontriviality.**
>
> For $n=4$, the unbroken $SU(2)$ tensor products of the candidate derivative and primitive contain $7$ and $5$ doublets, respectively.
>
> The central element $-1$ acts by $(-1)^7=(-1)^5=-1$, so both invariant-vector spaces vanish.
>
> Averaging over the compact flavor group takes any globally regular admissible primitive to an invariant primitive while preserving the linear admissibility conditions.
>
> Consequently, a global admissible primitive is excluded.
>
> Local D-term primitives exist on small Higgs patches; their global gluing gives the nontrivial higher F-term.
>
> $\square$

The dual Wess–Zumino theory at $N_f=3$ gives a tree calculation of this mechanism.

With the source $\lambda=\Lambda^{-3}$, the coupling factors of the six-point tree graph combine to

$$
\frac{\bar\lambda}{\lambda\bar\lambda}\lambda
\frac{\bar\lambda}{\lambda\bar\lambda}\lambda
\frac{\bar\lambda}{\lambda\bar\lambda}
=\lambda^{-1}=\Lambda^3.
$$

Under mass deformation, the source chirality differential becomes $\delta=\bar\partial+\iota_{dW}$.

Here $\iota_{dW}$ contracts a holomorphic tangent slot with $dW$, and holomorphy gives $\delta^2=0$.

The degree-$k$ mass term has coefficient $(-1)^k\binom{n-1}{k}$, with descent

$$
N_f:4\longrightarrow3\longrightarrow2\longrightarrow1,
\qquad
p:3\longrightarrow2\longrightarrow1\longrightarrow0.
$$

The resulting interactions are, in order, the eight-fermion term, six-fermion term, complex-structure deformation, and ADS superpotential.

Thus a systematic account of nonperturbative effects includes higher F-terms, with the ordinary scalar superpotential $W$ as one component.

[Beasley–Witten](https://arxiv.org/abs/hep-th/0409149)

<a id="np-4e3f9fd0ffa2c49b-1"></a>

### Explicit monopole lift on a finite circle

Compactification makes one spacetime direction periodic.

At distances larger than the circle, the fields form a three-dimensional theory.

A local monopole operator is specified by a magnetic-flux boundary condition on a small sphere surrounding its insertion.

The operator $Y$ below records such an insertion.

The extra finite-circle monopole configuration generates the term $\eta Y$, with the source parameter $\eta=\Lambda^{2N_c}$.

For compactified $SU(N_c)$ with $N_f=N_c\ge2$, the source effective superpotential is

$$
W=Y(B\widetilde B-\det M+\eta),
\qquad \eta=\Lambda^{2N_c}.
$$

Here $Y$ is a three-dimensional local monopole disorder operator, and $\eta Y$ arises from the finite-circle KK monopole.

The starting presentation at $\eta=0$ already includes the remaining three-dimensional nonperturbative dynamics.

Monopole-class removal means that the local operator $Y$ becomes a supercharge variation in this effective description.

The calculation compares the ring at zero circle parameter with the ring at nonzero circle parameter.

It retains the pre-existing three-dimensional interactions in both cases.

> **Proof: class removal and the full scalar quotient.**
>
> The Euler identity gives
>
> $$
> Y\frac{\partial W}{\partial Y}
> -B\frac{\partial W}{\partial B}
> -\frac1{N_c}M^i{}_{\tilde\jmath}
> \frac{\partial W}{\partial M^i{}_{\tilde\jmath}}
> =\eta Y.
> $$
>
> At $\eta=0$, setting $M=B=\widetilde B=0$ leaves an F-flat line with arbitrary $Y$.
>
> Evaluation into $\mathbb C[Y]$ detects all positive powers.
>
> At fixed $\eta\ne0$, $Y$ belongs to the F-term ideal, and
>
> $$
> \frac{\mathbb C[M,B,\widetilde B,Y]}{(\partial W)}
> =\frac{\mathbb C[M,B,\widetilde B]}
> {(\det M-B\widetilde B-\eta)}.
> $$
>
> In the antichiral Koszul convention, the explicit primitive is
>
> <a id="np-66ee4aae6c4219a5-1"></a>
>
> $$
> Q\left[
> \frac{-\bar Y\beta^Y+\bar B\beta^B
> +\bar M^{\bar i}{}_{\bar{\tilde\jmath}}
> \beta^M{}_{\bar i}{}^{\bar{\tilde\jmath}}/N_c}{\bar\eta}
> \right]=\bar Y.
> $$
>
> Over the formal parameter ring, this gives $\bar\eta[\bar Y]=0$.
>
> Over $\mathbb C[[\bar\eta]]$, the allowed coefficients are formal series with nonnegative powers of $\bar\eta$.
>
> The equation $\bar\eta[\bar Y]=0$ is torsion: multiplication by the parameter kills the class.
>
> In the separately specified theory with fixed nonzero $\bar\eta$, its inverse is an allowed scalar coefficient, and the displayed primitive makes $\bar Y$ exact.
>
> $\square$

Comparison of the monopole disorder sector with the four-dimensional polynomial jet sector requires the operator definitions in both sectors.

[Aharony–Razamat–Seiberg–Willett, §3.1](https://arxiv.org/abs/1305.3924)

<a id="multiplier-pages"></a>

### Which page removes each scalar in a multiplier superpotential?

**Proposition 19.1.** Take a polynomial $F$ in scalar coordinates $x^i$ and the superpotential $W=Y(F+\eta)$ with independent parameter $\eta$. Work over the scalar polynomial ring, or its localization at a specified critical point, with a domain of scalar coefficients and $F\ne0$. In each positive power of $Y$, multiplication by $\eta$ in the scalar F-term quotient is precisely multiplication by $-F$ in the Jacobian algebra of $F$, the quotient by its first derivatives. The degree-zero sector has no $\eta$-torsion.

> **Proof.** The two sets of derivatives of $W$ are $F+\eta$ and $Y\,\partial F/\partial x^i$. Eliminating the parameter by its monic equation gives a direct decomposition by powers of $Y$,
>
> $$
> \frac{\mathbb C[x,Y,\eta]}{(F+\eta,\,Y\,\partial F)}
> =
> \mathbb C[x]\ \oplus\
> \bigoplus_{j\ge1}Y^j\frac{\mathbb C[x]}{(\partial F)},
> \qquad \eta=-F .
> $$
>
> Here $(\partial F)$ denotes the ideal of all first derivatives. The exact exponent that kills a specified class in positive $Y$ degree is the least $r$ for which $F^r$ times its coefficient belongs to $(\partial F)$. In degree zero multiplication by $-F$ is injective. If the local Jacobian algebra is finite dimensional and $F$ acts nilpotently there, decompose the action of $-F$ into nilpotent Jordan blocks. A block of length $r$ is a copy of $\mathbb C[\eta]/(\eta^r)$. It gives one independent special-fiber scalar target on page $r$. The top generator is killed by $\eta^r$. Adding $\eta$ times any class preserves its unit coefficient modulo $\eta$ in this block. A scalar has no outgoing Koszul differential. The equation $\eta^r$ times a scalar equals a boundary supplies the incoming primitive to order $r$. These statements persist under parameter-adic completion. $\square$

In the circle theory $F=B\widetilde B-\det M$. Its derivatives generate $B,\widetilde B$ and all signed minors of $M$ of size $N_c-1$. Euler's determinant identity puts $F$ in this ideal. Thus the entire parameter-torsion ideal is $(Y)$, killed on page one. In degree $j\ge1$ its coefficients form

$$
Y^j\,\mathbb C[M]\big/\bigl(\text{all }(N_c-1)\text{-minors of }M\bigr).
$$

The scalar torsion is supported at $B=\widetilde B=0$ and $\operatorname{rank}M\le N_c-2$ on the special fiber. The finite-circle superpotential is the physical input. [Aharony–Razamat–Seiberg–Willett](https://arxiv.org/abs/1305.3924)

> **Example 19.2.** A first-order change in a superpotential can give a second-page scalar target. Take the local polynomial model at $x=y=0$ with $F=x^5+y^5+x^3y^3$. In its local Jacobian algebra,
>
> $$
> 5x^4=-3x^2y^3,\qquad 5y^4=-3x^3y^2,\qquad
> F-\tfrac15(xF_x+yF_y)=-\tfrac15x^3y^3 .
> $$
>
> Every monomial of degree seven has an exponent at least four, so the first two relations raise its degree by one. Nakayama's lemma gives $(x,y)^7=0$ in the local quotient. In degrees at most six the leading homogeneous generators $5x^4,5y^4$ have their first polynomial relation in degree eight and are independent up to degree six. Reduction therefore gives the independent basis $x^i y^j$ with $0\le i,j<4$. In particular $x^3y^3\ne0$, while multiplication by either $x$ or $y$ kills it. Multiplication by $F$ has rank one and its square is zero. For every $j\ge1$ the $Y^j$ sector has one block of length two and fourteen blocks of length one. In particular
>
> $$
> \eta[Y^j]=\tfrac15[Y^j x^3y^3]\ne0,\qquad
> \eta^2[Y^j]=0 .
> $$
>
> There is an explicit local primitive. Subscripts $x,y$ denote derivatives, and the polynomial needed for division below has constant term $25$.
>
> $$
> (25-9xy)x^4y^4=5y^4F_x-3x^2y^3F_y .
> $$
>
> Multiplication by $x^2y^2$ and $5F+x^3y^3=xF_x+yF_y$ give
>
> $$
> 25F^2=(5F-x^3y^3)(xF_x+yF_y)+
> \frac{5x^2y^6 F_x-3x^4y^5 F_y}{25-9xy}.
> $$
>
> Write $x,y,Y,\eta$ for the local coefficients entering the Koszul sign convention $Q\beta^i=-\partial_i W$. Using $\eta^2=(\eta-F)(\eta+F)+F^2$ gives the full cochain equation
>
> $$
> Q\left[
> -Y(\eta-F)\beta^Y
> -\frac{(5F-x^3y^3)(x\beta^x+y\beta^y)}{25}
> -\frac{5x^2y^6\beta^x-3x^4y^5\beta^y}{25(25-9xy)}
> \right]=\eta^2Y .
> $$
>
> This is a local superpotential benchmark with a tunable coefficient $\eta$. The term linear in $\eta$ and the page-two boundary are both exact statements within this model. [Independent quotient calculation](../checks/scalar_saturation.py)

<a id="np-aa0e649d863ebf40-1"></a>

### Three-dimensional and worldsheet benchmarks

The source superpotential of Aharony duality is

$$
W=M_i{}^{\tilde i}q^i\widetilde q_{\tilde i}
+V_+\widetilde V_-+V_-\widetilde V_+.
$$

Here $V_\pm$ are elementary magnetic singlets corresponding to electric monopoles, and $\widetilde V_\pm$ are magnetic monopole operators.

Differentiating with respect to the elementary singlets gives

$$
[\widetilde V_-]=[\widetilde V_+]=[q^i\widetilde q_{\tilde i}]=0.
$$

Composite monopoles retain their definitions as disorder operators when deriving F-equations from elementary coordinates.

For the one-flavor $U(1)$ mirror, $W=-MV_+V_-$ retains three coordinate branches, and the three pairwise products vanish.

[Aharony](https://arxiv.org/abs/hep-th/9703215); [Aharony et al.](https://arxiv.org/abs/hep-th/9703110)

The minimal $(0,2)$ $\mathbb{CP}^1$ model has no left-moving fermions.

The two-dimensional target $\mathbb{CP}^1$ is the space of complex lines in $\mathbb C^2$.

The labels $(0,2)$ count left- and right-moving supersymmetries.

In the expression below, $\phi$ is a complex target coordinate, $R_{\phi\bar\phi}$ is the target Ricci tensor, and $\bar\psi_+$ is the right-moving fermion.

The odd operator $\theta$ is the specified curvature–fermion composite; $t_0$ sets the instanton weight $e^{-t_0}$.

Here $\mathcal O(d)$ denotes the degree-$d$ holomorphic line bundle on $\mathbb{CP}^1$, and $h^p$ is the dimension of its $p$-th sheaf cohomology.

Yagi uses $Q=\bar Q_+$ and $\theta=-R_{\phi\bar\phi}\partial_z\phi\,\bar\psi_+/2$.

The degree-one fermion bundle is

$$
\mathcal O(-1)\otimes\mathcal O(2)=\mathcal O(1),
\qquad h^0=2,\quad h^1=0.
$$

The supercurrent and $\theta$ saturate the two modes.

The source gives

$$
\{Q,\theta\}=\Lambda^2e^{-t_0}O_1,
\qquad O_1\text{ is a nonzero constant operator}.
$$

The numerical normalization of $O_1$ remains unspecified in the source.

Candidate operators at higher instanton degrees have negative dimensions and vanish.

At a fixed nonzero instanton factor, a closed $\mathcal O$ satisfies

$$
[Q,\theta\{Q,\theta\}^{-1}\mathcal O\}=\mathcal O.
$$

Formal coefficients with nonnegative powers instead give the torsion relation $e^{-t_0}[\mathcal O]=0$.

[Yagi](https://arxiv.org/abs/1001.0118); [Tan–Yagi](https://arxiv.org/abs/0801.4782)

The Hori–Vafa $(2,2)$ $\mathbb{CP}^{N-1}$ model gives a different result.

Here $\Sigma$ is the field-strength twisted chiral field, $Y_i$ are periodic dual matter fields, and $t$ is the FI–Theta parameter:

$$
\widetilde W
=\Sigma\left(\sum_{\substack{i=1,\ldots,N\\\mathrm{matter\ species}}}Y_i-t\right)
+\mu\sum_{\substack{i=1,\ldots,N\\\mathrm{matter\ species}}}e^{-Y_i}.
$$

The F-equations give

$$
\Sigma=\mu e^{-Y_i},
\qquad
\Sigma^N=\mu^Ne^{-t}=\Lambda^N,
\qquad \Lambda=\mu e^{-t/N}.
$$

The lowest component $\sigma$ has $N$ vacua $\sigma=\Lambda e^{2\pi i\ell/N}$, with $\ell=0,\ldots,N-1$.

This is a vortex-generated ring deformation.

[Hori–Vafa](https://arxiv.org/abs/hep-th/0002222)
