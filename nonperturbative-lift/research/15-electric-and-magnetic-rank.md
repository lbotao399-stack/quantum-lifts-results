<a id="np-section-13"></a>

## 15. SQCD: electric rank identities and magnetic nonperturbative dynamics

Supersymmetric QCD (SQCD) couples gauge fields to matter multiplets, each with a complex scalar and a fermion.

For $SU(N_c)$, a fundamental scalar is a color column and an antifundamental scalar a dual row; $N_f$ counts pairs of these multiplets.

In the meson formula below, $Q,\widetilde Q$ denote these quark scalars.

Electric and magnetic label two descriptions related by Seiberg duality: their fields and gauge groups differ, while their low-energy gauge-invariant physics agrees.

A meson contracts one quark with one antiquark.

A baryon contracts $N_c$ quarks, or $N_c$ antiquarks, to a gauge-invariant scalar.

The superpotential $W$ is a holomorphic function of chiral matter fields, whose lowest scalar components are annihilated by the antichiral supercharges: its dependence is on complex scalar coordinates, independently of their complex conjugates.

The F-equations set its first derivatives to zero.

The gauge D-equations set the gauge auxiliary fields to zero.

Together with these equations and identification under gauge transformations, the F-equations determine supersymmetric vacua.

In the effective scalar operator ring, these first derivatives are set equal to zero as operator relations.

The scales $\Lambda,\widetilde\Lambda$ belong to the electric and magnetic descriptions; $\mu$ fixes their meson normalization.

The generic magnetic description in this section has $N_c+2\le N_f<3N_c$; its mass deformations use the flavor counts specified below.

Here $N_f$ counts the fundamental–antifundamental flavor pairs of $SU(N_c)$.

The electric meson is $M^i{}_{\tilde\jmath}=Q^{ia}\widetilde Q_{a\tilde\jmath}$; $a$ is an abstract color slot, and $i,\tilde\jmath$ are component labels for the two flavor spaces.

The electric polynomial identity therefore gives $\operatorname{rank}M\le N_c$.

In the magnetic theory, $M$ is an independent singlet, and the number of colors is $\widetilde N_c=N_f-N_c$.

A gauge singlet is unchanged by the magnetic gauge group.

Thus the magnetic matrix $M$ is an independent field before its F-equations are imposed.

The electric expression $Q\widetilde Q$ already has rank at most $N_c$ because its color contraction passes through an $N_c$-dimensional space.

The magnetic interactions must recover this restriction dynamically.

We use the dimension-two meson and matching scale $\mu$ of Intriligator–Seiberg:

$$
W_{\mathrm{tree}}=\frac1\mu M^i{}_{\tilde\jmath}q_i{}^A\widetilde q_A{}^{\tilde\jmath},
\qquad
\Lambda^{3N_c-N_f}\widetilde\Lambda^{3\widetilde N_c-N_f}
=(-1)^{\widetilde N_c}\mu^{N_f}.
$$

$A$ is an abstract magnetic color slot.

At vanishing holomorphic beta-function exponent, the source notation $\widetilde\Lambda^0$ in scale matching denotes a dimensionless instanton factor, distinct from the number $1$.

The zero-exponent endpoint of $SP$ below uses the same convention.

At $q=\widetilde q=0$, every $M$ satisfies the tree F-equations.

Perturbative superpotential nonrenormalization preserves this magnetic branch.

Magnetic strong dynamics restores the electric rank bound.

[Intriligator–Seiberg, §§4–5](https://arxiv.org/abs/hep-th/9509066)

<a id="np-a2081c7e96050b3d-1"></a>

### Rank patches and exact massive vacua

On the patch $\operatorname{rank}M=r$, $\det{}'(M/\mu)$ denotes the determinant of the invertible heavy flavor block.

A rank patch is the region where a specified rank-$r$ meson block is invertible.

Its entries divided by $\mu$ give masses to $r$ magnetic quark pairs.

Eliminating these massive pairs changes the strong scale through the displayed threshold relation.

The remaining stationary equations decide whether that region contains a supersymmetric vacuum.

Holomorphic threshold matching gives

$$
\widetilde\Lambda_L^{3\widetilde N_c-(N_f-r)}
=\widetilde\Lambda^{3\widetilde N_c-N_f}\det{}'(M/\mu).
$$

For $N_c<r<N_f$, the remaining light magnetic meson is $N_i{}^{\tilde\jmath}=q_i{}^A\widetilde q_A{}^{\tilde\jmath}$, with flavor labels in the light block.

The ADS superpotential is

$$
W_{\mathrm{eff}}
=\frac1\mu M^i{}_{\tilde\jmath}N_i{}^{\tilde\jmath}
+(r-N_c)
\left[
\frac{\widetilde\Lambda^{3\widetilde N_c-N_f}\det{}'(M/\mu)}{\det N}
\right]^{1/(r-N_c)}.
$$

> **Proof: rank exclusion.**
>
> The singlet F-equations require $N=0$.
>
> The integrated-in magnetic glueball $\widetilde S$ of the same source has stationary equation
>
> $$
> \widetilde S^{r-N_c}\det N
> =\widetilde\Lambda^{3\widetilde N_c-N_f}\det{}'(M/\mu)\ne0.
> $$
>
> The two equations are incompatible, so this rank patch contains no finite stationary vacuum.
>
> On the full-rank patch, all magnetic quarks are massive:
>
> $$
> \widetilde S^{\widetilde N_c}
> =(-1)^{\widetilde N_c}\frac{\det M}{\Lambda^{3N_c-N_f}},
> \qquad
> W_{\mathrm{eff}}=\widetilde N_c\widetilde S,
> $$
>
> $$
> \frac{\partial W_{\mathrm{eff}}}{\partial M^i{}_{\tilde\jmath}}
> =\widetilde S(M^{-1})^{\tilde\jmath}{}_i\ne0.
> $$
>
> The exact magnetic vacuum space therefore satisfies $\operatorname{rank}\langle M\rangle\le N_c$.
>
> At $r=N_c$, the remaining theory has equal numbers of flavors and colors:
>
> $$
> \det N-b\widetilde b=\widetilde\Lambda_L^{2\widetilde N_c},
> \qquad N=0,
> \qquad b\widetilde b=-\widetilde\Lambda_L^{2\widetilde N_c}.
> $$
>
> The source baryon map gives $\langle B\widetilde B\rangle=\det{}'\langle M\rangle$.
>
> $\square$

An invertible electric mass matrix $m^{\tilde\jmath}{}_i$ adds the superpotential $\operatorname{Tr}(mM)$.

The massive singlet equations now determine the vacua; the preceding rank bound pertains to the massless vacuum space.

> **Proof: $N_c$ massive vacua.**
>
> $$
> \widetilde S M^{-1}+m=0,
> \qquad M=-\widetilde S m^{-1},
> \qquad
> \det M=\frac{(-\widetilde S)^{N_f}}{\det m}.
> $$
>
> Substitution into the condensate equation gives
>
> $$
> (-\widetilde S)^{N_c}=\Lambda^{3N_c-N_f}\det m.
> $$
>
> The source electric glueball is $S=-\widetilde S$; $k=1,\ldots,N_c$ labels the vacua.
>
> $$
> S_k=e^{2\pi i k/N_c}
> (\Lambda^{3N_c-N_f}\det m)^{1/N_c},
> \qquad \langle M\rangle_k=S_km^{-1},
> $$
>
> $$
> W_{\mathrm{vac},k}
> =\widetilde N_c\widetilde S_k+\operatorname{Tr}(mM_k)
> =N_cS_k.
> $$
>
> For $\det m\ne0$, the roots are distinct, yielding exactly $N_c$ vacua.
>
> $\square$

> **Example 15.1. **
>
> **Rank exclusion and massive vacua.**
>
> $SU(3)$ with $N_f=5$ gives a concrete conformal-window example:
>
> $$
> \widetilde N_c=2,
> \qquad \Lambda^4\widetilde\Lambda=\mu^5.
> $$
>
> A rank-four meson leaves one-flavor magnetic $SU(2)$, whose ADS term excludes that patch.
>
> The two sheets of the full-rank magnetic condensate, together with the dynamical singlet equations, give
>
> $$
> \widetilde\Lambda_L^5=\frac{\mu\det{}'M}{\Lambda^4},
> \qquad
> W_{\mathrm{rank}\,4}
> =\frac1\mu M^5{}_{\tilde5}N_5{}^{\tilde5}
> +\frac{\mu\det{}'M}{\Lambda^4N_5{}^{\tilde5}},
> $$
>
> $$
> W_{\mathrm{vac},k}
> =3e^{2\pi i k/3}(\Lambda^4\det m)^{1/3},
> \qquad k=1,2,3.
> $$

These calculations establish vacuum-rank enforcement and the massive vacua.

The full nonreduced local-operator ring retains nilpotents; its identification requires matching operator products.
