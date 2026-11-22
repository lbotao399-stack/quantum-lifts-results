# Independent audit: the finite-circle SQCD monopole lift

**Verdict: yes, in the specified effective local chiral ring; not yet a four-dimensional microscopic polynomial-operator Q-lift.**

Source: [Aharony–Razamat–Seiberg–Willett, arXiv:1305.3924](https://arxiv.org/abs/1305.3924), §3.1, printed p.19, (3.3), and p.20, (3.5)–(3.6). Printed pages are PDF pages minus one. For $SU(N_c)$ with $N_f=N_c\ge2$, the source uses the meson matrix $M$, baryons $B,\widetilde B$, and monopole/Coulomb-branch coordinate $Y$:

$$
W=Y\left(B\widetilde B-\det M+\eta\right).
$$

The coefficient and sign of $\eta Y$ are exactly those of (3.5). The paper identifies its origin with the additional circle instanton-monopole in §2.2, printed p.15, (2.10), and identifies $\eta=\Lambda^{2N_c}$ in the normalization of (3.6). No separate numerical normalization of the ultraviolet disorder operator is inferred here.

## Exact ideal membership

With $M^i{}_j$ the meson entries, $i,j$ label the two flavor component spaces; differentiation with respect to an entry has the dual index type. Euler's homogeneous-polynomial identity gives

$$
M^i{}_j\frac{\partial\det M}{\partial M^i{}_j}=N_c\det M.
$$

> **Proof.** Direct substitution gives
>
> $$
> Y\frac{\partial W}{\partial Y}
> -B\frac{\partial W}{\partial B}
> -\frac1{N_c}M^i{}_j\frac{\partial W}{\partial M^i{}_j}
> =Y(B\widetilde B-\det M+\eta)
> -YB\widetilde B+Y\det M
> =\eta Y.
> $$
>
> For fixed $\eta\ne0$, this is a global polynomial F-term-ideal certificate for $Y=0$. It requires no inverse meson, determinant or baryon.

At $\eta=0$, set $M=B=\widetilde B=0$. Because $N_c\ge2$, every first derivative of $\det M$ vanishes there. All F-equations vanish for arbitrary $Y$:

$$
\left.\frac{\partial W}{\partial Y}\right|_{M=B=\widetilde B=\eta=0}=0,
\quad
\left.\frac{\partial W}{\partial B}\right|_{M=B=\widetilde B=\eta=0}=0,
\quad
\left.\frac{\partial W}{\partial\widetilde B}\right|_{M=B=\widetilde B=\eta=0}=0,
\quad
\left.\frac{\partial W}{\partial M^i{}_j}\right|_{M=B=\widetilde B=\eta=0}=0.
$$

Evaluation on this line maps the chiral ring onto $\mathbb C[Y]$. Thus $[Y]\ne0$, and in fact every positive power of $[Y]$ is nonzero, before the $\eta Y$ deformation. The restriction $N_c\ge2$ is necessary: a one-by-one determinant has nonzero constant derivative.

## Odd primitive and the parameter caveat

Use the manuscript's retained scalar/odd-letter convention, namely that the $\gamma$ letters are antichiral scalar coordinates and the dual odd $\beta$ letters obey $Q\beta=-\partial\bar W/\partial\gamma$ in the effective F-term Koszul model. Here the species labels on $\beta^Y,\beta^B,\beta^M$ refer to the coordinates $\bar Y,\bar B,\bar M$; barred flavor indices belong to the conjugate flavor spaces. Then

$$
Q\left[
\frac{-\bar Y\beta^Y+\bar B\beta^B
+\dfrac1{N_c}\bar M^{\bar i}{}_{\bar j}
\beta^M{}_{\bar i}{}^{\bar j}}
{\bar\eta}
\right]=\bar Y.
$$

The primitive is polynomial in dynamical fields and regular throughout the effective field space at fixed nonzero $\eta$. It is not regular at $\eta=0$: its coefficient contains $1/\bar\eta$. Before inverting the deformation parameter, the direct statement is

$$
\bar\eta[\bar Y]=0.
$$

The numerator is an odd closed cochain at $\eta=0$; its image under the term linear in $\bar\eta$ is $\bar Y$. This is the corresponding first deformation differential killing the scalar class. The nonzero special-fiber class and the exact nonzero-parameter class are compatible with parameter torsion; they are not a contradiction.

## Scope of the positive result

- This is a genuine removal of a nonzero operator class, not merely exclusion of a vacuum or deformation of a relation that was already zero.
- The undeformed input is the exact three-dimensional $\eta=0$ chiral-ring description (3.3). Other three-dimensional monopole effects have already been included. It is not a bare four-dimensional tree-level polynomial ring.
- $Y$ is a local monopole disorder operator in three dimensions, defined through a magnetic singularity. It is not an ordinary polynomial ultraviolet local operator on uncompactified four-dimensional spacetime. The source's printed p.27, footnote 6, also requires a UV completion for a standalone three-dimensional Lagrangian containing a monopole superpotential; the compactified four-dimensional origin supplies the relevant setting here.
- The explicit odd primitive is proved in the effective F-term Koszul model. Transport to the manuscript's one-supercharge microscopic holomorphic Ward complex requires its operator map and a compatible treatment of disorder sectors. Equation (3.5) alone does not supply that chain map.
- The assertion is global within the source's effective presentation and fixed nonzero-parameter fiber. It does not assert a primitive regular in the undeformed formal parameter ring.

Thus this is a valid finite-circle, monopole-extended effective chiral-ring benchmark for nonperturbative lifting, with an explicit ideal certificate and effective odd primitive. It should not be relabeled as a completed lift of an uncompactified four-dimensional polynomial-operator class.
