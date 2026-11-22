# Critical source audit: SU(5), Nf=4, k=2

**Verdict: OPEN for a microscopic nonperturbative lift.** Kutasov–Schwimmer–Seiberg (KSS) predict stronger nilpotence of the adjoint Casimir in their stated trace-sector presentation. They retain nilpotents, so their statement is stronger than equality of reduced vacuum sets. The cited passage does not independently establish the same square-zero identity in the full renormalized quantum chiral ring containing mesons, baryons and gauge-field-strength operators. The all-perturbative survival result for the specified square in the stated holomorphic Ward complex is unaffected by this source-scope refinement.

## 1. What the original source actually says

[KSS, hep-th/9510222](https://arxiv.org/abs/hep-th/9510222): printed page numbers are one less than PDF page numbers.

| Location | Source evidence | Consequence |
|---|---|---|
| p.4, footnote 3 | “Setting the quark fields to zero for the moment” | The opening adjoint-trace discussion has an explicit sector restriction. |
| pp.14–15, §2.3 | Lists generalized mesons, baryons and adjoint traces, then says of the traces “which we will focus on here” | The ensuing presentation is not an explicit complete presentation including all gauge-strength generators. |
| p.15, §2.3 | Adds magnetic classical characteristic relations to the electric ring because they eliminate unstable electric vacua | This supplies duality and vacuum evidence; it does not derive a microscopic composite-operator identity. |
| p.35, (6.10) | Gives the trace operator map | Its source normalization and positive sign must be retained. |
| p.36, (6.12)–(6.14) | Retains powers of the nilpotent trace; includes “setting $M_j=0$ as usual” | This is a nilpotent trace-sector prediction, not merely the reduced relation $\operatorname{Tr}X^2=0$. Its unrestricted full-ring promotion needs another argument. |
| p.37, case 1 | Discusses $N_c>N_f>\bar N_c$ and matching stable branches | The example $5>4>3$ lies in this case. |

Here $X,Y$ are the electric and magnetic adjoints, and $N_c,\bar N_c,N_f$ are their color ranks and common flavor number. For the source's cubic superpotential, $k=2$ and $\bar N_c=2N_f-N_c=3$. At zero quadratic deformation $s_1=0$, (6.10), (6.12), (6.14) give

$$
\operatorname{Tr}X^2=\frac53\operatorname{Tr}Y^2,
\qquad
(\operatorname{Tr}X^2)^3=0\quad\text{(classical electric trace sector)},
\qquad
(\operatorname{Tr}Y^2)^2=0\quad\text{(source magnetic trace presentation)}.
$$

The source prediction is consequently $(\operatorname{Tr}X^2)^2=0$ in that presentation. Its statement cannot be replaced by a vacuum-only assertion, but its full renormalized operator interpretation is not established by this citation alone.

KSS p.23, (3.39), explicitly allows meson-dependent trace-map terms for general $k$. Its index range $i=2,\ldots,k-1$ is empty at $k=2$: **(3.39) does not itself furnish a correction to (6.10).** A correction away from the mesonic origin must be derived in a definite common operator scheme.

## 2. A conditional repair of the dressed-glueball objection

KSS p.12, (2.14), fixes the magnetic superpotential, with the singlets $M_1,M_2$ representing the electric mesons $\widetilde QQ,\widetilde QXQ$ and with $q,\widetilde q$ the magnetic quarks:

$$
W_{\mathrm{mag}}
=-\frac{s_0}{3}\operatorname{Tr}Y^3
+\frac{s_0}{\mu^2}
\left(M_1\widetilde qYq+M_2\widetilde qq\right).
$$

Here $s_0$ and $\mu$ are the source coupling and matching scale. Flavor contractions in the superpotential are implicit, as in KSS. $\operatorname{Tr}$ is the color trace; $\operatorname{tr}$ below is the flavor trace. The following is conditional on the singlet F-equations and generalized Konishi identities being valid in one compatible composite-operator definition.

> **Proof.** The singlet equations give
>
> $$
> \widetilde qq=0,\qquad \widetilde qYq=0.
> $$
>
> Cyclic contraction then gives
>
> $$
> \operatorname{Tr}(qM_1\widetilde q)
> =\operatorname{tr}(M_1\widetilde qq)=0,
> $$
>
> $$
> \operatorname{Tr}\left[(qM_1\widetilde q)^2\right]
> =\operatorname{tr}\left(M_1\widetilde qqM_1\widetilde qq\right)=0.
> $$
>
> The traceless adjoint variation
>
> $$
> \delta Y=qM_1\widetilde q
> -\frac13\operatorname{Tr}(qM_1\widetilde q)\mathbf1
> $$
>
> is independent of $Y$, so its Konishi Jacobian with respect to $Y$ vanishes. Contracting with the adjoint F-equation yields
>
> $$
> -s_0\operatorname{Tr}(Y^2qM_1\widetilde q)
> +\frac{s_0}{\mu^2}\operatorname{Tr}\left[(qM_1\widetilde q)^2\right]=0,
> \qquad
> M_1\widetilde qY^2q=0.
> $$
>
> For the common variation $\delta q=Yq$ of all four flavors, the fundamental Konishi equation in the normalization of [Seiberg, hep-th/0212225, printed p.2, (2.3)](https://arxiv.org/abs/hep-th/0212225) gives
>
> $$
> \frac{s_0}{\mu^2}
> \left(M_1\widetilde qY^2q+M_2\widetilde qYq\right)
> =-\frac{4}{32\pi^2}\operatorname{Tr}(W_\alpha W^\alpha Y)=0.
> $$
>
> $W_\alpha$ is the magnetic gauge-field-strength spinor superfield; its spinor contraction is the source contraction. Therefore
>
> $$
> \operatorname{Tr}(W_\alpha W^\alpha Y)=0.
> $$

The electric variation $\delta Q=XQ$ similarly has zero superpotential variation and sets $\operatorname{Tr}(W_\alpha W^\alpha X)=0$ in such a definition. Thus the dressed-glueball objection has a concrete possible repair without assuming a dressed-glueball duality map. This does **not** fix the higher-trace definitions or justify combining classical Cayley–Hamilton with uncorrected Konishi relations. Seiberg's stated theory is $U(N)$ with fixed polynomial flavor masses; the present $SU(3)$ theory also requires traceless projection and dynamical singlets.

## 3. The unresolved characteristic-equation and operator-map bridge

[Ferrari, hep-th/0701220](https://arxiv.org/abs/hep-th/0701220), printed pp.6–7, §2, distinguishes two choices:

- Higher traces can be **defined** to obey classical characteristic identities in the quantum theory.
- A definition in which the generalized anomaly equations retain their perturbative form generally has quantum-corrected higher-trace identities.

These are not interchangeable conventions inside one derivation. Footnote 4 on p.6 separately acknowledges genuine deformation of nilpotents in the full chiral ring; higher-trace redefinition and full-ring deformation are different issues.

[Ferrari–Wens, arXiv:0710.2978](https://arxiv.org/abs/0710.2978), printed p.10, §2.3, states that nonperturbative corrections can be absorbed into operator definitions. It explicitly says the theorem was proved without flavors and takes its extension with flavors for granted. Their p.6, (2.17)–(2.21), uses a fixed matrix-valued mass polynomial, whereas $M_1,M_2$ here are dynamical singlets. Their polynomial parameter ring, p.9, (2.41), does not itself prove a statement about the full nilpotent special fiber of the present theory.

### Conditional coefficient diagnostic; not an established KSS correction

Use the source notation $u_j=\operatorname{Tr}Y^j$ for the adjoint moments, with $u_1=0$. Let $F(z)$ be the quantum characteristic function in Ferrari–Wens (2.12); $z$ is its spectral parameter. The first required terms are

$$
F(z)=z^3\exp\left(-\frac{u_2}{2z^2}-\frac{u_3}{3z^3}-\frac{u_4}{4z^4}+\cdots\right).
$$

Use their source mass determinant $U(z)=\det m(z)$, specialized formally to the magnetic flavor mass matrix:

$$
m(z)=\frac{s_0}{\mu^2}(M_1z+M_2),
\qquad
U(z)=\det\left[\frac{s_0}{\mu^2}(M_1z+M_2)\right].
$$

Suppose the quantum characteristic equation on p.12, (2.46), specializes to this $SU(3)$ theory and holds as a nonreduced operator identity:

$$
F(z)+\frac{\bar\Lambda^2U(z)}{F(z)}=H(z),
\qquad \deg H=3.
$$

Here $\bar\Lambda$ is KSS's magnetic strong scale, and $H$ is the source polynomial. Brackets $[z^j]$ extract the coefficient of $z^j$ in a Laurent expansion at infinity.

> **Proof.** Direct expansion gives
>
> $$
> [z^{-1}]F=\frac{u_2^2}{8}-\frac{u_4}{4},
> \qquad
> [z^{-1}]\frac{U}{F}
> =[z^2]U+\frac{u_2}{2}[z^4]U.
> $$
>
> The negative-power coefficient of $H$ vanishes, so
>
> $$
> u_4=\frac12u_2^2
> +4\bar\Lambda^2[z^2]U
> +2\bar\Lambda^2u_2[z^4]U.
> $$
>
> If compatible anomaly-clean F-relations also yield $u_4=u_2^2/3$, then
>
> $$
> u_2^2+24\bar\Lambda^2[z^2]U
> +12\bar\Lambda^2u_2[z^4]U=0.
> $$

Both extra terms vanish at $M_1=M_2=0$. Away from that locus, further singlet or nonperturbative relations might themselves set these determinant combinations to zero. **No nonzero correction and no full-ring counterexample have been proved here.**

The same characteristic equation also has

$$
[z]H=-\frac12u_2+\bar\Lambda^2[z^4]U,
\qquad
[z^4]U=\left(\frac{s_0}{\mu^2}\right)^4\det M_1.
$$

Thus even the relation between the quadratic trace and the degree-one characteristic coefficient must be specified. KSS's map (6.10) cannot silently be combined with a separately chosen resolvent scheme.

## 4. Symmetry and dimension check

Use KSS's conserved anomaly-free $R$ charge, pp.10 and 12, (2.6), (2.15), and engineering mass dimension, not an asserted infrared scaling dimension:

$$
R(X)=R(Y)=\frac23,\quad
R(Q)=R(\widetilde Q)=\frac16,\quad
R(q)=R(\widetilde q)=\frac12,
$$

$$
R(M_1)=\frac13,\quad R(M_2)=1,\quad R(W_\alpha)=1.
$$

The source normalization has dimensions $[X]=[Y]=1$, $[M_1]=2$, $[M_2]=3$, $[\mu]=[\Lambda]=[\bar\Lambda]=1$, $[s_0]=0$, and $[W_\alpha]=3/2$.

| Operator or candidate term | Engineering dimension | Conserved $R$ charge |
|---|---:|---:|
| $(\operatorname{Tr}Y^2)^2$ | $4$ | $8/3$ |
| $\operatorname{Tr}(W_\alpha W^\alpha Y)$ | $4$ | $8/3$ |
| $\bar\Lambda^2[z^2]U$ | $4$ | $8/3$ |
| $\bar\Lambda^2u_2[z^4]U$ | $4$ | $8/3$ |
| $\det M_1/\Lambda^6$ | $2$ | $4/3$ |
| $\operatorname{Tr}X^2$ | $2$ | $4/3$ |

KSS p.12, (2.16), retains its matching constant $C$:

$$
\Lambda^6\bar\Lambda^2=C\,s_0^{-8}\mu^8,
\qquad
\bar\Lambda^2\left(\frac{s_0}{\mu^2}\right)^4
=\frac{C}{s_0^4\Lambda^6}.
$$

Therefore symmetry and dimension allow mesonic contributions to both the squared relation and the trace map. They do not determine a coefficient or prove that any such term is nonzero.

## 5. What later literature verifies

[Ke Ye, arXiv:1706.02723](https://arxiv.org/abs/1706.02723), **PDF equation numbers** and printed pages:

- pp.17–19, (4.15)–(4.23): characteristic trace identities receive instanton corrections; the resolvent contains the flavor-mass determinant. The inference from all supersymmetric vacuum values on p.19 does not by itself fix nilpotents.
- p.32, §5.4: calls the KSS quantum Casimir relation a conjecture; checks mass-deformed $U(4)$ with $N_f=3,k=2$ by matching removed Coulomb vacua. This is not a verification of the undeformed $SU(5)$ square in the full operator ring.
- p.36, (A.8): explicitly distinguishes zero expectation values from zero chiral-ring elements and then omits generalized glueballs and photinos from the displayed ring.

The outstanding bridges are: the full renormalized trace map; compatible characteristic and Konishi definitions; $U(3)$ to traceless $SU(3)$; fixed mass coefficients to dynamical meson singlets; validity on the nilpotent special fiber; and the full two-supercharge physical chiral ring versus the manuscript's one-supercharge holomorphic Ward cohomology. The source prediction is a controlled candidate for nonperturbative lifting, not a completed microscopic lift calculation.
