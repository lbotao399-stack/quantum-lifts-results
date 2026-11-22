# Four-dimensional N=1 SQCD: nonperturbative relations and an explicit differential

The primary inputs are Seiberg, [hep-th/9402044](https://arxiv.org/abs/hep-th/9402044), equations (3.3), (3.4), (4.1), (4.8), (5.4), and Seiberg, [hep-th/9411149](https://arxiv.org/abs/hep-th/9411149), equations (4.13)–(4.14). Their PDFs and extracted text are saved in `sources/`.

The exact coefficients below use the first paper's holomorphic scale convention, fixed immediately before (3.3). No coefficient is transferred to the pure-SYM source's $C_0$ or instanton measure. The F-term/Koszul calculations are performed here; they do not assume a chain equivalence between the full electric and magnetic holomorphic-jet complexes.

## 1. Three distinct outcomes

For $SU(N_c)$ SQCD, $N_f$ counts pairs of fundamental/antifundamental quark multiplets. The source meson is $M^i{}_{\tilde\jmath}=Q^i\widetilde Q_{\tilde\jmath}$, and $B,\widetilde B$ are its baryons. The untilded and tilded flavor slots belong to the two flavor groups.

| Theory | Exact source input | What follows |
|---|---|---|
| $N_f=N_c$ | $\det M-B\widetilde B=\Lambda^{2N_c}$, (4.1) | The classical singular locus is excluded; a quantum moduli space remains. This is a ring deformation, not by itself a differential lifting theorem. |
| $N_f=N_c+1$ | $W_{\mathrm{eff}}=(B_iM^i{}_{\tilde\jmath}\widetilde B^{\tilde\jmath}-\det M)/\Lambda^{2N_c-1}$, (5.4) | The F terms reproduce the electric classical constraints. The massless moduli space is not removed. |
| $N_f<N_c$ | $W_{\mathrm{dyn}}=(N_c-N_f)(\Lambda^{3N_c-N_f}/\det M)^{1/(N_c-N_f)}$, (3.4) | The massless theory has no finite supersymmetric vacuum. Its runaway must not be called a stable vacuum with positive vacuum energy. |

These statements are the source's separate cases; a nonzero nonperturbative superpotential does not universally imply fewer scalar chiral-ring classes. [Seiberg 9402044, §§3–5](https://arxiv.org/pdf/hep-th/9402044).

## 2. SU(2), Nf=2: deformation and smoothness

In the enhanced $SU(4)$ flavor presentation of the same source, $V^{ij}=Q^iQ^j=-V^{ji}$, with flavor components $i,j=1,\ldots,4$. Its Pfaffian is

$$
\operatorname{Pf}V=V^{12}V^{34}-V^{13}V^{24}+V^{14}V^{23}.
$$

Equation (4.1) gives the two coordinate rings

$$
\frac{\mathbb C[V^{12},V^{13},V^{14},V^{23},V^{24},V^{34}]}
{(V^{12}V^{34}-V^{13}V^{24}+V^{14}V^{23})},
$$

$$
\frac{\mathbb C[V^{12},V^{13},V^{14},V^{23},V^{24},V^{34}]}
{(V^{12}V^{34}-V^{13}V^{24}+V^{14}V^{23}-\Lambda^4)}.
$$

**Proof of the changed singular locus.** In the displayed component order, the six derivatives of the Pfaffian are

$$
\left(V^{34},-V^{24},V^{23},V^{14},-V^{13},V^{12}\right).
$$

They vanish simultaneously exactly at $V=0$. For $\Lambda\ne0$, that point does not satisfy $\operatorname{Pf}V=\Lambda^4$. The quantum hypersurface is therefore smooth of complex dimension $6-1=5$. On its chart $V^{12}\ne0$,

$$
V^{34}=\frac{\Lambda^4+V^{13}V^{24}-V^{14}V^{23}}{V^{12}}.
$$

This removes the classical origin, not all vacua. No odd partner or nonperturbative differential follows merely from writing these two quotient rings.

## 3. One mass gives an exact ADS calculation

Use the source's meson/baryon presentation of (4.8). For $N_c=N_f=2$, $X$ is its Lagrange multiplier and $m_2$ is the coefficient of the second diagonal meson mass term:

$$
W=X\left(
M^1{}_{\tilde1}M^2{}_{\tilde2}
-M^1{}_{\tilde2}M^2{}_{\tilde1}
-B\widetilde B-\Lambda^4
\right)+m_2M^2{}_{\tilde2}.
$$

**Proof by integrating out the heavy fields.** Keep $M^1{}_{\tilde1}$ as the light meson. The other F terms give

$$
XM^1{}_{\tilde1}+m_2=0,
\qquad XM^1{}_{\tilde2}=XM^2{}_{\tilde1}
=XB=X\widetilde B=0.
$$

For $m_2\ne0$, $X$ and $M^1{}_{\tilde1}$ are nonzero. Consequently

$$
X=-\frac{m_2}{M^1{}_{\tilde1}},
\quad M^1{}_{\tilde2}=M^2{}_{\tilde1}=B=\widetilde B=0,
\quad M^2{}_{\tilde2}=\frac{\Lambda^4}{M^1{}_{\tilde1}}.
$$

Calling the sole remaining source meson $M$, and using the source's decoupled scale $\Lambda_L$,

$$
\Lambda_L^5=m_2\Lambda^4,
\qquad W_{\mathrm{eff}}=\frac{m_2\Lambda^4}{M}
=\frac{\Lambda_L^5}{M}.
$$

The light F term was not used when integrating out the heavy fields. It is

$$
\frac{\partial W_{\mathrm{eff}}}{\partial M}
=-\frac{\Lambda_L^5}{M^2}\ne0
\qquad(0<|M|<\infty,\ \Lambda_L\ne0).
$$

No omitted $M=0$ solution exists: the unintegrated equation $XM+m_2=0$ excludes it. Along a D-flat configuration with quark expectation values $Q=\widetilde Q=(v,0)$, the large-field canonical contribution has

$$
M=v^2,
\qquad
\left|\frac{\partial W}{\partial Q}\right|^2
+\left|\frac{\partial W}{\partial\widetilde Q}\right|^2
=\frac{2|\Lambda_L|^{10}}{|v|^6},
\qquad
\lim_{|v|\to\infty}\frac{2|\Lambda_L|^{10}}{|v|^6}=0.
$$

This last expression is the canonical-metric contribution, not an assertion of an exact quantum Kähler metric. The exact statement used here is the absence of a finite F-flat solution and the source's runaway. The superpotential and scale matching follow from (3.4) and §4.1 of [Seiberg 9402044](https://arxiv.org/pdf/hep-th/9402044).

## 4. The same example as an actual Koszul differential

The original manuscript `/Users/libotao/np_lifts_engine/work/allloop_repo/main.tex`, `eq-CY-completion-images`, gives, after gauge-singlet projection,

$$
Q\gamma^{(0)}=0,
\qquad Q\beta^{(0)}=-\frac{\partial\bar W}{\partial\gamma^{(0)}}.
$$

Here $\gamma^{(0)}=\bar M$ is the retained even meson letter, and $\beta^{(0)}$ is its retained odd partner. The complex conjugate is required because that source retains the antichiral superpotential. This uses its existing coefficient and sign; it does not identify Seiberg's holomorphic $W$ directly with $\bar W$.

For the exact ADS superpotential and the localized patch $\gamma^{(0)}\ne0$,

$$
\bar W=\frac{\bar\Lambda_L^5}{\gamma^{(0)}},
\qquad
Q\beta^{(0)}=\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2},
\qquad
Q\left((\gamma^{(0)})^2\beta^{(0)}\right)=\bar\Lambda_L^5.
$$

The Koszul complex means the even Laurent-polynomial functions together with one odd exterior generator $\beta^{(0)}$, with the displayed differential. Give the even functions degree $0$ and $\beta^{(0)}$ degree $-1$. Over the formal coefficient ring its entire differential is

$$
\mathbb C[[\bar\Lambda_L^5]][(\gamma^{(0)})^{\pm1}]\,\beta^{(0)}
\xrightarrow{\ Q\ }
\mathbb C[[\bar\Lambda_L^5]][(\gamma^{(0)})^{\pm1}],
\qquad
f\beta^{(0)}\longmapsto
f\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2}.
$$

**Proof of the cohomology.** This coefficient is a nonzero divisor, so the degree $-1$ kernel is zero. Since $\gamma^{(0)}$ is invertible, its image is exactly the ideal $(\bar\Lambda_L^5)$. Thus

$$
H^{-1}=0,
\qquad
H^0=
\frac{\mathbb C[[\bar\Lambda_L^5]][(\gamma^{(0)})^{\pm1}]}
{(\bar\Lambda_L^5)}
=\mathbb C[(\gamma^{(0)})^{\pm1}].
$$

This last cohomology is annihilated by $\bar\Lambda_L^5$; it is not a family of finite vacua at nonzero strong scale. At $\bar\Lambda_L=0$ the differential itself is zero, and the odd partner also survives. At fixed nonzero strong scale, or after inverting $\bar\Lambda_L^5$,

$$
Q\left(\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}\right)=1.
$$

For any cochain $\mathcal O$, odd Leibniz gives the explicit contracting homotopy

$$
Q\left(\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}\mathcal O\right)
+\frac{(\gamma^{(0)})^2\beta^{(0)}}{\bar\Lambda_L^5}Q\mathcal O
=\mathcal O.
$$

This is a true differential lifting in the stated effective, localized Koszul/BV complex. It does not by itself identify that complex with the full electric local-operator complex, allow inverse mesons in the UV polynomial algebra, or establish the full holomorphic-jet duality map.

## 5. A second mass restores two isolated vacua

Add the remaining source mass $m_1M$:

$$
W_{\mathrm{eff}}=\frac{\Lambda_L^5}{M}+m_1M,
\qquad
M^2=\frac{\Lambda_L^5}{m_1},
\qquad
M=\pm\sqrt{\frac{\Lambda_L^5}{m_1}},
$$

$$
W_{\mathrm{eff}}\big|_{\mathrm{vacua}}
=\pm2\sqrt{m_1\Lambda_L^5}
=\pm2\sqrt{m_1m_2\Lambda^4}.
$$

The same source differential now reads

$$
Q\beta^{(0)}
=\frac{\bar\Lambda_L^5}{(\gamma^{(0)})^2}-\bar m_1.
$$

**Proof.** At fixed nonzero $m_1,\Lambda_L$, multiplication by this Laurent polynomial has no kernel. Its degree-zero quotient is

$$
H^{-1}=0,
\qquad
H^0=
\frac{\mathbb C[(\gamma^{(0)})^{\pm1}]}
{((\gamma^{(0)})^2-\bar\Lambda_L^5/\bar m_1)}.
$$

The two roots are distinct and nonzero. This ring has dimension two over $\mathbb C$, with basis $1,\gamma^{(0)}$. Thus the massless runaway and the two massive vacua are different cases of the same exact differential, not interchangeable descriptions of one vacuum.

## 6. SU(2), Nf=3: the Pfaffian F terms are mandatory but do not erase the moduli space

Now the source $V^{ij}$ is antisymmetric with $i,j=1,\ldots,6$. Equation (5.4) fixes

$$
W_{\mathrm{eff}}=-\frac{\operatorname{Pf}V}{\Lambda^3}.
$$

**Proof of a displayed F term.** Expansion of the Pfaffian along its first flavor gives

$$
\begin{aligned}
\operatorname{Pf}V={}&
V^{12}(V^{34}V^{56}-V^{35}V^{46}+V^{36}V^{45})\\
&-V^{13}(V^{24}V^{56}-V^{25}V^{46}+V^{26}V^{45})\\
&+V^{14}(V^{23}V^{56}-V^{25}V^{36}+V^{26}V^{35})\\
&-V^{15}(V^{23}V^{46}-V^{24}V^{36}+V^{26}V^{34})\\
&+V^{16}(V^{23}V^{45}-V^{24}V^{35}+V^{25}V^{34}).
\end{aligned}
$$

Consequently

$$
\frac{\partial W_{\mathrm{eff}}}{\partial V^{12}}
=-\frac{V^{34}V^{56}-V^{35}V^{46}+V^{36}V^{45}}{\Lambda^3}=0.
$$

The other fourteen independent derivatives similarly impose all $4\times4$ Pfaffians. For an antisymmetric matrix these equations are equivalent to

$$
\operatorname{rank}V\le2.
$$

They reproduce the UV decomposable-bivector relations of $V^{ij}=Q^iQ^j$. The odd IR partner of $V^{12}$ therefore has a nonzero differential equal to the complex conjugate displayed quadratic divided by $\bar\Lambda^3$, with the positive sign fixed in §4. Omitting the dynamically generated superpotential leaves that quadratic incorrectly unpaired in the unconstrained IR singlet model. It was already zero in the electric scalar chiral ring.

At $V=0$, both the first and second derivatives of the cubic superpotential vanish. The origin remains a vacuum with all fifteen composite chirals massless. Thus this example supplies a necessary nonperturbative IR differential without claiming that the electric massless moduli space is destroyed. [Seiberg 9402044, (2.5), (5.4)](https://arxiv.org/pdf/hep-th/9402044).

## 7. Why the magnetic instanton term cannot be omitted

Start with electric $N_f=N_c+2$, whose magnetic gauge group is $SU(2)$. A flavor mass completely Higgses that magnetic gauge group. In the source variables, its broken-group instanton contributes the second term in

$$
W=\widehat M^i{}_{\tilde\jmath}\widehat q_i\widehat{\widetilde q}^{\tilde\jmath}
-\frac{\det\widehat M}{\Lambda^{2N_c-2}m}.
$$

After the source's baryon identification this is exactly

$$
W=\frac{B_i\widehat M^i{}_{\tilde\jmath}\widetilde B^{\tilde\jmath}
-\det\widehat M}{\Lambda_L^{2N_c-1}},
\qquad
\Lambda_L^{2N_c-1}=m\Lambda^{2N_c-2}.
$$

The coefficient and the baryon mapping here are those of [Seiberg 9411149, (4.13)–(4.14)](https://arxiv.org/pdf/hep-th/9411149), not a newly chosen normalization. The $SU(6)$ Pfaffian presentation in §6 and this meson/baryon presentation are kept in their respective source conventions; no unrecorded flavor-ordering sign is used to identify their components.

**Proof of the spurious branch without the term.** At $B=\widetilde B=0$, the tree term alone has zero derivatives for every $\widehat M$. The full superpotential instead requires

$$
\frac{\partial\det\widehat M}
{\partial\widehat M^i{}_{\tilde\jmath}}=0.
$$

For $N_c=2$ this is the vanishing of all $2\times2$ minors of the $3\times3$ meson matrix, hence $\operatorname{rank}\widehat M\le1$ on this baryon-zero locus. In particular a full-rank diagonal $\widehat M$ is excluded. This is an explicit removal of false classical magnetic vacua by a required nonperturbative term; it restores the electric moduli-space constraints. The full baryonic locus obeys the cofactor relation, not the baryon-zero rank bound everywhere.

## 8. Most actionable extension of the pure-SYM mechanism

The controlled sequence is the $SU(2)$ s-confining theory, its quantum-deformed $N_f=2$ descendant, and the mass-deformed ADS theory. It supplies three different tests:

1. Construct the effective singlet BV differential from the known Pfaffian superpotential and check the exact Pfaffian F-term relations, including their syzygies. This is an explicit target for an electric/magnetic chain map, not just a matching of scalar dimensions.
2. Under holomorphic mass decoupling, transport the explicit relation
   $$
   Q\big((\gamma^{(0)})^2\beta^{(0)}\big)=\bar\Lambda_L^5
   $$
   and distinguish formal strong-scale torsion from cohomology after inverting the scale.
3. With both masses nonzero, compare the two-dimensional effective cohomology with the pure-SYM sector using the existing source operator and scale dictionary. This does not determine the current pure-SYM $C_2$ partner coefficient without the physical composite/BV and measure matching.

The new calculation here is the explicit differential/decoupling comparison. The SQCD superpotentials and quantum moduli spaces are established primary-source inputs. Extending the comparison to the full one-supercharge holomorphic-jet complex, its products, and its chain homotopies remains a research task.
