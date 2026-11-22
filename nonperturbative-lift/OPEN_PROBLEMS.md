# Open problems

The first two problems ask for the variation of a specified microscopic local operator in the instanton sector, modulo perturbative variations. The five numbered problems are the open entries in [claims.json](claims.json). The former zero-mass torsion problem has been solved. [All-power proof](research/16-determinant-lift.md#massless-all-powers)

<a id="su2-ward-matching"></a>

## 1. Pure SU(2): finish the local Ward identity · `su2-ward-matching`

**Next equation.** In [the regulated product comparison](research/10-instanton-measure.md), $C_0$ is the gaugino bilinear, $\star$ the local quantum product, and $x_2$ the tree primitive of its ordinary square. The product and the supercharge variation have respective instanton coefficients $\nu_{\rm prod}$ and $\nu_Q$ relative to the same renormalized square. At one instanton the common subtraction cancels:

$$
C_0\star C_0-Qx_2=(\nu_{\rm prod}-\nu_Q)\Lambda^6.
$$

Here $\Lambda^6$ is the pure $SU(2)$ one-instanton weight. Start with $(Q_{\rm physical}-Q_L)c_\alpha$, the normal two-fermion variation of the gaugino insertion relative to the finite-cutoff collective-coordinate transformation; carry it through [the local insertion equation](research/13-connection-matching.md#np-b0509e2b7b2fc17a-1).

**Inputs.** Use [the Green response](research/11-green-response.md), [the transported measure and faces](research/12-finite-cutoff.md), and [the already integrated connection term](research/13-connection-matching.md#np-c209f4b75ff549bb-1). Include the second-order response of the pure-ghost component $M_0$, the remaining normal and antifield variations, action and measure derivatives, both instanton-size faces, the spatial face, and coincident contact terms. Compare the product and variation in the same local prescription.

**Completion test.** Give the invariant difference with its sign and normalization, or prove its cancellation, after the common subtraction and regulator limits. For a nonzero pairing, identify its nonzero target class modulo the [secondary-operation ambiguity](research/09-secondary-pairing.md). **Insufficient:** the nonzero separated integral, Wilson-line bound, finite connection term, or a vacuum condensate alone.

<a id="sqcd-microscopic-image"></a>

## 2. Massive SQCD: act on the polynomial current · `sqcd-microscopic-image`

**Next equation.** Electric $SU(2)$ with four flavors and one nonzero mass $m$ has a completely Higgsed magnetic $SU(2)$ description. In [the source construction](research/16-determinant-lift.md#np-f028df5dda216e10-1), $h$ is the heavy quark contraction, $A$ the three-light-row flavor current, and $\bar hK$ its polynomial Higgs-frame correction. Bars indicate the conjugate fields in the retained supercharge convention; $g$ is the gauge coupling, $\hbar$ counts loops, $\mu$ is the meson matching scale, and $\widehat M$ is the light meson matrix. The statement to test, with the magnetic instanton weight $\overline{\widetilde\Lambda}^{\,2}$ **included** on both sides, is

$$
\left[Q\left(\bar h A+\frac{3i\hbar g}{64\pi^2}\bar hK\right)\right]_{\text{one instanton}}
\ \stackrel{?}{=}\ -\frac{3\overline{\widetilde\Lambda}^{\,2}}{\bar\mu^3}\,[\det\overline{\widehat M}]
\quad\text{in perturbative cohomology}.
$$

**Inputs.** Start with the [local source derivative](research/16-determinant-lift.md#np-4f8eecadb9b2d576-1) at fixed magnetic light fields and the [triangle coefficient](checks/original/sqcd_duality_survey/fundamental_triangle_coefficient_audit.md). Match the source-dependent descendants, quantum products, contact subtractions, and any polynomial primitive for the difference in the same microscopic complex.

**Completion test.** Compute the local image and exhibit a finite polynomial representative for any perturbative boundary subtracted from it, or detect a nonzero residual class that changes the proposed equality. **Insufficient:** differentiating the effective superpotential gives the correct Wilsonian coefficient, but supplies no local operator map. The effective determinant result and the perturbative nonexactness of the incoming current already have separate proofs.

The [mass-pole Ward test](research/16-determinant-lift.md#mass-pole-ward) fixes additional input for this calculation. With nonnegative loop powers, the one-instanton, zero-loop action on $A$ must evaluate to $3\det\overline{\widehat M}/(\bar\mu^4\bar m)$. A mass-regular polynomial action on $A$ is excluded by the massive-family evaluation. A local calculation on the Higgs patch can test the residue and then the remaining, unevaluated contact terms.

<a id="su4-k6-primitive"></a>

## 3. SU(4): solve the missing first correction · `su4-k6-primitive`

**Next equation.** For the [explicit leading degree-six trace polynomial](evidence/su4_cartan/first_lift_bridge.md) $a_6$, solve $Q_0u=-Q_1a_6$ over the rationals at weight $(7,7)$. Here $u$ is the actual first correction. The historical search stores a *positive primitive* with the opposite sign.

**Inputs.** Use the [Cartan and nonabelian calculations](research/05-leading-modules.md#np-7d31a48e28a59909-1), [bounded-search record](evidence/su4_cartan/k6_first/bounded_status.json), and [source algorithm](sources/engines/originals/su4_cartan/lift_higher_first.py). **Completion test:** recompute the full rational polynomial $Q_0u+Q_1a_6$ and obtain zero with the proper grading and symmetries. A negative result requires an obstruction on the *whole* allowed source, for example a functional vanishing on its image but nonzero on $Q_1a_6$. **Insufficient:** vanishing on diagonal fields or exhaustion of the 169 searched directions and 133 retained image columns.

<a id="su5-kss-matching"></a>

## 4. Cubic-adjoint SU(5): retain the dynamical mesons · `su5-kss-matching`

**Next equation.** In [the Kutasov–Schwimmer–Seiberg comparison](research/20-cubic-adjoint-su5.md#np-1073436cbe3d7b68-1), $u_2=\operatorname{Tr}Y^2$ for the magnetic adjoint $Y$; $U(z)$ is the determinant of its flavor-mass matrix with meson entries, and $[z^j]$ extracts a polynomial coefficient. The bar on $\bar\Lambda$ here labels the **magnetic scale**. At zero quadratic deformation the proposed joint characteristic and anomaly relations would imply

$$
u_2^2+24\bar\Lambda^2[z^2]U+12\bar\Lambda^2u_2[z^4]U=0.
$$

**Inputs.** Establish the [operator map and singlet equations](sources/notes/duality_nonpert/new_critical_source_audit.md) in one quantum prescription before using this identity; reduce the last two terms with the mesons still dynamical. **Completion test:** determine the operator represented by the electric $(\operatorname{Tr}X^2)^2$ in the full ring. A microscopic instanton lift further requires its primitive and its instanton order. **Insufficient:** imposing the mesonic origin at the outset or inspecting only values on ordinary vacua, which cannot see a nonzero nilpotent class.

<a id="n4-nonscalar-image"></a>

## 5. N=4 Yang–Mills: a nonscalar instanton image · `n4-nonscalar-image`

**Next equation.** Choose a definite class $[\mathcal O]$ with fermions or derivatives, certify survival of the complete perturbative differential, then compute $d[\mathcal O]$. Here $d$ is the coefficient of $e^{2\pi i\tau}$ in the transferred supercharge, with $\tau$ the complex gauge coupling in [the source convention](research/22-n4-sym.md). No such class–target pair has yet been certified in this project.

**Inputs.** Fix the rank, charges, perturbative source and target bases, possible boundaries and mixing, and the complete instanton insertion. Exact diagonal vacua detect the specified pure-scalar ring. The [descendant calculation](research/22-n4-sym.md#n4-descendant-detector) additionally excludes all total holomorphic translations of the one-scalar half-BPS primaries. **Completion test:** exhibit a nonzero target class, or prove the proposed map vanishes on a further stated complete sector. A zero in one candidate sector does not settle the general existence question. **Insufficient:** a long-operator dimension shift, a mismatch of classical counts, or a vanishing linear instanton correction to the energy. A linear supercharge correction can first affect the energy quadratically.

## Two unfinished representatives

The following two tasks construct representatives of established perturbative results.

<a id="su2-weight-nine"></a>

### SU(2), weight $(9,9)$ · `su2-weight-nine`

The [leading block has stable dimension one](research/04-complex-and-products.md#np-9b603520cf852d74-1). In the [compact-chain record](evidence/rational_c4/m6_lift/STAGE3_STATUS.md), component $j$ is the coefficient of $\hbar^j$. The next right-hand side for $Q_0(\text{component }4)$ is $-Q_1(\text{component }3)-Q_2(\text{component }2)-Q_3(\text{component }1)-Q_4(\text{component }0)$. Choose its closed freedom to meet the fifth-order obstruction; then solve the fifth equation, including $Q_5(\text{component }0)$. **Inputs:** use the [published third polynomial](verification/README.md#external-witnesses) and the earlier components in the [calculation archive](artifacts/historical-files.json). Rebuild intermediate differential images as needed. The [historical fourth-component script](sources/engines/originals/rational_c4/m6_lift/advance_fourth.py) *expects* a fourth primitive; its existence is not evidence one was found. **Completion:** both full residuals vanish on independently reopened components; the assembled polynomial is closed. The dimension proof alone does not produce those components.

<a id="su3-k6-e3"></a>

### SU(3), degree six · `su3-k6-e3`

At weight $(7,7)$ the [joint solve](research/06-su3-lifts.md#np-0dbe41c29df07e6c-1) reaches $E_3$, the page after the first two corrections. The next target for $Q_0(\text{third correction})$ is $-Q_1(\text{second correction})-Q_2(\text{first correction})-Q_3(\text{leading term})$. **Inputs:** `archive:su3_rational/k6_two_loop/joint_two_loop_Q.pkl` and its [independent readback](evidence/su3_rational/k6_two_loop/independent_joint_two_loop_readback.json). Allow the earlier corrections to vary within their admissible freedom. **Completion:** a full rational primitive or a nonzero obstruction after quotienting that freedom decides the next page. All remaining orders require their own argument before claiming full survival; a diagonal residual alone does not decide it.

<a id="sqcd-massless-torsion"></a>

## Resolved SQCD massless torsion problem

For every $k\ge0$, $\hbar^k[C]\ne0$. In the complete source, the triangle row vanishes on the tree kernel, and derivative order excludes all higher loops from this projection. Flavor permutations give a free rank-sixteen submodule. The remaining instanton question for the massive current is [microscopic matching](#sqcd-microscopic-image). [Theorem 16.3](research/16-determinant-lift.md#massless-all-powers) · [Exact matrix](checks/massless_source.py)
