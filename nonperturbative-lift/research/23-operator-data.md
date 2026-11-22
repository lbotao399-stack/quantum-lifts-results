<a id="np-section-21"></a>

## 23. Discussion: nonperturbative effects and operator data

The protected-operator problem requires retaining the differential, product, and coefficient ring together.

An ordinary effective superpotential gives the scalar F-equations at zero derivative order.

Transport to the microscopic one-supercharge jet complex requires a local operator map.

| Object | Calculated effect | Conditions |
|---|---|---|
| Pure SYM quantum product | Nonzero hidden product, nilpotency, and explicit partners | Anomaly-free perturbative product construction |
| Pure SYM secondary operation | Conditional formula for the instanton differential | Exact multiplication homotopy and local composite matching |
| Magnetic SQCD determinant | The entire scalar torsion ideal is the determinant times the meson polynomials | Fixed magnetic quark coordinates on the massive Higgs patch |
| Massless magnetic SQCD | Sixteen microscopic current components survive all loop powers | Cubic relative jet complex, two opposite color evaluations of the full one-loop source |
| Symplectic Pfaffian | All submaximal Pfaffians and their polynomial multiples; $105$ first fermionic targets at six flavors | Fixed electric meson coordinates; the fermionic statement uses the ungauged holomorphic effective model |
| SQCD rank patches | ADS dynamics excludes excess-rank stationary vacua | Holomorphic thresholds and correct condensate sheets |
| Quantum constraints | Deformation of classical relations with retained vacuum branches | Specified gauge group, flavor count, and scale normalization |
| $SO$ branches | Instanton and condensate cancellation retains a flat branch | Branch signs and low-rank scale conventions |
| Higher F-terms | Nontrivial interaction class and mass descent | Globally admissible primitives and fermion zero modes |
| Circle monopole ring | Every scalar torsion class is in $(Y)$, supported at rank at most $N_c-2$ and zero baryons | Three-dimensional disorder operators and the finite-circle term |
| $\mathcal N=2$ periods | Central charges, singular loci, and prepotential corrections | Source curve, period basis, and coupling coordinate |
| Pure $\mathcal N=2$ bulk jets | Two full derivative towers are protected | Relative basic complex, spin and ghost grading, holomorphic translation Ward identity |
| $\mathcal N=4$ scalar quotient | Exact diagonal vacua detect nonzero classes | Reduced invariant quotient and vacuum Ward identity |
| $\mathcal N=4$ long multiplets | Nonzero instanton anomalous dimension | Specified multiplet and complete zero-mode saturation |

When a positive power of a formal parameter times a class becomes a boundary, the class has torsion.

Its vanishing after inverting that parameter is a statement over the localized coefficient ring.

Mass-threshold elimination also takes place on a specified open patch.

The zero-mass fiber requires a new calculation.

Each exact identity therefore retains its theory, operator sector, parameter domain, and completed perturbative page.

<a id="complete-tests"></a>

### Complete tests when the operator differential is known

For a polynomial scalar presentation let $I$ be the ideal of F-equations in the parameter and field coordinates and let $\eta$ denote the parameter being tested. Saturation means taking every polynomial $f$ for which $\eta^r f\in I$ at some $r\ge0$. The kernel after allowing the inverse parameter is

$$
\ker\!\left(\mathbb C[\eta,\mathrm{fields}]/I
\longrightarrow \bigl(\mathbb C[\eta,\mathrm{fields}]/I\bigr)[\eta^{-1}]\right)
=\frac{I:\eta^\infty}{I}.
$$

> **Proof.** The defining equivalence of fractions sets $f/1=0/1$ precisely when an allowed denominator annihilates $f$ modulo $I$. At finite exponent $r$ the analogous kernel is $(I:\eta^r)/I$. Eliminating an auxiliary variable $u$ from $(I,u\eta-1)$ computes $I:\eta^\infty$, since its quotient is the localization. The determinant calculation proves equality of the entire saturated ideal at every matrix size. The two confined presentations show why the regular field coordinates must be fixed before this operation. $\square$

For a finite free operator sector over the one-parameter power series ring $\mathbb C[[\eta]]$, there is a corresponding exact test at every fermion degree. Use its outgoing differential to take the kernel in the target degree. Express the incoming differential in any basis of this kernel. Diagonalizing this matrix by invertible row and column operations over $\mathbb C[[\eta]]$ is Smith normal form. Each nonzero diagonal entry is a unit times $\eta^r$ for a unique exponent $r$ up to ordering. A zero row left in the target gives a free class.

> **Proof.** The ring is a discrete valuation ring, and a submodule of a finite free module is free. Nilpotency puts the incoming image in the outgoing kernel. The quotient by a diagonal image decomposes into free summands and the quotients $\mathbb C[[\eta]]/(\eta^r)$. An exponent $r=0$ gives a boundary already at zero weight. For $r>0$ the isolated incoming generator has differential $\eta^r$ times its target. There is no earlier differential between this pair, and its coefficient on page $r$ is a unit. Equivalently, if the least parameter exponent of a nonzero $j$-rowed minor is $v_j$, with $v_0=0$, then the ordered exponents of the nonzero diagonal entries are $v_j-v_{j-1}$. Invertible row and column operations preserve the ideals generated by minors, and the diagonal form gives the formula. Localizing makes every finite-exponent target exact. $\square$

This test requires the matrix of the actual differential, or enough exact coefficients and bounds to determine the relevant minors. Two independent parameters allow nonprincipal ideals such as $(\hbar,\eta)$ in $\mathbb C[[\hbar,\eta]]$, precluding this general diagonalization. The mass-pole calculation instead tests the two expansions before specialization. The multiplier superpotential calculates invariant exponents one and two explicitly. The bulk translation calculation excludes an incoming image without knowing its matrix element. [Mass-pole test](16-determinant-lift.md#mass-pole-ward) · [Multiplier calculation](19-higher-f-terms-and-monopoles.md#multiplier-pages) · [Translation test](21-seiberg-witten.md#bulk-translation)

The microscopic SQCD family realizes a nonprincipal annihilator for the quadratic ghost trace, generated by the mass and the loop parameter. The relation between these two variations is a free incoming class with a nonzero image when cohomology is reduced modulo both parameters. Its cochain becomes exact in the jointly specialized differential. The second differential of the two-parameter Koszul resolution computes that discrepancy. [Microscopic base-change calculation](16-determinant-lift.md#mass-loop-base-change)

<a id="index-pairs"></a>

**Index test.** In a finite sector graded by charges commuting with the differential, the signed count of states is the even dimension minus the odd dimension. For the six-flavor-slot Pfaffian model assign degree one to each scalar and degree two to each odd partner. The differential preserves this degree because it replaces one odd partner by two scalars. At degree four, the products of two odd partners span an even space of dimension $105$ and their one-instanton images span an odd space of dimension $105$. Together these two spaces form a subcomplex, with signed dimension given on either fiber by

$$
\chi_{\eta=0}=105-105=0,\qquad
\chi_{\eta\ne0}=0-0=0 .
$$

Here $\eta$ denotes the weight in the index test. At zero weight all $210$ states are classes of the subcomplex. At fixed nonzero weight the map between its two halves is an isomorphism. Thus even exact knowledge of this signed character is compatible with $105$ independent pairings. In general an odd differential removes one even and one odd state per independent pair, giving $\dim H_{\rm even}+\dim H_{\rm odd}\ge|\chi|$. The operator matrices determine the pairings compatible with that bound. [Explicit Pfaffian images](18-sp-and-so.md#pfaffian-jet-images)

Seiberg's Higgsed instanton fixes the determinant interaction. Intriligator and Pouliot fix the symplectic interaction. Gray and collaborators compute classical scalar invariants and their relations. The results here calculate the kernels as the magnetic weight passes through zero, including their dependence on the regular field coordinates. The one-loop transformation uses the generalized Konishi framework studied by Budzik and Kulp. The paired color evaluations control all possible microscopic primitives of the stated class. Holomorphic Seiberg duality at the level of local observables remains a separate physical identification in the formulation of Gwilliam and Williams. [Seiberg](https://arxiv.org/abs/hep-th/9411149) · [Intriligator–Pouliot](https://arxiv.org/abs/hep-th/9505006) · [Gray et al.](https://arxiv.org/abs/0803.4257) · [Budzik–Kulp](https://arxiv.org/abs/2512.07771) · [Gwilliam–Williams, §7](https://doi.org/10.1112/blms.70152)

Two local calculations connect the main results: the full Ward matching of the original $SU(2)$ secondary partner, and the instanton image of the SQCD polynomial incoming class.

They would realize the established quantum-product structure and exact effective determinant relation, respectively, as nonzero matrix elements of a microscopic operator differential.
