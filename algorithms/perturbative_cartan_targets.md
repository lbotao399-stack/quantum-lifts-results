# Cartan restriction of perturbative supercharge images

At fixed finite rank in undeformed four-dimensional $\mathcal N=4$ SYM, let $E_1=H(Q_0)$ be the classical cohomology and $d_r$ the order-$r$ supercharge obstruction, with $E_{r+1}=H(E_r,d_r)$. Write $\rho_C$ for evaluation of the leading representative with every letter in the Cartan subalgebra. In the multi-trace sectors of $SU(n+1)$, $SO(2n+1)$, $USp(2n)$ and $SO(2n)$, this evaluation descends to every perturbative page and

$$
\rho_C(d_r[F])=0 \qquad (r\geq 1).
$$

**Proof.** Take a classical cocycle $F$ representing an element of $E_r$. Write $\Delta$ for the characteristic-polynomial discriminant of the undifferentiated scalar $\Phi=\phi^1(0)$. At this rank there are an integer $k\geq 0$, a polynomial in multi-gravitons $G$ and a polynomial multi-trace $H$ such that

$$
\Delta^k F=G+Q_0H.
$$

This identity holds also when $F$ has nonzero Cartan restriction. The exponent may depend on $F$. Both $G$ and $\Delta$ are polynomials in gravitons. Replacing each product by the quantum product gives closed completions with those leading terms. Classical restriction kills $Q_0$-boundaries, since the classical letter variations are graded commutators. [Kibe, van Leuven and Mathieson, v3, §6.2, equations (98), (109)–(112)](https://arxiv.org/html/2610.22332v3#S6.SS2) · [Closed completions of multi-gravitons](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L4703)

Suppose Cartan evaluation is defined on $E_r$, as it is on $E_1$. The product on $E_r$ multiplies leading terms, and $d_r$ is its derivation. The two sides of $[\Delta]^k[F]=[G]$ agree first on $E_1$ and hence on each successive page reached by $F$. A closed completion has zero outgoing differential wherever its leading class is defined. Since $\Delta$ is even, on $E_r$ this gives

$$
[\Delta]^k d_r[F]=d_r\big([\Delta]^k[F]\big)=d_r[G]=0,
\qquad
\rho_C(\Delta)^k\,\rho_C(d_r[F])=0.
$$

On Cartan-valued jets the algebra is polynomial in the even variables with exterior generators for the odd ones. Each exterior monomial has its own polynomial coefficient. The nonzero polynomial $\rho_C(\Delta)$ acts injectively on each coefficient, also after restricting to Weyl invariants. Thus $\rho_C(d_r[F])=0$. Restriction on $\ker d_r$ descends through the quotient by $\operatorname{im}d_r$ to $E_{r+1}$. Induction proves the assertion at every finite order. The argument uses the quantum master equation through the product theorem. [Page product and derivation](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L3495)

In the trace sectors, the same calculation identifies the entire Cartan kernel on $E_r$. If $\rho_C(F)=0$, the classical identity has $G=0$ and $[\Delta]^k[F]=0$ already on $E_1$. Conversely, if $[\Delta]^k[F]=0$ on $E_r$, injectivity on Cartan polynomials gives $\rho_C(F)=0$. In particular all differentials vanish after inverting $[\Delta]$. A Cartan-visible class can still have a nonzero *outgoing* differential before this inversion. [Classical identity for zero Cartan restriction, equation (111)](https://arxiv.org/html/2610.22332v3#S6.SS2)

## Reflection-odd operators in the even orthogonal groups

For $SO(2n)$, let $P=\operatorname{Pf}(\Phi)$, odd under a reflection. Its classical variation is zero. Each $Q_r$ for $r\geq1$ lowers each of the three flavour degrees by $r$; two of those degrees of $P$ are zero. Thus $P$ itself is quantum closed. In the oriented Cartan basis with blocks $\left(\begin{smallmatrix}0&x_i\\-x_i&0\end{smallmatrix}\right)$ for $1\leq i\leq n$, its restriction is $x_1\cdots x_n$, again a nonzero scalar polynomial. [Supercharge grading](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L427)

If $F$ is reflection odd, $PF$ is reflection even. After polarizing the adjoint letters, the orthogonal invariant theorem expresses each reflection-even invariant as pairwise contractions. Each letter has two colour slots, so each pairing is a product of trace cycles. Equivalently, an even number of orthogonal volume tensors reduces in pairs to determinants of metric contractions. Antisymmetrization in odd letters preserves this span. Reflection commutes with the perturbative supercharge and product, since their colour contractions preserve the bracket and invariant form. The even part of every page is therefore the page of the trace sector. Inducting on the first page at which the odd restriction is needed gives

$$
0=\rho_C\big(d_r([P][F])\big)=\rho_C(P)\,\rho_C(d_r[F])
\quad\Longrightarrow\quad \rho_C(d_r[F])=0 .
$$

At the induction step the odd restriction exists on $E_r$ by the lower-order statements, and the displayed identity makes it descend to $E_{r+1}$. This includes the Pfaffian sector of $SO(8)$.

## The exceptional groups $G_2$ and $F_4$

The Cartan argument also applies to all polynomial colour invariants of these two groups. Here is the classical input for this extension. Work over characteristic zero. On the regular semisimple locus of $\Phi$, choose a Cartan and invert the root values $\alpha(\Phi)$. On the root space for $\alpha$, the bracket with $\Phi$ acts by that scalar. In the classical transformations this pairs the positive Taylor modes of the first scalar with the gauge ghost modes, the two other matter fermions with the two other scalars, and the field strength with the remaining matter fermion. The zero mode of the first scalar has already been put in the Cartan. Replace the image variable in each pair by the classical variation of its partner divided by the signed nonzero root value, retaining the numerical coefficients of the classical variations. This change is triangular: at fixed Taylor order recover the ghost first, then the two scalars, then the remaining fermion. Terms containing the nonconstant part of $\phi^1(z)$ have smaller Taylor order on the field being recovered. This is the root-space version of the four bracket identities and their inverse. Only the Lie bracket and invertibility of $\operatorname{ad}\Phi$ on the root spaces enter this step. [Classical transformations and triangular inversion, §6.2, equations (99)–(103)](https://arxiv.org/html/2610.22332v3#S6.SS2)

For each remaining Cartan letter $Y$, take the differentiated single-graviton trace components whose terms linear in $Y$ are nonzero multiples of $\operatorname{Tr}(\Phi^{d-1}Y)$. Choosing the coefficient with $d-1$ copies of the undifferentiated first scalar in the polarized trace gives precisely this linear term. The other terms on the Cartan involve at least two positive-weight letters, individually of lower weight, in the grading of equation (104). The matrix of the linear terms is the matrix of differentials of the scalar traces $\operatorname{Tr}(\Phi^d)$, up to nonzero row factors. Its invertibility gives closed replacements for the Cartan letters, recursively subtracting the expressions in already replaced lower-weight letters. These replacements agree with the original letters on the Cartan. The determinant check for each group follows below. [Construction of the closed coordinates, §6.2, equations (104)–(107)](https://arxiv.org/html/2610.22332v3#S6.SS2)

In the new coordinates the non-Cartan variables are differential pairs. On their monomials of positive total pair degree, replacing each image by its partner as a graded derivation and dividing by the total pair degree contracts the complex. Consequently a closed invariant equals a polynomial in the closed Cartan replacements modulo a classical boundary. Reexpressing the replacements by the invariant graviton traces gives coefficients regular on the regular Cartan. A torus-invariant primitive exists by taking its zero-weight part. Averaging over the finite Weyl group preserves the identity and makes the coefficients Weyl invariant. Restriction to the Cartan slice and extension by conjugation identify invariants over this locus: the regular semisimple locus is the associated bundle of the regular Cartan under its normalizer, with the remaining letter variables carried along. Regular invariant coefficients are polynomials in the basic scalar traces after inversion of the Weyl discriminant, the product of squared positive roots. Multiplication by a sufficient power clears every denominator, including those in the primitive, for each finite jet expression. This proves the classical identity used above for arbitrary invariant inputs. The graviton completions then give $\rho_C d_r=0$ by the same page induction. [Chevalley restriction and fundamental degrees](https://doi.org/10.1017/CBO9780511623646) · [Single-graviton traces and their quantum closure](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L4703)

For $G_2$ the seven-dimensional weights on $\Phi$ are $0,\pm x,\pm y,\pm z$ with $z=-x-y$. Write $e_2=xy+yz+zx$ and $e_3=xyz$ for the elementary symmetric polynomials. The two basic traces and their Jacobian are

$$
\operatorname{Tr}_7\Phi^2=-4e_2,\qquad
\operatorname{Tr}_7\Phi^6=-4e_2^3+6e_3^2,\qquad
\det\frac{\partial(\operatorname{Tr}_7\Phi^2,\operatorname{Tr}_7\Phi^6)}{\partial(x,y)}
=-48xyz(x-y)(x-z)(y-z).
$$

The six linear factors are the six positive-root hyperplanes. The characteristic-polynomial discriminant in the seven equals $64e_3^{10}(-4e_2^3-27e_3^2)^2$.

For $F_4$ use the Euclidean root coordinates of the level-30 calculation. The 26-dimensional weights are the 24 short roots and two zero weights. Take the scalar traces of degrees $2,6,8,12$. Their Jacobian is alternating under each Weyl reflection and homogeneous of degree $1+5+7+11=24$. Each of the 24 positive-root linear forms divides it. For positive roots chosen by the first nonzero coordinate, at $(1,2,3,5)$ the determinant and root product are $329986597248000000$ and $-5893965000$. Thus their ratio is the constant $-55987200$, and the Jacobian is invertible everywhere on the regular Cartan. The independent homogeneous traces have precisely the fundamental degrees in each group, so their polynomial subalgebra has the full invariant-ring Hilbert series. In $F_4$ the repeated zero weight is retained throughout; the discriminant used for localization is the Weyl discriminant. [Exact trace and Jacobian replay](exceptional_discriminants.py) · [Rational values](../data/N4_exceptional_discriminants.json) · [Root coordinates used in the campaign](dual_singlets_20261118/f4_bmn_root_restriction.py)

## Three second-order target components

Here $j=2(E+J_L)$ is the physical level and $Y$ counts letters. For the three components below, independent Cartan product evaluations attain the rational upper bounds from the two adjacent classical maps. Thus $\ker\rho_C=0$ there already on $E_1$. Earlier images restrict to zero, so every later page in the component still injects into Cartan polynomials. If the filtered Taylor selection identifies the physical target with the listed component, its incoming $d_2$ has rank zero.

| Gauge group | $j$ | $Y$ | Classical dimension | Cartan rank | Possible incoming $d_2$ rank |
|---|---:|---:|---:|---:|---:|
| $SU(3)$ | 30 | 9 | $3683-2174-1491=18$ | 18 | 0 |
| $SU(3)$ | 32 | 10 | $6828-4185-2615=28$ | 28 | 0 |
| $USp(4)$ | 36 | 12 | $46860-28039-18760=61$ | 61 | 0 |

Cartan restriction commutes with the classical Taylor projections. A projection of a zero Cartan polynomial is again zero. The table accounts for $18+28+61=107$ target directions in three distinct sectors. The proof of the all-orders restriction itself uses no Taylor selection rule. [Twenty-seven Cartan matrices and classical bounds](dual_singlets_20261118/taylor_hooks/cartan_products/README.md)

At $SU(3)$ level 30 the independent Fortuity analysis identifies eight Cartan-detected multiplets after its full lowest-weight test. Its one Cartan-kernel direction fails that test. Combining its selection rule with the result above excludes $d_2$ for its BMN threshold class R-007. Its analysis separately excludes orders three and four and higher targets, and reports survival at order one. For level 32 the full first-correction equation remains open. The formal perturbative product and discriminant argument gives no instanton matrix element. [Fortuity, perturbative lifting of the $SU(3)$ threshold class](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex#L455)

For the three source flavour weights $(7,4,4)$, $(8,4,4)$ and $(10,4,4)$ the third-order tail entries in the same polynomial Taylor highest weight are $(2,-1,-1)$, $(3,-1,-1)$ and $(5,-1,-1)$. One negative entry excludes the ordered partner, and its value decreases at each subsequent order. An incoming perturbative image in a BMN sector would require a negative spacetime derivative degree in its source. Hence for the indicated $SU(3)$ level-32 and $USp(4)$ level-36 BMN sources, a positive first differential is the only remaining perturbative lifting mechanism under this Ward selection. [Taylor highest weights](dual_singlets_20261118/b2_schur_complex.py) · [Symplectic source weights and lengths](dual_singlets_20261118/b2_level36/README.md)
