# The $F_4$ fifth-power channels

The terminal fifth gaugino power has level $j=2(E+J_L)=30$. Its incoming spin and flavour singlets have doubled charges $(1,1,10,10,10)$. At derivative bidegree $(t,t)$ their three flavour degrees are $5-t$ and their length is $14-t$. The first, second and third perturbative sources for the pure fifth power have $t=4,3,2$, respectively. The following are full invariant cochain dimensions before taking classical cohomology.

| $t$ | Length $13-t$ | Source, length $14-t$ | Target, length $15-t$ |
|---:|---:|---:|---:|
| 5 | 31 | 20 | 12 |
| 4 | 2053 | 745 | 144 |
| 3 | 18533 | 6862 | 1278 |
| 2 | 46192 | 17499 | 3092 |
| 1 | 31394 | 12353 | 2333 |
| 0 | 4295 | 1862 | 428 |

The integral character algorithm multiplies adjoint Adams operations by reflecting shifted highest weights into the dominant chamber. Symmetric and exterior invariant dimensions through degree 12 agree with the independent generator-degree formulas. Nine $B_4$ cochain calculations agree with the orthogonal symmetric-group calculation, including each separate fermion degree. [Character algorithm](../f4_singlet_characters.py) · [All cochain degrees and checks](characters.json)

A8 job `5670599` used all 52 independent derivations of the Albert algebra, acting on its 27-dimensional space of Hermitian octonionic matrices. The scalar identity line contributes zero to each nonempty product. Each generator satisfies the derivation identity on all $27^2$ basis pairs and preserves the trace metric. [Colour matrices](../f4_colour.py) · [Albert-algebra realization](https://arxiv.org/abs/math/0105155)

At $t=4$, evaluating the first classical map on these trace expressions gives a 134-dimensional image at the recorded prime. The full source has dimension $666+79=745$, so its rational classical kernel has dimension at most $745-134=611$. In formal orthogonal traces the outgoing map has rank 358 over the rationals. The direct sum of the finite-colour source evaluation and the formal classical map gives a minor of order 948. It follows that at least $948-358=590$ independent rational closed directions are restrictions of formal classical cycles. The kernel modulo these cycles has dimension at most $611-590=21$. Counting possible lifts after removing these directions uses one-loop protection of universal multi-graviton classes. The quotient estimate itself uses only classical maps. [Projection and determinants](ranks.json) · [A8 receipt](a8_5670599.json) · [First-order protection](https://arxiv.org/abs/2610.22332)

From the repository root, using the formal source plan from the [$SO(9)$ calculation](../rank_four/README.md), run

```sh
python3 algorithms/dual_singlets_20261118/f4_singlet_characters.py --out /tmp/f4-characters.json
python3 algorithms/dual_singlets_20261118/f4_formal_cycles.py --plan /tmp/so9-s5/plan.json.gz --work /tmp/f4-singlets --workers 32
```

## The classical target and its pure-gaugino predecessor

Write $S=\operatorname{Tr}_{27}(\lambda_+\lambda_-)$ for the bilinear in the two gaugino components in the Albert representation. Its scalar identity line contributes zero. At derivative bidegree $(5,5)$, with all three flavour degrees zero, the singlet dimensions and classical ranks over the rationals are

| Number of letters | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cochains | 1 | 2 | 18 | 62 | 99 | 74 | 31 | 20 | 12 |
| Incoming rank | 0 | 1 | 1 | 17 | 45 | 52 | 22 | 9 | 11 |
| Classical cohomology | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 1 |

Every differential meets a length at which its rank lower bound and that of the adjacent differential sum to the character dimension. The identity $Q_0^2=0$ forces equality at both maps. At length ten, integer evaluations of all 45 formal predecessor images have rank 11. Appending $S^5$ raises the rank to 12. The predecessor evaluations have full cochain rank 20. The saved rational functional vanishes on all 45 images and takes value one on $S^5$. The 96 target evaluations also agree with independent exterior Pfaffians. At length nine the entire classical singlet source has dimension $20-9-11=0$. [Integer target and functional](pure.json.gz) · [Character dimensions and minor certificates](pure_chain.json)

## Exceptional colour contractions in the first-order source

At derivative bidegree $(4,4)$ two open words in the 27-dimensional representation can be paired through either the 52 derivations or the 26 traceless Jordan multiplication operators. The latter span an invariant subspace because a derivation acts on multiplication by $x$ through multiplication by its action on $x$. For each subspace the contraction uses the inverse trace Gram matrix. Reversing an open word includes the transpose sign of its insertion, negative for a derivation and positive for Jordan multiplication. Each sampled parent is a spin and flavour singlet. Its classical variation was formed before evaluating any colours.

| Insertions | Rank detected in the length-eleven image | Lower bound on formal cycles plus nine-letter boundaries | Upper bound on the remaining length-ten kernel |
|---|---:|---:|---:|
| 52 derivations | 134 | 590 | 21 |
| 26 Jordan multiplications together with the ordinary trace rows | 139 | 599 | 7 |

The outgoing determinant of order 139 is 345573 modulo 1000003. The determinant of order 957 for the source evaluation augmented by the formal differential is 118669 at the same prime. The formal differential has rational rank 358. Each new nine-letter boundary is entered with zero formal coordinate. The rational dimension of the resulting closed span is at least $957-358=599$. Hence its codimension in the complete classical kernel is at most $745-139-599=7$. Direct source evaluations attain rank 738 of 745. These bounds use classical maps. A quantum pairing requires evaluating the physical first differential on the surviving quotient. [Outgoing minor](jordan_outgoing.json) · [Incoming minor](jordan_incoming.json) · [A8 job 5671196](a8_5671196.json) · [Adjoint outgoing test](adjoint_outgoing.json) · [Adjoint incoming test](adjoint_incoming.json)

From the repository root, using the source plan and ranks above, run

```sh
python3 algorithms/dual_singlets_20261118/f4_pure_target.py --out /tmp/f4-pure.json.gz --workers 4
python3 algorithms/dual_singlets_20261118/f4_pure_chain.py --pure /tmp/f4-pure.json.gz --work /tmp/f4-pure-chain
python3 algorithms/dual_singlets_20261118/f4_joined_boundaries.py --plan /tmp/so9-s5/plan.json.gz --previous /tmp/f4-singlets/ranks.json --work /tmp/f4-jordan --channel jordan --count 450 --workers 48
python3 algorithms/dual_singlets_20261118/f4_joined_incoming.py --plan /tmp/so9-s5/plan.json.gz --previous /tmp/f4-singlets/ranks.json --outgoing /tmp/f4-jordan/jordan_outgoing.json --work /tmp/f4-jordan --channel jordan --count 480 --source-count 300 --workers 48
```

## A necessary and sufficient classical equation for the first differential into $S^5$

At derivative bidegree $(4,4)$ the 666 seven-fermion directions have three separate scalar flavours. The first-order pair rule requires all three flavours at two sites, hence vanishes. In the other 79 directions the two distinguished letters are a scalar and a matter fermion. The remaining eight letters are gauginos. Projecting formal cycles and nine-letter boundaries to this component gives rank at least $436-358=78$ over the rationals. The determinant of order 436 is 710453 modulo 1000003. The full nine-fermion evaluation reaches rank 79. [Projection](pair_support.json) · [A8 job 5671321](a8_5671321.json)

To evaluate the coefficient of $[S^5]$ after acting on an ordinary trace word, insert a colour generator and its trace dual through the two commutators. In pure terminal cohomology the two gauginos anticommute, and each squares to zero. Cutting at the two insertions leaves at most two open strings of gauginos, each of length at most two. A nonvanishing term therefore has at least three spectator bilinears. After removing three, the four-letter classical singlet has dimension $2-1=1$. Two commuting colour generators give top exterior coefficient 4608 on $S^2$. Division by this value extracts the normalized coefficient. The pure fifth power is independently nonzero.

Of the 1,904 ordinary trace singlets, only rows 717, 727 and 1348 have nonzero first-action coefficients, respectively $-18,1,-1/2$. On every row this coefficient equals three times formal classical image coordinate 11 plus one half of each of coordinates 47 and 236. All coordinates are zero-based. All 1,904 residuals vanish over $\mathbb Q$. The first action therefore vanishes on the restrictions of all formal cycles. It vanishes on classical boundaries in cohomology by $Q_1Q_0=-Q_0Q_1$. [Coefficients and rational row identity](first_functional.json) · [A8 job 5671727](a8_5671727.json)

For the three-term cochain in the saved canonical trace order, $Q_1$ is $18 S^5$ and $Q_0$ is $6S^4\operatorname{Tr}_{27}(\phi_1\phi_2\phi_3)$. The first-action functional is consequently nonzero on the 79-dimensional component. The preceding 78-dimensional projected closed span is its entire kernel. Subtracting the corresponding closed representatives from any putative preimage and rescaling reduces existence of a first-order pairing to the single equation

$$
Q_0 X=-6S^4\operatorname{Tr}_{27}(\phi_1\phi_2\phi_3),
\qquad X\text{ in the seven-fermion singlet of dimension }666 .
$$

Conversely a solution added to the saved cochain has zero classical action and first action $18[S^5]$. A8 job 5671407 reached seven-fermion source rank 659 and classical image rank 138. Adjoining the prescribed right-hand side gives rank 139 at the same prime. The unrestricted rational equation is open. [Cochain and integer actions](scalar_source.json) · [Minors](scalar_correction.json) · [Slurm receipt](a8_5671407.json)

```sh
cd algorithms/dual_singlets_20261118
python3 f4_first_source.py --out /tmp/f4-scalar.json
python3 f4_pair_support.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --jordan f4/jordan_incoming.json --out /tmp/f4-pair
python3 f4_scalar_correction.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --work /tmp/f4-correction --count 420
python3 f4_first_functional.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --work /tmp/f4-first
```

A8 job 5671812 also evaluated the ordinary word recipes in the 52-dimensional adjoint representation at the same colour coordinates. Exact reconstruction verifies all $52^2$ brackets in the 27-dimensional representation. The new seven-fermion source, image and augmented ranks are again 659, 138 and 139 modulo 1000003. [Minors](adjoint_traces.json) · [Slurm receipt](a8_5671812.json) · [Replay](../f4_adjoint_traces.py)

A8 job `5672967` also tested 227 seven-fermion parents formed by projecting two adjoint-representation words through the symmetrized Albert-matrix product. On the enlarged family the source rank is 659, the image rank is 138 and adjoining the prescribed image gives rank 139 modulo 1000003. Seven source dimensions remain outside the detected span. [Source minor](adjoint_product_source.json) · [Image and augmented minors](adjoint_product_image.json) · [A8 receipt](a8_5672967.json) · [Algorithm](../f4_adjoint_product.py)

## Restricting the scalar correction to simple factors

For a simple subalgebra, write the scalar trilinear as the invariant trace pairing of the first scalar with the commutator of the other two. The equation for a seven-fermion correction permits brackets of two gauginos and brackets of a gaugino with a scalar. Fix the power of the bilinear and restrict this classical map to the spin and flavour singlet with three scalar flavours. Each source evaluation spans its independently computed character dimension. Each terminal evaluation is injective by the same test. The last two columns test membership in the image over the rationals.

| Group | Bilinear power | Source dimension | Terminal dimension | Image rank | Rank with trilinear times the bilinear power |
|---|---:|---:|---:|---:|---:|
| $G_2$ | 2 | 21 | 11 | 10 | 11 |
| $G_2$ | 3 | 119 | 35 | 35 | 35 |
| $USp(6)$ | 2 | 32 | 18 | 16 | 17 |
| $USp(6)$ | 3 | 254 | 79 | 79 | 79 |
| $SO(7)$ | 2 | 32 | 18 | 16 | 17 |
| $SO(7)$ | 3 | 256 | 81 | 81 | 81 |
| $SU(3)$ | 2 | 36 | 10 | 10 | 10 |

The orthogonal, symplectic and exceptional calculations use the ordered trilinear trace. For orthogonal or symplectic generators it equals half the trilinear commutator trace. The $SU(3)$ calculation keeps both cyclic orientations and uses the commutator trace. The three strict increases have rational separating functionals. The four equalities have rational preimages in 25, 34, 80 and 10 source rows, respectively. A8 job `5672818` computed the first six rows. [Coefficients and bounds](subgroups.json) · [Slurm receipt](a8_5672818.json) · [Unitary coefficients](a2.json) · [Subgroup algorithm](../f4_subgroup_equation.py) · [Unitary algorithm](../f4_a2_equation.py)

In a product of simple subalgebras the restricted trace pairing is a sum of the invariant forms of the factors. Mixed trilinear commutator traces vanish. In $G_2+A_1$ and $C_3+A_1$ the only potentially surviving distribution of four gaugino bilinears is $(3,1)$. In $A_2+A_2$ it is $(2,2)$. This follows from the classical nilpotence exponents $4,4,2,3$ in $G_2,C_3,A_1,A_2$. In $A_1$, imposing the gaugino–scalar commutators makes the antisymmetric colour tensor of each gaugino and each scalar vanish. Multiplying the alternating scalar triple by any gaugino component gives zero by exchanging that component with the first scalar and then antisymmetrizing the gaugino with the second scalar. The table thus supplies a classical preimage after restriction to each of these three subalgebras of $F_4$. The unrestricted seven-fermion equation still needs an independent calculation. [Classical chiral-ring theorem](../../../supplement/supplement.tex) · [Maximal subgroups and the Albert algebra](https://arxiv.org/abs/math/0105155)

## Sources divisible by the gaugino bilinear

Remove each displayed factor of $S$ before computing the classical cohomology of the remaining factor. At level 30 the possible second and third sources have bidegrees $(3,3)$ and $(2,2)$. A8 job `5672668` computes the ordinary trace images. Job `5673042` adds contractions of traceless Jordan multiplication operators in the two adjacent maps. A four-letter open word projected to the Jordan component can be graded-symmetrized. Interchanging two letters produces a Lie bracket. Its remaining cubic projection vanishes because the quadratic projection vanishes and transposition reverses the sign of a symmetrized cubic. The exact Albert-algebra checks cover $26\cdot52$ linear and $26\cdot52^2$ quadratic traces.

| Incoming order | Factors of $S$ removed | Remaining weight | Cochains | Incoming rank at least | Outgoing rank at least | Rational cohomology dimension at most |
|---:|---:|---|---:|---:|---:|---:|
| 2 | 3 | $(0,0,2,2,2)$ | 3 | 0 | 3 | 0 |
| 2 | 2 | $(1,1,2,2,2)$ | 48 | 28 | 20 | 0 |
| 2 | 1 | $(2,2,2,2,2)$ | 625 | 466 | 156 | 3 |
| 3 | 2 | $(0,0,3,3,3)$ | 16 | 6 | 10 | 0 |
| 3 | 1 | $(1,1,3,3,3)$ | 633 | 466 | 156 | 11 |

The three zeroes exclude all classical source classes in the corresponding product images. Before the additional contractions, the two positive upper bounds were 11 and 35. Both improvements follow from nonzero determinants modulo 1000003 and the consecutive classical maps having zero composite. [Ordinary trace calculations](dressed_sources.json) · [Second-order incoming](t2_incoming.json) · [Second-order outgoing](t2_outgoing.json) · [Third-order incoming](t1_incoming.json) · [Third-order outgoing](t1_outgoing.json) · [First A8 receipt](a8_5672668.json) · [Second A8 receipt](a8_5673042.json) · [Classical algorithm](../f4_dressed_jordan.py)

In the last row the Weyl-invariant Cartan singlet has dimensions 1 and 5 at fermion degrees 1 and 3. On four commuting Albert-matrix assignments two products of stress-multiplet components give the minor

$$
\det\begin{pmatrix}0&31850496\\95551488&-47775744\end{pmatrix}
=-3043362286338048.
$$

All 124 classical letter variations vanish on the four assignments. Under the exact current and primary-product input used for the rank-three groups, these two classes have closed completions. This leaves at most nine base directions after quotienting by their span. Multiplication by the protected bilinear gives closed completions of their products. At $(2,2,2,2,2)$ the Cartan source dimension is zero, as the central Weyl element acts by minus one on nine letters. [Cartan characters and evaluations](cartan.json) · [Replay](../f4_cartan_sources.py) · [Current and product input](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

For the derivative-degree and flavour-degree order used above, an order-$r$ differential shifts the five degrees by $(r,r,-r,-r,-r)$ and the letter number by $1-r$. Acting in the order of the archived six spin-to-flavour substitutions takes $S^5$ from $(5,5,0,0,0)$ to $(2,2,2,2,2)$ at length 10. Applying the six reverse substitutions in reverse order gives exactly $7200S^5$ in formal trace polynomials. Each elementary charge has the sign of (167) and (169). The twelve minus signs multiply to one. For $S^6$ the corresponding integer is 25200. [All trace coefficients](charges.json) · [Forward and reverse word actions](../gaugino_primary_order.py) · [Physical charges, Appendix D](https://arxiv.org/html/2610.22332v3#A4)

Consequently the lowered target for $S^5$ has its second-order source at $(0,0,4,4,4)$ and length 11, in the derivative-free BMN sector at $j=24$. A third-order source for this lowered target would have derivative bidegree $(-1,-1)$. For $S^6$ the third-order source of the lowered target is at $(0,0,5,5,5)$ and length 14, at $j=30$. The application to an incoming differential on the original gaugino power uses the perturbative Ward identity for a filtered superconformal action. Such an action induces its classical leading action on each page and commutes with the page differential with the graded sign. If a gaugino power survives to that page, the reverse polynomial identity forces its lowered class to survive there too. An incoming differential into the power would then give one into the lowered class. The kinematic exclusions for the lowered targets are independent of this higher-order Ward input. [Superconformal selection rule and its perturbative status](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex) · [One-loop Ward identity](https://arxiv.org/html/2512.07771#S3.SS4)

For the lowered $j=24$ test, the complete character dimensions at lengths 10, 11 and 12 are 203, 141 and 62. Length ten splits into 11 scalar field-strength and 192 two-fermion singlets. The ordinary trace images give ranks 74 and 51. Albert-algebra contractions raise these to 79 and 57. Hence the complete rational BMN source cohomology satisfies

$$
0\leq \dim H(Q_0)_{(0,0,4,4,4),\,11}\leq 141-79-57=5.
$$

The terminal squared Gram determinant is nonzero on three independent Cartan scalars. This gives the complementary upper bound 61 on the outgoing rank. All seven derivative-free letter identifications intertwine the two saved classical actions. [Both minors, frame seeds and contraction recipes](primary.json) · [A8 job 5673281](a8_5673281.json) · [Replay](../f4_bmn_primary.py)

A8 job `5673532` replayed both minors, then added 1,152 flavour contractions at each source length in both the 27-dimensional Albert and 52-dimensional adjoint traces. The enlarged image ranks are 79 and 57 at the recorded prime. [Determinants and rows](trace_images.json) · [A8 receipt](a8_5673532.json) · [Algorithm](../f4_bmn_trace_images.py)

From the repository root, use the environment of the rank-four calculation.

```sh
cd algorithms/dual_singlets_20261118
python3 gaugino_primary_order.py --out /tmp/charges.json
python3 f4_bmn_primary.py --work /tmp/f4-bmn-primary --tasks 640 --workers 16
python3 f4_subgroup_equation.py --group G2 --power 3 --plan mixed33/outgoing.json.gz --out /tmp/g2-three.json
python3 f4_a2_equation.py --out /tmp/a2-two.json
python3 f4_cartan_sources.py --formal mixed11/classical.json.gz --out /tmp/f4-cartan.json
```

## Two more classical boundaries at level 24

Allow a spectator trace alongside two words contracted through traceless Jordan multiplication. At length eleven seeds 3712362 and 3712377 add two independent images. The outgoing minor of order 59 is 446538 modulo 1000003. The length-ten minor of order 79 is 579125 at the same prime. This bounds the rational second-order source dimension by $141-79-59=3$. Two Jordan projections in a trace, two adjoint bridges between traces, and cyclic products of three to five adjoint-projected Albert words reproduce these ranks in the tested families. [Both minors and recipes](joined_images.json) · [A8 job 5673641](a8_5673641.json) · [Two Jordan projections](double_jordan.json) · [Two adjoint bridges](adjoint_bridges.json) · [Adjoint-projected words](adjoint_compression.json) · [Projection algorithm](../f4_bmn_adjoint_compression.py)

## The Cartan restriction of the level-30 third-order source

The lowered sixth power has derivative-free incoming third-order weight $(0,0,5,5,5)$ and length 14. A Cartan cochain has one matter fermion $\psi_i$ and thirteen scalars. Its Weyl- and flavour-invariant space has dimension two. Write $h^i$ for the three Cartan scalars, $v$ for the oriented normal to their ordered triple in the four-dimensional Euclidean root realization, and $p_6(v)$ for the sum of $(u\cdot v)^6$ over the 24 short roots $u$. Two candidate coefficients of $\psi_i$ are $(v\cdot v)^2h^i$ and the Hessian of $p_6$ at $v$ applied to $h^i$. Both transform as the scalar flavour vector. Under a Weyl reflection $v$ transforms with an additional orientation sign, which cancels in these even fourth-degree expressions.

At three Cartan elements on a root hyperplane, vary the scalars by $sH,sE,sF$ in its root $\mathfrak{sl}_2$, using the usual Cartan, raising and lowering matrices. The scalar commutators start at order $s^2$. For any invariant cochain linear in a matter fermion, the differential at order $s^3$ pairs each derivative of its fermion coefficient with the commutator of the other two variations. The projection of this derivative from the root $\mathfrak{sl}_2$ to itself is equivariant, hence scalar. Its scalar is determined by its derivative in the Cartan root direction. The three cyclic diagonal contractions all give $\operatorname{Tr}(H[E,F])=2$ and the six off-diagonal contractions vanish. A closed cochain therefore obeys one trace equation at each root hyperplane.

Using normals $(0,0,0,1)$ and $(0,0,1,1)$ and multiplying each trace equation by its root norm squared gives the following two rows on the displayed polynomials,

$$
\begin{pmatrix}3&225/2\\24&-360\end{pmatrix},
\qquad \det=-3780.
$$

The rows establish the independence of the two Cartan polynomials and the vanishing of the full classically closed Cartan image. A third-order source in this block must vanish on the Cartan. Its full cohomology dimension is open. [Rational calculations and root probes](root_restriction.json) · [Replay](../f4_bmn_root_restriction.py)
