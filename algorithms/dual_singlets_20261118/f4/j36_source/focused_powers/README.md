# The sixth and seventh glueball powers of F₄

The sixth power of the glueball $S$ is a perturbative permanent cycle of the relative holomorphic jet complex. At $j=36$ its unique potential second-differential source has a classical kernel of dimension 3017, spanned by 3016 classical boundaries and one product with zero obstruction. The fixed order-3017 minor has residues 449 and 60 at primes 1009 and 1013. The seventh power at $j=42$ survives the first differential. A possible first hit at order two, three or four requires a source outside the corresponding product subspace. [Sixth-power proof and replays](sixth_survival.md) · [Seventh-power exchange identities](seventh_exchange_kernel.md)

## Source ranks

At level 42 a second-order hit requires a source combining matter-fermion and differentiated-scalar terms. In the complete family of 518 marked Albert words with at most three exchanges the classical split has rank 191 and all 327 closed combinations give zero classes. The global mixed Ward equation remains open. [Reduction on the full source](seventh_ward.md) · [Native rational kernel](triple_exchange/README.md)

For each potentially incoming order $r\ge2$, the Taylor diagram is $(p-r,p-r,3^r)$ at target $S^p$. Let $n$ be the size of the diagram. The cochains in the table have $n-2,n-1,n$ letters, respectively. The putative source is in the middle. [Selection and grading](../first_page/README.md#the-fifth-glueball-survives-to-perturbative-infinity) · [Exact rational character calculation](source_sizes.py)

| Power | Order | Predecessors | Source cochains | Classical target cochains | Upper bound on classical source cohomology |
|---|---:|---:|---:|---:|---:|
| $S^6$ | 2 | 10007 | 3720 | 704 | **1** |
| $S^6$ | 3 | 4295 | 1862 | 428 | **0** |
| $S^7$ | 2 | 217297 | 62866 | 8553 | |
| $S^7$ | 3 | 333841 | 101602 | 14062 | |
| $S^7$ | 4 | 88602 | 29305 | 4707 | |

For the third-order $S^6$ source, use the three-scalar polarization of shape $(5,5,5)$. Its first, second and third Taylor jets obey $Q_0X_i=0$, $Q_0Y_{ij}=[X_i,X_j]$ and $Q_0Z_{ijk}=[X_i,Y_{jk}]+[X_j,Y_{ki}]+[X_k,Y_{ij}]$. Indices $i,j,k$ label the three scalar flavours. Alternating five columns of three labels projects to the required determinant character. Gauge invariance follows by using defining traces, the native adjoint and Albert colour frames, and the inverse of the trace form between open legs.

| Classical map | First minor size | Determinants mod 1009 and 1013 | Extended minor size | Determinants mod 1009 and 1013 |
|---|---:|---:|---:|---:|
| Into the source | 1396 | 680, 780 | 1428 | 191, 544 |
| Out of the source | 419 | 200, 723 | 426 | 101, 944 |

The extension uses colour traces with two separated open insertions. Koszul and Jacobi boundaries vanish on the substitution $Y_{ij}=[X_i,X_j]$. The selected square minors are regenerated independently at the second prime. These trace minors alone give $\dim H(Q_0)\le 1862-1428-426=8$. The Albert minors below exhaust the classical kernel. [Trace contractions](three_axis.py) · [Separated insertions](sandwich.py) · [A8 receipts](cluster.json) · [Previously certified bound of 15 at order two](../global_boundaries/README.md)

## Incoming differential at the sixth power

In the relative full-jet complex, the only Taylor diagram for an incoming $d_3$ on $S^6$ has 1,862 source cochains. Antisymmetrize five triples of scalar labels. For the outgoing map split one second jet. For the incoming map split two second jets by the graded product rule or split one third jet by the Jacobi formula above. The resulting gauge-invariant contractions use four open lines of the native 26 representation and the Albert cubic. The two-second-jet difference and the Jacobi insertion each vanish when all remaining marks are replaced by their commutators. A8 job 5718997 regenerated the indicated rows at both primes. [Native contractions and replay](third_albert.py)

| Classical map | Minor size | Determinant modulo 1009 | Determinant modulo 1013 |
|---|---:|---:|---:|
| Into the third-order source | 1434 | 767 | 849 |
| Out of the third-order source | 428 | 862 | 258 |

A nonzero reduction gives a rational lower rank bound. Since $Q_0^2=0$ and $1434+428=1862$, both ranks are exact. The classical cohomology of this source is zero. Hence an incoming $d_3$ on $S^6$ is impossible. The Taylor selection rule $r\le p-3$ for $S^p$ excludes all $r\ge4$ here. [Incoming pivots](third-search-2.json) · [Outgoing pivots](third-search-1.json) · [Second-prime incoming](third-replay-2.json) · [Second-prime outgoing](third-replay-1.json)

For the sole possible incoming $d_2$ on $S^6$, the outgoing classical rank is 703. An independently checkable fourteen-column extension of the order-3002 incoming minor has determinant $587\bmod1009$. The extended determinants on the same ordered rows and columns are $833\cdot587=615\bmod1009$ and $615\cdot953=581\bmod1013$. Consequently the rational source cohomology has dimension at most $3720-703-3016=1$. Both recorded crossings have exactly zero corrected quantum pairing. The [exceptional product minor](sixth_survival.md) accounts for the remaining possible direction and excludes the entire incoming $d_2$. [New classical boundary columns](albert_boundaries.py) · [First prime](boundary-1009.json) · [Second prime](boundary-1013.json) · [A8 accounting](sixth-source-jobs.json) · [Rebuild from the previous native rows](sixth_boundary_minor.py)

No positive-degree classical product can hit $S^7$ at order two. At orders three and four, a contributing product on that page would have a surviving pure spectator $S^k$, $k\ge1$. Multiplying the required image $S^{7-k}$ by $S^{k-1}$ would kill the permanent $S^6$. [Seventh-power product proof](seventh_products.md)

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -B third_albert.py replay --order 1 --prime 1013 --workers 4 --certificate third-search-1.json
python3 -B third_albert.py replay --order 2 --prime 1013 --workers 4 --certificate third-search-2.json
python3 -B sixth_boundary_minor.py --prime 1013 --workers 4 --fixed boundary-1009.json
~~~

## Which sources can be seen on commuting scalars

In the ordinary scalar polarization, removing the two slots of the alternating mark $Y_{ij}$ removes one box from each of two different rows. The Weyl character average over all 1152 elements gives the following dimensions for scalar polynomials linear in a marked Cartan colour. [Weyl and Levi calculation](support.py) · [Exact results](support.json)

| Power and order | Ordinary shape | Cartan cochains | Cartan restriction of classical cycles |
|---|---|---:|---|
| $S^6$, 2 | $(4,4,4,2)$ | 0 | zero |
| $S^6$, 3 | $(5,5,5)$ | 2 | **zero by the root equations below** |
| $S^7$, 2 | $(4,4,4,2,2)$ | 0 | zero |
| $S^7$, 3 | $(5,5,5,2)$ | 1 | dimension at most one |
| $S^7$, 4 | $(6,6,6)$ | 0 | zero |

For the third source of $S^6$, take a Cartan triple $X_1,X_2,X_3$. In the native four Cartan coordinates, let $N$ be specified by $B(N,v)=\det(X_1,X_2,X_3,v)$, with $B=\operatorname{Tr}_{26}/6$. Set $(Y^1,Y^2,Y^3)=(Y_{23},Y_{31},Y_{12})$ for the three cyclic components. This is the test polynomial convention used in the replay. Two Weyl-invariant test polynomials are

$$
 B(N,N)^2\sum_{i\in\{1,2,3\}\;\mathrm{flavours}} B(X_i,Y^i),
 \qquad
 \sum_{i\in\{1,2,3\}\;\mathrm{flavours}}
 \operatorname{Tr}_{26}(X_iY^iN^4).
$$

Both transform by the fifth power of the flavour determinant. To test a classical cycle, choose three elements spanning the Cartan hyperplane of a root $\alpha$ and perturb them by $t u_i$ in its root $\mathfrak{sl}_2$. The background commutes with this $\mathfrak{sl}_2$, so the three commutators of the perturbed scalars are $t^2[u_j,u_k]$. A derivative of the marked-colour coefficient in an active $\mathfrak{sl}_2$ direction, paired with an active mark, is an invariant bilinear form on the adjoint of $\mathfrak{sl}_2$. It is fixed by evaluating both arguments on the Cartan vector dual to $\alpha$. The scalar coefficient of $t^3$ in $Q_0$ is therefore the cyclic sum of these three derivatives times a common nonzero normalization of $B(u_1,[u_2,u_3])$. Classical closure forces this sum to vanish for every root. With the long and short native roots 27 and 10, the rows on the two displayed polynomials are

$$
\begin{pmatrix}24&-12\\3&15/4\end{pmatrix},
\qquad
\det=126 .
$$

These two rows establish independence, exhaust the character dimension two, and exclude both Cartan restrictions of a cycle. The explicit derivative and exact degree-five interpolation agree on all four simple roots. [Rational calculation](wall_test.py)

The same four zero entries in the table vanish on **every reduced commuting tuple**, with arbitrary mark. Indeed, put all semisimple parts in the centre of their common Levi algebra and the commuting nilpotent parts in its positive nilradical. The centre torus projects the mark to the Levi. If its semisimple rank is two, three or four, the minimum total number of positive-root entries in every allowed scalar diagram exceeds the maximum root height $3,5,11$, respectively. One mark cannot balance this weight. For semisimple rank one the commuting nilpotents are multiples of one element and are limits of simultaneous semisimple tuples. On the Cartan the torus projects the mark to the Cartan itself. The chart homotopy for a regular linear combination of scalars therefore provides a local classical primitive in these four cases. For $S^6$ at order three the denominator-clearing exponent is 10, for $S^7$ at order two it is 12 or 14, and at order four it is 12. These restrictions leave open the cohomology on singular configurations. [Explicit homotopy](../first_page/regular_source.md) · [All shape and height bounds](support.json)

On $E_2=H(E_1,d_1)$ the marked product is compatible with $d_2$. A source hitting $S^6$ at this order would give a source hitting $S^7$ by multiplication by $S$. The sixth-power global kernel test excludes this mechanism. The seventh-power independent source remains open in the relative complex. Comparison with untwisted physical operators uses Appendix B of the manuscript.

~~~sh
python3 -B source_sizes.py
python3 -B support.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B three_axis.py replay --order 1 --prime 1013 --workers 4 --certificate three-1-1009-429.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B three_axis.py replay --order 2 --prime 1013 --workers 4 --certificate three-2-1009-1664.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B sandwich.py replay --stage 3 --order 1 --prime 1013 --workers 4 --certificate sandwich-3-1.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B sandwich.py replay --stage 3 --order 2 --prime 1013 --workers 4 --certificate sandwich-3-2.json
~~~

## The first correction for the seventh power

In the intermediate spin and flavour singlet at dotted and flavour degree $(6,6\mid1,1,1)$, the four terminal Taylor diagrams $(9,6),(8,6,1),(7,6,2),(6,6,3)$ have cochain dimensions $112,573,1167,1107$. A covector killing the same-species brackets lies in the three harmonic tensor spaces of dimensions $43,16,9$. On exactly these 68 columns the homogeneous differentiated-scalar bracket matrix has determinant $909\bmod1009$ and $962\bmod1013$. Hence its rational kernel is zero. The image of the complete classical map fills all $2959$ intermediate terminal cochains. Every classical cycle in the $62866$-dimensional second source of $S^7$ therefore admits a first correction. [Native matrix and character replay](seventh_intermediate.py) · [Three modular reductions](seventh-intermediate.json)

For either power, write $\lambda$ for its classical target covector and $\eta$ for the rational solution of $\lambda Q_1=\eta Q_0$. The bare $Q_2$ vanishes on the selected sources by flavour counting. For $Q_0X=0$, choose $Q_0X_1=-Q_1X$. Then

$$
\lambda(d_2[X])=\lambda(Q_1X_1)=-\eta Q_1X.
$$

An incoming second differential is equivalent to the row $-\eta Q_1$ increasing the rank of $Q_0$ by one on the selected source. For $S^7$ this classical matrix has 8553 rows and 62866 columns. The preceding rank 2959 pertains to the intermediate equation for $X_1$.

For a level-36 source $X$, multiplicativity gives $d_2([X][S])=d_2[X]\,[S]$. The ratio of the seventh and sixth target covectors on their native powers is $194893931232000/15313094596800=140/11$. Thus the two numerical obstruction pairings on a source and its product with $S$ have this fixed ratio. A seventh-power source outside this product image can still hit the target independently. The comparison here is in the relative jet complex.

## Both orders of Dot's crossed graph

Both orders of the two specified Albert exchanges now have zero corrected pairing over the rationals. For order A, 18 left-mark and 30 right-mark graph entries each give a total of zero after the exterior-order conversion. The classical section is zero independently. For order B, both the ordered graph and its section give $-299015516814300/49$. [Both calculations, coordinates and replay](rational_section/README.md)

## Further second-differential tests

For $S^6$, nine ordered exchange choices and nineteen combinations of three scalar insertions give zero after classical subtraction over the rationals. The two choices of first marked endpoint agree in every case. Six of the nine unsubtracted scores and eight of the nineteen are nonzero. The computation evaluates 312 distinct pairs of extremal graph contractions in the first family and 670 in the second. These counts refer to source recipes within the specified families. Their independence in classical cohomology is undetermined. [Exact subtractions and reconstruction](rational_section/README.md#varying-the-exchanges-and-scalar-insertions)

For $S^7$, consider the defining-trace presentation with one matter fermion, four scalars and five gauginos of each dotted species, imposing the inactive classical splittings and the polarized quartic trace identities. Direct enumeration gives 162 reduced monomials. They contain respectively three, four and five factors of the defining glueball trace with multiplicities 33, 99 and 30. On the 19 spin and flavour singlets the classical split has rank 14 before the quartic identities and rank 5 afterwards. Multiplying the fourteen previously constructed level-36 cycles by the defining glueball trace gives a level-42 column space of rank 14. Appending the independently calculated level-42 kernel keeps its rank 14. Each of these products has zero second image by the product theorem and the exact level-36 pairings. Further finite-colour identities and sources outside this trace presentation still enter the full level-42 test. [Calculation](seventh_trace.py) · [Ranks](seventh-trace.json) · [Level-36 cycles and quantum scores](../first_page/trace_obstruction.md)

For the full selected level-42 source, a positive-degree classical product can hit $[S^7]$ under $d_2$ exactly when $[S^6]$ is hit at level 36. The degree-six terminal quotient has 16 independent tensors. The two nonempty first targets of a bilinear complement have dimensions 238 and 209 and are classical boundaries. These checks allow the Leibniz rule to be applied to all classical products in this diagram. An independent seventh-power hit would require a class outside their span. [Proof, minors and replay](seventh_products.md)
