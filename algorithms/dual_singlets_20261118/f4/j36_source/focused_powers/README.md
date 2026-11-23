# The sixth and seventh glueball powers of F₄

Both pure powers of the glueball $S$ survive the first supercharge differential $d_1$ in the relative holomorphic jet complex. Their target covectors take the respective nonzero values $-15313094596800$ and $-194893931232000$. The withdrawn crossed eight-Albert source for $S^6$ has **exactly zero** corrected $d_2$ pairing over the rationals. [First-page Ward certificates](../eight_albert_review/README.md) · [Exact crossed-source subtraction](rational_section/README.md)

## Source ranks

For each potentially incoming order $r\ge2$, the Taylor diagram is $(p-r,p-r,3^r)$ at target $S^p$. Let $n$ be the size of the diagram. The cochains in the table have $n-2,n-1,n$ letters, respectively. The putative source is in the middle. [Selection and grading](../first_page/README.md#the-fifth-glueball-survives-to-perturbative-infinity) · [Exact rational character calculation](source_sizes.py)

| Power | Order | Predecessors | Source cochains | Classical target cochains | Upper bound on classical source cohomology |
|---|---:|---:|---:|---:|---:|
| $S^6$ | 2 | 10007 | 3720 | 704 | 15 |
| $S^6$ | 3 | 4295 | 1862 | 428 | **8** |
| $S^7$ | 2 | 217297 | 62866 | 8553 | |
| $S^7$ | 3 | 333841 | 101602 | 14062 | |
| $S^7$ | 4 | 88602 | 29305 | 4707 | |

For the third-order $S^6$ source, use the three-scalar polarization of shape $(5,5,5)$. Its first, second and third Taylor jets obey $Q_0X_i=0$, $Q_0Y_{ij}=[X_i,X_j]$ and $Q_0Z_{ijk}=[X_i,Y_{jk}]+[X_j,Y_{ki}]+[X_k,Y_{ij}]$. Indices $i,j,k$ label the three scalar flavours. Alternating five columns of three labels projects to the required determinant character. Gauge invariance follows by using defining traces, the native adjoint and Albert colour frames, and the inverse of the trace form between open legs.

| Classical map | First minor size | Determinants mod 1009 and 1013 | Extended minor size | Determinants mod 1009 and 1013 |
|---|---:|---:|---:|---:|
| Into the source | 1396 | 680, 780 | 1428 | 191, 544 |
| Out of the source | 419 | 200, 723 | 426 | 101, 944 |

The extension uses colour traces with two separated open insertions. Koszul and Jacobi boundaries vanish on the substitution $Y_{ij}=[X_i,X_j]$. The selected square minors are regenerated independently at the second prime. Nilpotence and the two lower bounds give $\dim H(Q_0)\le 1862-1428-426=8$. The third differential still requires earlier-page survival and the corrected quantum pairing. [Trace contractions](three_axis.py) · [Separated insertions](sandwich.py) · [A8 receipts](cluster.json) · [Previously certified bound of 15 at order two](../global_boundaries/README.md)

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

On $E_2=H(E_1,d_1)$ the marked product is compatible with $d_2$. A source hitting $S^6$ at this order would give a source hitting $S^7$ by multiplication by $S$. Each power remains an open higher-order target in the relative complex. Comparison with untwisted physical operators uses Appendix B of the manuscript.

~~~sh
python3 -B source_sizes.py
python3 -B support.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B three_axis.py replay --order 1 --prime 1013 --workers 4 --certificate three-1-1009-429.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B three_axis.py replay --order 2 --prime 1013 --workers 4 --certificate three-2-1009-1664.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B sandwich.py replay --stage 3 --order 1 --prime 1013 --workers 4 --certificate sandwich-3-1.json
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B sandwich.py replay --stage 3 --order 2 --prime 1013 --workers 4 --certificate sandwich-3-2.json
~~~
