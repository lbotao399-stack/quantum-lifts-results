# Perturbative survival of the sixth F₄ glueball

**Proposition.** In the relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM with gauge algebra $\mathfrak f_4$, the pure glueball $[S^6]$ at $j=36$ is nonzero on $E_\infty$ of the perturbative supercharge spectral sequence. In particular the entire incoming $d_2$ has rank zero at this target. Here $S$ uses the native trace convention of the first-page certificate. The comparison with untwisted physical operators requires the hypothesis of Appendix B. [First-page Ward calculation](../eight_albert_review/README.md) · [Manuscript](https://arxiv.org/abs/2611.24160)

## A product completes the classical kernel

The sole Taylor diagram for the second source is $(4,4,3,3)$. The source cochain dimension is 3720, and the outgoing classical rank is 703. The full $Q_0$-kernel therefore has dimension 3017. Previously certified length-twelve predecessors supply 3016 independent boundaries in this kernel. [Classical source and third-order ranks](README.md#incoming-differential-at-the-sixth-power)

In its ordinary four-scalar polarization the required highest weight is $(4,4,4,2)$. Use the previously constructed [nine-scalar Albert polynomial](../../../bmn_primary/f4_scalar/README.md), the scalar bilinear and the one-mark bilinear. The latter two factors are represented by $\operatorname{Tr}_{26}(X_iX_j)$ and $\operatorname{Tr}_{26}(Y_{ij}X_k)$ in the native colour matrices. The first jets $X_i$ have $Q_0X_i=0$, and the alternating second jets have $Q_0Y_{ij}=[X_i,X_j]$. Hence the classical derivative of the marked factor is the alternating three-form $\operatorname{Tr}_{26}([X_i,X_j]X_k)$. The indices $i,j,k$ label the four auxiliary scalar flavours.

Polarize the nine-scalar factor in three triples. Each triple alternates three of the four flavours and carries the sign of its omitted flavour as in the supplied construction. Its three triple factors are symmetrized. The resulting twenty coefficients, the ten scalar bilinears and the twenty-four marked bilinears give 96 formal products **at the specified weight**. Impose the three simple raising equations and the vanishing of the alternating three-form for each fixed nine-scalar and scalar-bilinear factor. Rational elimination leaves one vector with 87 nonzero entries. The equations impose closure separately on the marked factor for each fixed pair of spectators. Consequently the vector is a sum of products of positive-degree classical cycles. The [product theorem at level 36](../first_page/glueball_factors.md) gives it zero $d_2$ image at $[S^6]$. The construction also checks the finite flavour shears and an independent polarization of the cube. [Rational coefficients and native evaluator](albert_product_cycles.py) · [Controls and fixed-minor reconstruction](sixth_product_minor.py)

Adjoin this product as column 49 (zero-based) to the 49-column extension of the original order-3002 boundary minor. Keep fourteen boundary columns and the product, using the same ordered rows at both primes.

| Prime | Order-3002 minor | Order-14 complement on boundaries | Order-15 complement including the product | Order-3017 minor |
|---:|---:|---:|---:|---:|
| 1009 | 833 | 587 | 737 | 449 |
| 1013 | 615 | 953 | 346 | 60 |

The two full residues are $833\cdot737\bmod1009=449$ and $615\cdot346\bmod1013=60$. All columns and point assignments have a common rational definition. Either nonzero minor proves independence over $\mathbb Q$. Thus 3016 boundaries and one product span the entire classical kernel, and every classical cycle has zero second image at the target. This argument allows the product itself to be an additional boundary. [First minor](sixth-product-1009.json) · [Same minor at the second prime](sixth-product-1013.json) · [A8 accounting](sixth-product-jobs.json)

The native Ward covector takes the value $-15313094596800$ on $S^6$ and proves first-page survival. At order three the unique source has classical cohomology of dimension $1862-1434-428=0$. The Taylor bound $r\le 6-3$ exhausts incoming orders. The pure target is strictly closed. These facts prove the proposition. [Ward covector](../first_page/README.md) · [Order-three certificates and Taylor selection](README.md#incoming-differential-at-the-sixth-power)

A8 job 5719653 independently regenerated the fixed order-3017 minor modulo 1013 from the native evaluation points. The exact rank and determinant agree with the saved certificate. [Independent rebuild and scheduler readback](triple_exchange/accounting.json)

## Seventh power

Both $[S^6]$ and $[S^7]$ are nonzero on $E_2$. The [seventh-power product theorem](seventh_products.md) now excludes the entire subspace of positive-degree classical products from the incoming $d_2$ at level 42. On each later page a product of positive-degree classes on that page also gives zero at the pure target. A contributing pure spectator would be $S^k$ with $k\ge1$, and multiplying an incoming image $S^{7-k}$ by $S^{k-1}$ would kill the permanent $[S^6]$. A seventh-power hit, if present, first occurs at order two, three or four. The independent primitive source spaces at level 42 still require a test. [Independent classical test at level 42](seventh_exchange_kernel.md)

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -B sixth_product_minor.py --prime 1009 --controls
python3 -B sixth_product_minor.py --prime 1013 --controls
python3 -B sixth_product_minor.py --prime 1013 --workers 40
~~~
