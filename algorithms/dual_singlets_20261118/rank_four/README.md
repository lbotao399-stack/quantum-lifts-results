# The fifth gaugino power in $SO(9)$ and $USp(8)$

At level $j=2(E+J_L)=30$, put $S=\operatorname{Tr}(\lambda_+\lambda_-)$ with the defining trace and the two gaugino components of the source convention. The terminal singlet $S^5$ has ten letters and doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(0,0,10,10,10)$. Its first-order source has charges $(1,1,10,10,10)$.

| Terminal singlet | Formal dimension | Finite-colour dimension | Classical boundary rank | Rank with $S^5$ |
|---|---:|---:|---:|---:|
| $SO(9)$ | 20 | 20 | 19 | 20 |
| $USp(8)$ | 20 | 19 | 19 | 19 |

The formal boundary rank is 19 over the rationals. Each nine-letter classical primitive contains at most four traces. Replacing one letter by its classical variation preserves trace number, so the coefficient of the product of five bilinears vanishes on its image. The finite-colour cochain dimensions follow from independent singlet characters. The integral evaluations attain the indicated dimensions and ranks. Thus extraction of the canonical five-bilinear coefficient, divided by its sign in $S^5$, defines a rational functional on the entire $SO(9)$ target. It evaluates to one on $S^5$. [Target minors](pure.json)

At derivative bidegree $(4,4)$ and flavour degrees $(1,1,1)$, the ten-letter source has 1,904 formal singlets. Its two finite-colour fermion degrees contribute $1707+197=1904$. The classical supercharge $Q_0$ maps them to a singlet space of dimension 357 at length eleven. One first-order contraction changes the number of traces by at most one. Of 52,608 words appearing in the source basis, 1,576 have at least four traces. Exact pair contraction on these words gives only three nonzero singlet coordinates for the normalized fifth-power coefficient of $Q_1$. Their zero-based row numbers and coefficients are $(717,-84),(727,12),(1348,-6)$. [Recoupling rule](../../su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-protected/q1d_exact.py)

| Prime | Classical minor order | Determinant | Minor order after adjoining the $Q_1$ coefficient | Determinant |
|---:|---:|---:|---:|---:|
| 1000003 | 357 | 739063 | 358 | 973909 |
| 1000033 | 357 | 375392 | 358 | 637559 |

Each modular minor lifts to a nonzero rational minor by lifting the independent entries of the antisymmetric matrices. The classical rank upper bound is 357. The augmented map therefore has rational rank 358, and some rational vector in the classical kernel has a nonzero fifth-power coefficient under $Q_1$. The identity $Q_0Q_1+Q_1Q_0=0$ makes this coefficient vanish on classical boundaries in the source. The one-dimensional target cohomology then gives

$$
\operatorname{im}d_1=\operatorname{span}_{\mathbb Q}\{[S^5]\}\quad\text{in this terminal }SO(9)\text{ singlet},\qquad E_2=0.
$$

The source's independence of previously counted products and descendants is open. [All pivot indices, determinants and quantum coordinates](first.json) · [A8 job 5670470](a8_5670470.json)

From the repository root, with the Python requirements of the independent-prime runtime and a C++17 compiler, run

```sh
python3 algorithms/dual_singlets_20261118/rank_four_s5.py --work /tmp/so9-s5 --workers 24
cmp /tmp/so9-s5/pure.json algorithms/dual_singlets_20261118/rank_four/pure.json
cmp /tmp/so9-s5/first.json algorithms/dual_singlets_20261118/rank_four/first.json
```
