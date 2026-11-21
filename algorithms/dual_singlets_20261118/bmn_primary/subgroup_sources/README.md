# Classical restrictions of the two $G_2$ candidates

An odd BMN generator accounting for the level-36 primary index has flavour partition $(6,6,4)$ at level 32 or $(6,6,6)$ at level 36. The remaining possibilities for its letter length are respectively $11,13$ and $11,13,15$. The classical differential $Q_0$ is a Lie-bracket differential, so restriction to a Lie subalgebra defines a map on classical cohomology. [Source index and full-group edge ranks](../README.md#two-remaining-g2-weights)

| Flavour | Length | Cochains on $SU(2)\times SU(2)$ | Cohomology there | Principal $SU(2)$ cochains | Cohomology there |
|---|---:|---:|---:|---:|---:|
| $(6,6,4)$ | 11 | 4348 | 0 | 39 | 0 |
| | 13 | 1706 | 0 | 9 | 0 |
| $(6,6,6)$ | 11 | 2946 | 0 | 12 | 0 |
| | 13 | 3218 | 0 | 12 | 0 |
| | 15 | 752 | 0 | 3 | 0 |

For the two simple factors, the cochain algebra is their tensor product. Taking invariants and applying the Künneth formula over the rationals gives the cohomology as the tensor product of their cohomologies. The two letter lengths add. For each pair of highest weights the Littlewood–Richardson coefficient gives its multiplicity in the indicated flavour partition. Every contributing rank-one cochain dimension was checked by the $SU(2)$ adjoint Weyl character. For a permutation cycle of length $k$ the adjoint weights give $z^{-k}+1+z^k$, with $z$ recording multiples of the positive root. Taking the coefficient at zero minus that at one after multiplication gives the trace on colour invariants. The published modular ranks bound rational cohomology from above. All five upper bounds in each restriction are zero. [Rank extraction and replay](../g2_subgroup_sources.py) · [Compact rank table](upstream_ranks.json) · [Fortuity rank-one data](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su2/bmn_hw_cohomology_j44.json)

The long roots form an $A_2$ subsystem. The Weyl group of $G_2$ contains inversion on the Cartan. Inversion normalizes $A_2$ and acts on $SU(3)$ by charge conjugation up to an inner automorphism. A $G_2$-invariant therefore restricts to the even charge-conjugation space. At $(6,6,4)$ the upper bounds before imposing colour parity are

$$
11649-5413-6235=1,\qquad 6366-4434-1932=0
$$

at lengths eleven and thirteen. At length eleven a product of four closed graviton traces has a highest-weight Cartan coefficient $18/5$. Every classical boundary vanishes on this Cartan evaluation. The four traces have lengths $2,3,3,3$ and charge-conjugation eigenvalues $+1,-1,-1,-1$. The product spans the entire rational line, so its even subspace is zero. [Factor inputs and rational minor](su3_664_product11.json) · [Fortuity rank-three defining-matrix data](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/hw_cohomology_j32_records/su3_c_6_6_4.jsonl)

At $(6,6,6)$ and length fifteen, direct projection of every trace evaluation by $1+C$, with $C$ charge conjugation, gives

$$
1918-1471-447=0 .
$$

The adjacent even cochain dimensions are $4365,1918,560$. At prime $2147483647$ the two pivot minors have determinants $1118222283$ and $776053499$. Both matrices come from integral trace evaluations after clearing the highest-weight denominators. Their ranks therefore bound the rational ranks from below, and $Q_0^2=0$ makes their sum exact. The same engine returns ranks $(93,38)$ in the 131-dimensional $(4,4,4)$ length-nine control and $(87,93)$ in its 181-dimensional length-eight neighbour. [A6 job 7658790](../7658790/even_15.json) · [Receipt](../7658790/receipt.json) · [Independent length-nine calculation](su3_444_independent.json)

The maximal proper connected compact subgroups of $G_2$ have Lie algebras $A_2$, $A_1+A_1$, and the principal $A_1$. Thus the restrictions at $(6,6,4)$, lengths eleven and thirteen, and $(6,6,6)$, length fifteen, are classical boundaries on **every proper connected compact subgroup**. At $(6,6,6)$, lengths eleven and thirteen, the long-root $SU(3)$ remains to be determined. The unrestricted $G_2$ cohomology and its quantum differential require separate calculations. [Thomas, Theorem 3.1 and Table 1](https://arxiv.org/pdf/1501.04858) · [Numerical bounds](bounds.json)

There is also a fixed-scalar consequence at **all five** candidate lengths. Fix any scalar linear combination to a value with nonzero semisimple part. Its semisimple centralizer has Lie algebra a Cartan or a root $A_1$ with a one-dimensional centre. Every $G_2$ root has a perpendicular root of the other length, so each such Levi lies inside a conjugate of $A_1+A_1$. The restricted source is already a $Q_0$ boundary there. On the nonzero eigenspaces of the semisimple part, the adjoint action of the fixed scalar is invertible. The linear part of $Q_0$ contracts the transverse pairs. Each nonlinear term raises the number of unfixed letters, which is bounded by the remaining total flavour degree. Iterating the linear contraction therefore terminates. Hence a possible class is exact on every such fibre, including for arbitrary other colour entries. This statement allows the primitive to depend on the fixed value. [Fixed-scalar contraction](../../f4/levi/README.md#fixing-a-scalar)

From the repository root, the ordinary and parity-resolved calculations can be replayed with Python, NumPy, SciPy, python-flint, and C++17.

```sh
d=algorithms/dual_singlets_20261118
p="$d/bmn_primary"
python3 "$p/g2_subgroup_sources.py" --out /tmp/g2-subgroup-bounds.json
cmp "$p/subgroup_sources/bounds.json" /tmp/g2-subgroup-bounds.json
c++ -O3 -std=c++17 -pthread -I"$d/taylor_hooks" -I"$d/taylor_hooks/native" -I"$d/b2_level36" "$p/su3_even_bmn.cpp" -o /tmp/su3-even-bmn
/tmp/su3-even-bmn 6 6 6 15 847151 2147483647 8 "$p/subgroup_sources/su3_even_666.txt" > /tmp/su3-even-666.json
cmp "$p/7658790/even_15.json" /tmp/su3-even-666.json
```

Passing `--upstream /path/to/fortuity-ai/data` checks both Git blobs at the pinned revision and re-extracts the rank table.
