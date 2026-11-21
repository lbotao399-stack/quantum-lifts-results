# Instanton partners of the two fermionic $SU(3)$ seeds

The $F_1$ seeds at levels 30 and 32 have flavour partitions $(7,4,4)$ and $(8,4,4)$, letter lengths 10 and 11, and colour charge conjugation $-1$ and $+1$. Both have doubled spins $(2J_1,2J_2)=(5,5)$. An incoming supercharge partner has doubled spins $(6,6)$, an outgoing one $(4,4)$. The three doubled index charges are $(14,8,8)$ and $(16,8,8)$ respectively. [Public seeds](https://arxiv.org/html/2412.08695v2#S2) · [First-order protection](https://arxiv.org/html/2610.20340v2#S4.SS2)

Let $t$ be the equal number of Taylor derivatives in each spacetime direction. At level 32 the three flavour degrees are $(8-t,4-t,4-t)$. The incoming and outgoing lengths are $10-t$ and $12-t$. At level 30 replace the first flavour degree by $7-t$ and the lengths by $9-t$ and $11-t$. The charge equations permit $0\le t\le4$. Under the filtered Taylor Ward selection on the leading instanton map, the odd highest-weight degrees in the last two slots are $2-t$. Thus the selected component can occur only for $t=0,1,2$. This use of the Ward selection for an instanton is an additional assumption. [Charge equations](../../../su3_independent_prime_20261118/runtime/work/sp3-native.nosync/g2_full_classical_screen_v1.py) · [Taylor selection](../../b2_schur_complex.py)

| Level | $t$ | Incoming rational classical dimension | Incoming rank of protected products | Outgoing rational classical dimension | Outgoing rank of protected products |
|---:|---:|---:|---:|---:|---:|
| 30 | 0 | 5 | 5 | 1 | 1 |
| 30 | 2 | 23 | 23 | 18 | 18 |
| 32 | 0 | 7 | 7 | 1 | 1 |
| 32 | 2 | 43 | 43 | 28 | 28 |

At $t=0$ the public finite-colour classical ranks bound the level-30 dimensions by five and one. Rational Cartan product minors attain both. At level 32 the upper bounds at lengths $8,9,10,11,12$ are $1,3,8,8,1$, with zero at all other lengths. The even colour-parity Euler index is 8. At lengths 8 and 12 there is an even product, and at lengths 9 and 11 there are three and seven independent odd products. Therefore the even dimension at length 10 is at most $8-1-1+(3-3)+(8-7)=7$. Seven independent even products attain it. The length-12 product attains the full outgoing bound. Each product uses components of the defining quadratic and cubic graviton traces. [Rational Cartan minors](bmn/) · [Fortuity level-30 data](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/bmn_hw_cohomology_j30.json) · [Fortuity level-32 ranks](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/hw_cohomology_j32_records/su3_c_8_4_4.jsonl)

At $t=2$ the two selected diagrams are $(5,5,1,1,1)$ and $(5,5,1,1,1,1)$. The incoming lengths are seven and eight. On each of the five strip components the classical matrices are evaluated at prime $2147483647$. The colour-parity dimensions are computed independently. The Schur projector commutes with reduction because the Taylor degrees 13 and 14 have invertible factorials. Consequently the projected modular ranks give the following rational upper bounds. The other colour parity has dimension zero in both cases. [A6 job 7658864](7658864/receipt.json) · [Classical rank driver](../../su3_f1_incoming_ranks.py)

$$
2059-467-1569=23,\qquad 4501-1156-3302=43.
$$

The 18 complete product enumerations at the same prime give the ordinary and spin-raised Cartan ranks. Subtracting the raised ranks isolates spin highest weights. The horizontal-strip formula extracts the selected Schur multiplicities,

$$
62-54+15-0+0=23,\qquad 107-88+24-0+0=43.
$$

Every classical boundary evaluates to zero on these Cartan letters. The two lower bounds equal the classical upper bounds, so these are rational product spans. The earlier outgoing product calculations give 18 and 28. [A8 job 5692604](5692604/receipt.json) · [Outgoing ranks and projection proof](../cartan_products/README.md) · [Resolved ranks](resolved.json)

If the graviton components and their cohomological primary products have exact closed completions, every incoming class at grades zero and two is represented by a closed product. Under the filtered Ward selection an incoming leading instanton source must then have a grade-one component outside the protected span. Its selected cochain dimensions are 5365 and 10783 at levels 30 and 32. The outgoing grade-one dimensions are 11043 and 18911. The grade-one cohomologies and the instanton matrix elements remain open. Closure of the outgoing product representatives alone allows a possible incoming differential into their span. [Exact currents and primary product](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

From the repository root, using NumPy, SciPy and python-flint, replay the ranks and determinants with

```sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_f1_partners"
python3 "$d/su3_f1_partners.py" --out /tmp/su3-f1-partners.json
cmp "$p/resolved.json" /tmp/su3-f1-partners.json
```

Pass `--upstream /path/to/fortuity-ai/data/su3` to verify both upstream Git blobs. Pass `--rebuild-products` to regenerate the complete 18 modular Cartan matrices. To recompute the ten classical rank records at their saved seeds, build the native engine and run

```sh
c++ -O3 -std=c++17 -pthread -I"$d/taylor_hooks/native" "$d/taylor_hooks/unitary_parity.cpp" -o /tmp/su3-parity
python3 "$d/su3_f1_incoming_ranks.py" /tmp/su3-parity /tmp/su3-f1-incoming --threads 8 --job 7658864
diff -r "$p/7658864" /tmp/su3-f1-incoming
```
