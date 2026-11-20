# The singlet at doubled target charges (0,0,8,8,8)

At derivative bidegree (4,4) all letters are gauginos. Here the classical supercharge increases length by one. At lengths 6, 7, 8 the complete rational result is:

| Group | Cochain dimensions | Incoming ranks at lengths 6, 7, 8 | Outgoing ranks | Cohomology |
|---|---|---|---|---|
| SO(7) | 10, 12, 9 | 6, 4, 8 | 4, 8, 0 | 0, 0, 1 |
| USp(6) | 10, 12, 8 | 6, 4, 8 | 4, 8, 0 | 0, 0, 0 |

The last SO(7) line is the fourth power of the gaugino bilinear. Appending this trace polynomial to the eight classical images adds rank one for SO(7), zero for USp(6). The SO(7) class is the target of the previously known first differential. This particular classical difference is therefore absent on the second page. Other derivative bidegrees at the same physical charges are separate calculations.

The replay enumerates all formal trace words at lengths 2 through 8, applies the rational highest-weight projector and evaluates them on integral antisymmetric or symplectic matrices in an exterior algebra. At each length the selected minor has the full independent character dimension. This makes evaluation injective on the cochains over the rationals. The full differential, its consecutive zero products, and all coordinates of the fourth power are checked over the rationals. For length eight the integral trace-product bound is `8! * 7^8 < 2^63`; all sums with projector coefficients are taken in FLINT over the rationals.

From the repository root, using the same Python requirements as the level-34 computation:

```sh
python3 algorithms/dual_singlets_20261118/pure_gaugino_ranks.py --output /tmp/pure_gaugino_ranks.json
cmp /tmp/pure_gaugino_ranks.json algorithms/dual_singlets_20261118/pure_gaugino_ranks.json
```

[Exact rational matrices](pure_gaugino_ranks.json) · [Cochain dimensions for all five derivative bidegrees](charge_maps.json) · [Known incoming target and next SO(7) pairing](../so7_cartan_target_20261118/README.md).

## Zero derivative degree

At length 12 all letters are scalars. The square of the determinant of the three-by-three matrix of pairwise scalar traces is a flavour singlet of multidegree (4,4,4). Every classical boundary in this block evaluates to zero on three commuting scalar matrices: the only contributing letter replacement is a scalar commutator. On three independent Cartan generators the matrix of pairings is diagonal, with entries −2 for SO(7) and 2 for USp(6). The squared determinant is 64 in both cases. This gives a rational class independent of all boundaries. [Integer Cartan calculation](scalar_upper_bound.py)

| Group | Length-10 image: minor order, determinant | Length-11 image: minor order, determinant | Length-11 cohomology | Length-12 cohomology |
|---|---:|---:|---:|---:|
| SO(7) | 189, 956309 | 87, 508222 | 276 − 189 − 87 = 0 | 88 − 87 = 1 |
| USp(6) | 186, 910259 | 86, 455275 | 272 − 186 − 86 = 0 | 87 − 86 = 1 |

The determinants are modulo 1000003. Since the entries come from integral trace evaluations, they bound the rational ranks from below. Consecutive images compose to zero, so the two lower bounds, whose sum is the entire length-11 cochain dimension, are exact. Derivative bidegrees (1,1), (2,2) and (3,3) are separate blocks.

```sh
python3 algorithms/dual_singlets_20261118/scalar_upper_bound.py
python3 algorithms/dual_singlets_20261118/scalar_neighbors.py --replay algorithms/dual_singlets_20261118/scalar_minors.json.gz
python3 algorithms/dual_singlets_20261118/replay_scalar_signs.py
```

[Four minors and trace recipes](scalar_minors.json.gz) · [Integral sign check on 278 distinct recipes](scalar_sign_replay.json).

## Terminal singlets and the mixed length-eight classes

For derivative bidegrees (0,0) through (4,4), the terminal lengths are 12, 11, 10, 9, 8. The five weight-one letters are closed. The 12 weight-two letter rules generate their allowed graded commutators; the predecessor contains exactly one weight-two letter. This was checked for all 40 feasible dominant and shifted weights in the character calculation. After imposing the commutators and before finite-colour trace identities, the singlet bounds are 4, 0, 6, 0, 1. At (2,2), five independent finite-colour relations reduce the upper bound to one: the integer minors are 8192 for SO(7) and −4 for USp(6). The Cartan witness evaluates to 64 and −64, respectively. Together with the separate results above the **terminal classical dimensions are (1,0,1,0,1) for SO(7) and (1,0,1,0,0) for USp(6)**.

[Letter rules](terminal_letters.json) · [Commutator quotient](terminal_commutators.json) · [Orthogonal integer relations](orthogonal_replay.json) · [Symplectic Pfaffian relations](symplectic_relations.json) · [Cartan evaluation](mixed_cartan.py).

At bidegree (3,3), one length **below** the terminal length, the exact formal incoming rank is 225. Its 335 integer dependencies, along with an integer minor, are checked directly. Formal outgoing rank is 85. The formal representatives contain 30 and 18 trace words. The zero terminal quotient gives the full outgoing rank in each finite group. Combining these ranks with the finite-colour minors gives:

| Group | Classical cohomology at length 8 | Ranks: boundaries, two formal classes, then the known $F$ | Augmented minor mod 1000003 | Augmented minor mod 1000033 |
|---|---:|---|---:|---:|
| SO(7) | 309 − 81 − 225 = 3 | 225, 227, 228 | 392950 | 372431 |
| USp(6) | 306 − 79 − 225 = 2 | 225, 227 | 447499 | 659562 |

The nonzero modular minors are rational lower bounds because the evaluated rows are integral after clearing denominators. The preceding upper bounds make the dimensions exact. The third SO(7) line is accounted for by the previously counted $F$, whose first obstruction hits the fourth power of the gaugino bilinear. This identifies the source of the classical mismatch, without producing a new generator or a nonperturbative obstruction. The two common classes' higher differentials are not determined by this comparison.

```sh
python3 algorithms/dual_singlets_20261118/mixed33_exact.py
python3 algorithms/dual_singlets_20261118/mixed33_classes.py --group SO_odd --prime 1000003 --output /tmp/so7-classes.json
```

For the other certificates set `--group Sp` or `--prime 1000033`. [Formal dependencies, finite-colour minors and seeds](mixed33/) · [Stored 560 source images](mixed33_incoming.json.gz) · [The independently certified first differential](../so7_cartan_target_20261118/README.md).

## The remaining mixed singlets at doubled charges $(1,1,8,8,8)$

The spin-singlet source at derivative bidegree $(t,t)$ has length $11-t$. The classical differential raises that length by one. At $t=1,2$ the terminal cohomology dimensions are $0,1$, respectively, fixing the outgoing ranks. A8 job `5670267` evaluated the complete formal incoming singlet images in both groups. Each modular minor is a rational rank lower bound because the source recipes and signed trace evaluations lift to integers.

| $t$ | Group | Middle dimension | Incoming minor order | Outgoing rank | Rational cohomology upper bound |
|---:|---|---:|---:|---:|---:|
| 1 | $SO(7)$ | 1,419 | 1,105 | 310 | 4 |
| 1 | $USp(6)$ | 1,400 | 1,091 | 305 | 4 |
| 2 | $SO(7)$ | 1,411 | 1,084 | 327 | 0 |
| 2 | $USp(6)$ | 1,393 | 1,075 | 318 | 0 |

The two input spaces comprise 3,153 and 2,889 formal singlets. The eight pivot lists, determinants, seeds and matrix hashes cover both primes. Run `mixed_middle.py --t 1 --work /tmp/dual-t1 --workers 16` and repeat with `--t 2` and a separate directory. [A8 rank records](mixed_middle/a8_5670267.json.gz) · [Exact singlet recipes and finite-colour rank algorithm](mixed_middle.py)

At $t=3$, let $S$ be the gaugino bilinear and use the stress-multiplet currents of (4.83)–(4.85). Multiplication by $S$ takes the two six-letter products in the [current calculation](../G2_current_product_span_20261120.md) to this eight-letter sector. In the order of the two archived formal cycles, their superfield coefficients are $(0,-1/12)$ and $(-1,0)$. The equalities hold term by term in the 18- and 30-term polynomials. The physical current products, evaluated directly on two rank-three Cartan jets in each group, give

$$
\det\begin{pmatrix}-32/3&0\\8/3&48\end{pmatrix}=-512.
$$

Both evaluations annihilate every classical boundary. The two products are closed under the exact supercharge by the supersymmetric primary product. They span the two possible outgoing sources after the previously established first differential on the third $SO(7)$ line. Possible incoming boundaries can only reduce the surviving source dimension. [Source identifications and Cartan matrices](mixed33/current_products.json) · [Exact replay](mixed33_currents.py) · [Physical currents](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

At $t=1$ the ten-letter formal trace space gives $1597-1201-389=7$ classical singlets before finite-rank trace identities. A8 job `5670358` saved all seven rational representatives. Four Cartan evaluations on these classes give minor $3145728$ in each group. On the same jets four products of graded-symmetric single traces have determinant $-75497472$ in each group. All 124 evaluated classical letter variations vanish. The product minor attains the two finite-colour upper bounds above. [Representatives](mixed11/classical.json.gz) · [Products, contractions and probes](mixed11/current_product_test.json) · [Slurm receipt](mixed11/a8_5670358.json)

| Group | $t=0$ | $t=1$ | $t=2$ | $t=3$ | $t=4$ |
|---|---:|---:|---:|---:|---:|
| $SO(7)$ classical source dimension | 0 | 4 | 0 | 3 | 0 |
| $USp(6)$ classical source dimension | 0 | 4 | 0 | 2 | 0 |

The four $t=1$ products use single-graviton representatives in protected multiplets. Assuming their exact completions and the supersymmetric primary product, these four classes admit closed completions. At $t=3$ two protected products span the possible outgoing survivors in each group. This exhausts the tested spin and flavour singlet sources at doubled charges $(1,1,8,8,8)$. The statement concerns outgoing lifts and uses no instanton integral. [Single-graviton construction](https://arxiv.org/abs/1305.6314) · [Primary product](https://arxiv.org/html/2509.16737#S2)

```sh
python3 algorithms/dual_singlets_20261118/mixed_formal.py --t 1 --first-plan /tmp/dual-t1/plan.json.gz --work /tmp/dual-t1-out --workers 16
python3 algorithms/dual_singlets_20261118/mixed11_products.py --formal algorithms/dual_singlets_20261118/mixed11/classical.json.gz --out /tmp/products-t1.json
cmp /tmp/products-t1.json algorithms/dual_singlets_20261118/mixed11/current_product_test.json
```
