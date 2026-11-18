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
