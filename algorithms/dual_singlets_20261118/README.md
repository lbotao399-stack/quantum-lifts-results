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
