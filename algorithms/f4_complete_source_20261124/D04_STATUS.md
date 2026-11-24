# Seventh-power second map: open

The complete source and first-correction counts are

```text
217297 -> 62866 -> 8553
23895 - 2959 = 20936 homogeneous first corrections.
```

The first-page target pairing is `-194893931232000`. Its complete rational
Ward replay checked all 2704 ordered scalar inputs, 1378 symmetric scalar
columns, 1117 differentiated-scalar components and 18948 mixed-ghost
components. Of the differentiated-scalar equations, 1074 have nonzero terms
that cancel exactly. The 68 homogeneous coefficients are retained.

The intermediate classical map has full terminal rank 2959. The freshly
rebuilt order-68 homogeneous matrix has determinants

| Prime | Determinant |
|---:|---:|
| 1009 | 909 |
| 1013 | 962 |
| 1000003 | 176568 |

Thus every classical second-source cycle admits a first correction, and all
20936 homogeneous first corrections remain free.

The public scalar Ward reduction has exact character counts

```text
129428650 - 40998203 = 88430447
641 + 1614 + 1387 = 3642
3642 - 1 = 3641.
```

The scalar quotient has 35 distinct irreducible colour types. Every type and
its complete invariant-tensor multiplicities in the three dotted channels are
listed in `seventh_ward_complete_inventory.json`. This is a character
inventory. It is not an explicit native basis for those multiplicity spaces.
The replay uses the recorded 54 scalar-weight rank saturations; those large
native matrices were not regenerated here.

The remaining differentiated-scalar colour shapes are `(3,2,1)` and
`(2,2,1,1)`, of unreduced dimensions 438534720 and 223835430. The
`(2,2,2)` summand is the scalar part of a total derivative; its unreduced
dimension is 129428650. These large numbers are colour-module dimensions
before the full bracket and invariance conditions, not dimensions of physical
cohomology.

The complete unresolved calculation is to impose every mixed-ghost equation
and both remaining differentiated-scalar types on the full Ward extension,
including all invariant multiplicities and homogeneous corrections. It must
test the full intersection with the scalar-bracket image. The previously
tested Albert-word families cannot replace this system.

No nonzero or zero full second-map coefficient is certified here. D04 remains
open; the first-page and intermediate calculations are its verified inputs.
