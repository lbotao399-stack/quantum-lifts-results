# Exact symplectic enumeration from identities and stable abelianization

This addendum upgrades the earlier existence certificates in the parent
package. The field, supercharge, Taylor coefficients, trace conventions,
and physical representatives are unchanged. All assertions concern
classical Q0. Quantum survival has not been calculated.

## Results

| Group | Paper sector | j | Exact fortuitous multiplicity | Exact full sector cohomology dimension |
|---|---|---:|---:|---:|
| Sp(4)=USp(8) | (7,0,0;3,3), Y=11 | 32 | 1 | 1 |
| Sp(5)=USp(10) | (9,0,0;3,3), Y=13 | 36 | 1 | 1 |

The full-sector conclusion also uses the parent package's independently
certified stable-centre acyclicity: 11627-6940-4687=0 and
37085-22495-14590=0. The fortuitous count itself follows from the new
upper certificates and the existing rational nonboundary witnesses.
The earlier frozen discovery files record the historically correct
lower bound; this addendum establishes equality.

The two ranks can belong to the same extending family. This is not a
claim of two unrelated primitive families, an all-rank symplectic theorem,
or a complete census outside the displayed highest-weight and length.

## The two missing identities

The previous proof used 181 genuine target identities in an exactly
183-dimensional relation space. A deterministic highest-weight
construction supplies the missing two directions:

1. One four-letter path with the remaining marked endpoints paired.
2. A first-degree Pfaffian relation multiplied by a separate quadratic
   trace, projected to the required highest weight.

Both recipes contain only two terms. Their exact trace expansions have
262 and58 terms at Sp(4), and674 and110 terms at Sp(5). The Sp(5) recipes
are re-expanded with twelve antisymmetrized fundamental indices in
dimension ten; this transports identity recipes only, not cohomology.
The complete183 rows are independent, and exact native characters give
dimension183. The source relation dimension remains188 and its Q0-image
dimension is121.

The separate trace factor in the second recipe matters: considering
connected path topologies alone left the relation space incomplete.
Neither of the added directions changes the fortuitous multiplicity.

## Which Cartan map is used

The input is a finite-rank trace identity regarded as a polynomial in the
**universal trace covering**, before imposing its finite-rank vanishing.
Replace its letters by supercommuting letters, with the original
Grassmann signs, and keep the resulting formal trace monomials. In the
symplectic Cartan, odd-length traces vanish. Even traces retain their
original normalization; no trace is divided by2.

This is stable abelianization. Evaluating the same polynomial on the
physical finite-rank Cartan would give zero for every trace identity and
would not provide this certificate.

The native rule `{Q0,Psi}=Psi^2` becomes zero in a supercommutative
algebra because Psi is odd. Equivalently, pair the two ordered splits in
`Q_letter`: their supercommutative products cancel. The executable checks
all127 admissible letters with Taylor exponents at most3 in each
direction, including every letter that can occur in these sectors.
Consequently stable abelianization annihilates Q0 of **every** universal
cochain in the relevant grade. This argument does not assume that Cartan
restriction detects all stable cohomology.

## The exact upper bound

Apply two maps to the complete183-dimensional target identity space:
actual Q0 and stable abelianization. Integer coefficient matrices and
rational reconstruction give, for both groups,

```
rank(Q0 on target identities) = 45,
rank(stable abelianization on target identities) = 24,
rank(the combined map) = 61.
```

The selected61-by-61 integer minors have determinants

```
Sp(4):  2^156 * 3^16 * 5
      = 19660266632258300811874950676834221032552013803846369280,
Sp(5): -2^113 * 3^27
      = -79188731731498164979739428589463199233557397504.
```

Thus the common kernel has dimension at most183-61=122. The parent
rational certificates already supply122 independent vectors in that
kernel: the121-dimensional image of source identities, and the defect
of one nonboundary operator. Both types are actual universal Q0 images,
so Q0 squared and stable abelianization annihilate them exactly. The
common kernel therefore has dimension exactly122.

Every target identity that is an actual universal Q0 image lies in this
common kernel. Conversely, the122 known independent universal Q0 images
exhaust it. For these two specified sectors, therefore,

```
target identities intersect universal Q0 image
  = target identities killed by both Q0 and stable abelianization.
```

After dividing by the121-dimensional image of source identities,

```
fortuitous dimension = (183 - 61) - 121 = 1.
```

This is a certificate of sufficiency in these sectors, established by
matching a rigorous upper bound with explicit rational lower witnesses.
It is not a universal assertion that Cartan-zero implies Q0-exact.

The individual ranks have exact rational upper certificates: every one
of the183 Q0 images is reconstructed from45 basis images, and every
abelianized row from24 basis rows, with equality checked on every trace
coefficient. The61-dimensional lower bound is an integer determinant,
not an inference from a modular nullity.

## Why closed trace identities overcount fortuity

In the target degree the relation ideal has cohomology dimension

```
183 - 45 - 121 = 17.
```

Stable abelianization has rank

```
61 - 45 = 16
```

on these relation cohomology classes. Those16 classes are already
nonzero in stable cohomology: the abelianization annihilates every
universal boundary. The remaining one is the connecting obstruction
of the finite-rank fortuitous class.

This identifies the actual splitting needed for enumeration: retain the
kernel of the map from relation cohomology to stable cohomology. Merely
finding a closed trace identity would retain17 directions here. No
Hilbert-space metric or orthogonal complement to gravitons is used.

## Computation and its limit

The new upper certificate needs183 identity polynomials, their Q0 images,
and their stable abelianizations. It does not require enumerating the
11627- or37085-dimensional stable centre. The existing rational witness
still supplies the indispensable lower bound.

As an independent cross-check, the previous full stable-image algorithm
was rerun with all183 target identities. It returned intersection
dimension122 at both ranks, and the original rational operators were
reconstructed with exact finite closure. These runs are recorded in
`full_stable_image_crosscheck_sp4.json` and
`full_stable_image_crosscheck_sp5.json`.

Exact characters have also been extended beyond the certified classes:

| Group | Native weight | Adjacent lengths | Relation dimensions |
|---|---|---|---|
| Sp(4) | (3,3,7,0,0) | 10,11,12 | 67,188,183 |
| Sp(5) | (3,3,9,0,0) | 12,13,14 | 67,188,183 |
| Sp(6) | (3,3,11,0,0) | 14,15,16 | 67,188,183 |
| Sp(7) | (3,3,13,0,0) | 16,17,18 | 67,188,183 |

Sp(6) and Sp(7) are character results only; no new representative or
cohomology multiplicity is certified at those ranks. The machine-readable
scan additionally records two completed Sp(8) degrees; its third degree
hit a time cap. No all-rank conclusion is inferred from this scan.

## Reproduce

Place this folder at `extension/` inside the parent four-certificate
package. Its existing Python requirements suffice. Run

```bash
python extension/verify_extension.py --rank 4
python extension/verify_extension.py --rank 5
```

Each command replays the parent's rational lower witness, reconstructs
the two added identities from antisymmetrizer recipes, verifies the
complete183-dimensional identity frame, checks the127 letter identities,
checks every rational map coordinate, and recomputes the exact61-by-61
integer determinant. To reuse a previously verified lower witness, add
`--upper-only`. To keep the extension elsewhere, supply `--base PATH`.

The parent full identity-generation replays remain available as
`python replay.py identities --rank 4` and
`python replay.py identities --rank 5 --seconds 300`.

`discovery_scripts/` preserves the actual exploratory programs with
hashes in `DISCOVERY_PROVENANCE.json`. They retain their historical
workspace dependencies. `verify_extension.py` is the portable proof
entry point and does not import the original discovery workspace.

The included successful upper-only receipts reuse the parent witnesses
already independently replayed during the original run. This addendum
does not change the four original physical representatives or their
coefficient normalization.
