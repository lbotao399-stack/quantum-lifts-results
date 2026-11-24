# Exact Sp(7)–Sp(16) fortuitous rank members and bounded character scans

## Result and precise scope

This package certifies **ten additional classical fortuitous rank/charge
classes**, one at each Sp rank7 through16, relative to the inspected
repository baseline below. They are ten members of **one growing-weight
family**. They are not ten independent families or ten newly identified
superconformal primary multiplets. No all-rank existence theorem or
loop-corrected survival claim is made.

The native convention is `Sp(r)=USp(2r)`. Native weights are
`(n_plus,n_minus,theta1,theta2,theta3)`; the paper orders them as
`(theta1,theta2,theta3;n_plus,n_minus)_Y`. `Y` counts native Taylor
letters, and `j=3(n_plus+n_minus)+2(theta1+theta2+theta3)`.
No physical letter, differential or trace normalization is changed.

| Group | Paper weight | Y | j | Exact fortuitous dimension |
|---|---|---:|---:|---:|
| Sp(7) | (13,0,0;3,3) | 17 | 44 | 1 |
| Sp(8) | (15,0,0;3,3) | 19 | 48 | 1 |
| Sp(9) | (17,0,0;3,3) | 21 | 52 | 1 |
| Sp(10) | (19,0,0;3,3) | 23 | 56 | 1 |
| Sp(11) | (21,0,0;3,3) | 25 | 60 | 1 |
| Sp(12) | (23,0,0;3,3) | 27 | 64 | 1 |
| Sp(13) | (25,0,0;3,3) | 29 | 68 | 1 |
| Sp(14) | (27,0,0;3,3) | 31 | 72 | 1 |
| Sp(15) | (29,0,0;3,3) | 33 | 76 | 1 |
| Sp(16) | (31,0,0;3,3) | 35 | 80 | 1 |

At each of these ranks, the independently computed identity dimensions
at lengths `(2r+2,2r+3,2r+4,2r+5)` are exactly

```text
67, 188, 183, 45.
```

The complete source-boundary map has rank121; the combined target map
`(Q0, stable abelianization)` has rank61. Thus

```text
dim fortuitous quotient = (183 - 61) - 121 = 122 - 121 = 1.
```

The distinction between a closed relation and an actual fortuitous
defect is numerically substantial here. Closure alone leaves
`183-45-121=17` relation-cohomology directions. The additional stable
abelianization obstruction has rank `61-45=16` on that quotient, so
only `17-16=1` direction has a universal Q0 preimage. Counting closed
relations alone would overcount by16 at each of these ranks.

Each frozen certificate contains an integer combination of genuine
Pfaffian identity recipes and a rational trace-coefficient functional.
The defect has exact `Q0=0` and full stable abelianization0; the
functional pairs to0 with all188 source-relation boundaries and to1
with the defect. Complete identity dimensions and nonzero minors make
the selected closure probes injective on the entire outgoing space.

**The compact certificate proves existence of a universal preimage;
the physical operators at these ten ranks have not been expanded.**
The preimage criterion is proved in `STABLE_ABELIANIZATION.md` from the
published stable classification. These records must not be represented
as explicit fully expanded physical-operator files.

At the same weight and field length, the fortuitous quotient at every
larger Sp rank is zero: at rank `r+1`, `Y=2r+3=2(r+1)+1`, inside the
proved exclusion range in `../sp6_extension/FIRST_LAYER_EXCLUSION.md`.
Every proper positive-length product has factors of length at most
`2r+1`, where all cohomology is stable. Consequently these classes
cannot be products of lower-length cohomology classes. This does not
classify their entire superconformal multiplets.

## What changed algorithmically

### Stable abelianization gives a sufficient preimage test

The map sorts letters inside each trace with their Grassmann signs and
keeps independent formal trace generators. It is applied before the
finite-rank quotient; it is not finite Cartan evaluation. The argument
in `STABLE_ABELIANIZATION.md` proves, for the ordinary SU/SO/Sp trace
covering with all native Taylor letters,

```text
ker Q0 intersect ker(stable abelianization) = image Q0.
```

The input is the complete graded-symmetrized single-graviton derivative
family of Chang–Yin2013, Section4, equations(4.5)–(4.7), also stated in
Chang2026, Definition2 and Proposition5. For SO/Sp the proof takes the
primitive trace-reversal quotient first, then the graded-commutative
algebra on the resulting primitive complex. It does not assume an
unstated super-relative extension of an ordinary dihedral theorem.

This replaces the construction of a huge universal inverse-Q matrix
when the immediate task is to certify existence and nontriviality.
It does not replace inverse-Q when an explicit physical operator or
correlation function is required.

### Compute the missing character directly

`direct_character_deficit.py` uses the Littlewood identity
`h_m[h_2] = sum_(lambda partition m) s_(2 lambda)` for the largest
bosonic species. The remaining small plethysm is treated in native
power sums, using exact skew Murnaghan–Nakayama characters. Only
paired-row Schur shapes excluded by the finite Sp rank are summed.
All factorial denominators and the native highest-weight subtraction
are retained.

Twenty independent species/rank controls agree with the original
`stable dimension - finite dimension`. Ranks2 and3 reproduce their
different deficits `[47,123,107,21]` and `[65,181,174,42]`. Every rank4
through16 is computed separately; no constant sequence is extrapolated.
See `DIRECT_CHARACTER_DEFICIT.md` for the complete enumeration argument.

### Transport coefficients instead of expanding higher-rank words

`SPECTATOR_TRANSPORT.md` proves the exact scalar-pair coefficient rule
for adding `k=r-5` native `(A,A)` open paths to a Sp5 recipe. If the
baseline recipe has M such paths and a probe has t quadratic AA traces,
the coefficient multiplier is

```text
(-1)^k binom(M+k,k) / binom(t+k,k).
```

Three target rows,174,176,178, have a distinguished closed AA trace.
Their exact correction is included; omitting it would invalidate the
transport. Connected outgoing probes have rank39. Adding the necessary
two-trace probes raises the rank to45, equal to the complete outgoing
identity dimension. This certifies closure of the full polynomial.

Stable abelianization is still evaluated in full at each actual rank.
Its output support grows from2,466 trace words at Sp7 to174,068 at Sp16.
The method is not claimed to have constant cost at arbitrary rank.

## Executed verification and measured runtime

All certificate arithmetic is exact integer or rational arithmetic.
The prime1000003 only selects and certifies nonzero minors; positive
classes are never inferred from a modular kernel alone.

| Check | Exact completed amount |
|---|---:|
| Complete source and target seed rows, fully re-expanded | 371 |
| Nonzero seed polynomial coefficients compared | 1,797,734 |
| Seed matrix entries independently recomputed | 76,812 |
| Highest-weight recipe rows checked | 371 |
| Abstract differential/native comparisons | 105 |
| Coefficient-oracle comparisons | 1,150 over62 recipes |
| Full commutative/native-expansion comparisons | 125 |
| Independent transport-coefficient controls | 165 |
| Stable-cohomology control blocks | 58 over1,389 cochains |
| Primitive reversal/involution/chain checks | 588 |
| Frozen frontier certificates replayed | 10 |
| Exact source-boundary pairings checked in frontier replay | 1,880 |

The core Sp7–Sp16 certification invocation records
`elapsed_seconds=380.5390892730002`, after preparing the finite seed.
The separate frozen-certificate replay records
`elapsed_seconds=386.3446145959997`. Full seed re-expansion records
`elapsed_seconds=203.11253998099983`; seed-map verification records
`elapsed_seconds=5.136695172000145`.
These are scoped wall-clock receipts on this workspace, not an end-to-end
development budget or a price/speed ratio against the entire cluster
campaign. Research and failed exploratory runs are also preserved.

## Bounded carpet scans and exclusions from novelty counts

`one_flavour_character_scan.json` contains all2,590 declared character
sectors at Sp ranks2–6: `theta2=theta3=0`, `n_plus>=n_minus`,
`n_plus+n_minus<=8`, `theta1<=2r+2`, `j<=42`, and source lengths
`2r+2` and `2r+3`. The recorded statuses are:

| Character result | Sectors |
|---|---:|
| Requires differential computation | 957 |
| No target identity | 271 |
| Empty source by letter weight | 1,344 |
| Empty finite source | 18 |
| Total | 2,590 |

A further declared differential subset has Sp ranks2–4, source length
`2r+3`, and identity dimensions at most700. All167 planned sectors were
attempted;161 have exact results and6 retain explicit time caps.
The capped Sp4 native weights are
`(4,4,5,0,0)`, `(6,2,5,0,0)`, `(4,3,6,0,0)`, `(5,2,6,0,0)`,
`(5,3,5,0,0)`, `(3,3,8,0,0)`. No capped entry is recorded as zero.

The ten positive small-rank sectors have total dimension12. They are
not added to the ten frontier rank members: several are already
catalogued, and remaining multiplet/descendant novelty is unresolved.
In particular the new exact multiplicity certificate at Sp3 native
`(4,4,3,0,0)`, `Y=9`, `j=30` has dimension2; this is an exact
enumeration result, not a claim of two new independent families.

The even-orthogonal first scalar layer was evaluated for
SO4,SO6,SO8,SO10,SO12,SO14,SO16:122 scalar partitions,30 positive
certificates. Their existence follows from prior Pfaffian-marker and
odd-trace-primitive constructions, so all30 are explicitly excluded
from discovery counts. SU/SO/Sp/G2 native character API controls and
the broader frontier inventory are supplied, but there is no claim of
a complete differential census across every classical and exceptional
group. Epsilon sectors and exceptional invariant generators require
their own covering.

## Prior-art baseline and sources

The earlier publication already contains Sp4/5/6 members. The bounded
primary-source audit found no positive all-rank Sp omega theorem or
published Sp7–Sp16 certificates in the inspected sources. This is a
repository-relative and source-bounded novelty statement, not a claim
about unpublished work.

- Authorized research-branch base:
  `b5907e39ac2c30dc03aed6aa5041a5d65fad532b`.
- Latest inspected repository main:
  `19cd8718b752b04104d6ff0becee42d580b27ca3`.
- Inspected review branch:
  `c778cd7c3ceb72ad8082d791b2aee5f0df082164`.
- Upstream fortuity-ai:
  `4256f207986029540060fb44dd494746828f65ce`.
- Detailed source inventory and exclusions:
  `inventory/FRONTIER_INVENTORY.md`, `inventory/sector_lineage_20261124.md`.
- Chang–Yin stable representatives:
  https://arxiv.org/html/1305.6314v2#S4 .
- Chang, Proposition5 and Definition2:
  https://arxiv.org/html/2604.23549v1#S4.SS1 .
- Symplectic second fundamental theorem:
  https://arxiv.org/abs/1504.02586 .
- Littlewood identity:
  https://arxiv.org/abs/1506.02755 and https://arxiv.org/abs/2108.12766 .

## Reproduce

Clone the authorized research branch and enter this directory. Keep the
published parent package in place. Its `requirements.txt` pins the
dependencies; Python3.12 was used. `package_runtime.py` discovers
`../replay.py`, or uses an explicit `FORTUITY_BASE_PACKAGE`.

```bash
python -m pip install -r ../requirements.txt
python restore_data.py
python verify_seed_maps.py
python verify_seed.py
python verify_transport.py
```

`verify_transport.py` consumes the frozen vectors and physical duals;
it does not solve replacement witnesses. The receipts are under
`transport/`. To reconstruct the finite seed and rerun discovery:

```bash
python prepare_transport.py
python validate_transport_formula.py
python direct_character_deficit.py --controls
python certify_transport.py --ranks 7 8 9 10 11 12 13 14 15 16
```

The release stores every small-rank frozen relation recipe, its row
order and its exact polynomial hash in `relation_recipes.tar.xz`.
It contains497 frame records with45,065 rows, representing16,516,222
nonzero polynomial coefficients, in3,171,928 compressed bytes.
Expanded relation matrices and their Q caches are reproducible and
are omitted from the compact archive. To recover a sector:

```bash
python rebuild_relation_frames.py --sector sp2_4_2
```

Use `--all` to regenerate every saved frame, including capped partial
frames. `one_flavour_relations.py` regenerates Q caches when a sector is
reanalyzed. `run_relation_scan.py` resumes the declared subset; a capped
sector can be continued by invoking `one_flavour_relations.py` directly
with a larger explicit `--seconds` value. The previously attempted
Sp7 full expansion is represented by its log and input provenance;
the unchanged parent `sp6_extension/build_rank.py --rank 7` regenerates
that alternative attempt.

`restore_data.py` verifies archive chunks and every restored file hash.
The archive inventory records the omitted generated-cache hashes as
well. No separate downloaded textbook PDF or page-image assets are
needed for any replay.

## Next research boundary

The exact transport formulas now permit systematic searches at larger
rank without expanding every noncommutative trace word. A next
mathematical objective is to prove an all-rank existence statement;
the present finite-rank certificates do not establish it. To reach ten
independent new families requires new highest-weight/multiplet data,
not merely more members or descendants of this family. Computing
physical observables and loop corrections also requires additional
input beyond these classical quotient certificates.
