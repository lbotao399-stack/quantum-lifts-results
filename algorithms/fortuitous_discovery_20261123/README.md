# Classical fortuitous certificates and exact enumeration

Independent research track, 2026-11-23. Native repository field,
divided-power Taylor, and Q0 conventions are unchanged.

## Updated results and publication status

The [Sp(6) addendum](sp6_extension/README.md) supplies a further exact
classical result in paper sector `(11,0,0;3,3), Y=15`, at `j=40`:
fortuitous multiplicity `(183 - 61) - 121 = 1`, a 152833-term
highest-weight operator, its 98882-term defect, and a rational
111-probe single-trace dual. Its frozen operator and exact certificate
have passed the packaged portable verifier. The full finite-sector
cohomology dimension is not inferred.

The addendum includes the cycle-based Pfaffian implementation,
defect-specific inverse solver, complete 1961-block verification,
scalar-spectator reuse comparison, and first-identity-layer exclusion
argument. At the same native weight and operator length, the adjacent
Sp(7) sector has fortuitous dimension zero. This is a new certified
rank member relative to the preceding calculation; ten independent
new families have not been established. No loop-corrected survival
claim is made.

The two symplectic multiplicities are now **exactly one**, over Q, in
the displayed sectors. [The enumeration addendum](extension/README.md)
contains the complete183-dimensional relation frames and exact integer
upper certificates using Q0 and stable abelianization.
[The first-relation stability proof](extension/stability/README.md)
proves dimension67 for the first relation degree at every symplectic rank
at least4. This does not assert an all-rank fortuitous family.

The initial 172-file release was published at commit
`2dc1dec8eee723c44e9da10ac71b7648d5c20fe1`. It includes the research source, data, exact
certificates, pinned dependencies, included upstream source, enumeration
addendum and stability proof. The user explicitly authorized public
publication to [`research/fortuitous-dual-discovery-20261123`](https://github.com/lbotao399-stack/quantum-lifts-results/tree/research/fortuitous-dual-discovery-20261123/algorithms/fortuitous_discovery_20261123)
in `lbotao399-stack/quantum-lifts-results`. `PUBLICATION_STATUS.json`
records the release scope and current file count, including the Sp(6)
addendum. The final Git commit and branch verification
receipt is recorded separately from this immutable payload.

The original sections below retain the construction and audit of the
four original rational representatives; the addenda supply the new upper
bounds without changing their normalization.

## Results

| Gauge group | Paper highest weight and length | j | Certified fortuitous multiplicity | Exact witness |
|---|---|---:|---:|---|
| Sp(4)=USp(8) | (7,0,0;3,3), Y=11 | 32 | exactly1 | 17430 rational trace terms; 103-probe rational separating dual |
| Sp(5)=USp(10) | (9,0,0;3,3), Y=13 | 36 | exactly1 | 48231 rational trace terms; 40-probe rational separating dual |
| SO(15) | (7,7,3;0,0), Y=16 | 34 | exactly1 in the specified quotient | complete3-dimensional relation frame; two-probe dual |
| SO(15) | (9,5,3;0,0), Y=16 | 34 | exactly1 in the specified quotient | complete4-dimensional relation frame; two-probe dual |

The native computational order is `(n_plus,n_minus,q1,q2,q3)`, so the two
symplectic files use `(3,3,7,0,0)` and `(3,3,9,0,0)`. The scalar shapes are
the existing flavour highest weights. The SO(15) irreps are [0,4] and
[4,2], respectively, with dimensions15 and60.

These are classical Q0 results. Quantum lifting and all-loop BPS survival
have not been calculated for these four representatives.

## Exact symplectic evidence

| Quantity | Sp(4) | Sp(5) |
|---|---:|---:|
| Stable centre dimension | 11627 | 37085 |
| Stable incoming rank | 6940 | 22495 |
| Stable outgoing rank | 4687 | 14590 |
| Complete preceding relation dimension | 67 | 67 |
| Complete source relation dimension | 188 | 188 |
| Exact source-relation image rank | 121 | 121 |
| Genuine target identities used | 181 | 181 |
| Complete target relation dimension | 183 | 183 |
| Nonboundary augmented rank | 122 | 122 |
| Operator trace terms | 17430 | 48231 |
| Q0 defect trace terms | 15556 | 42480 |
| Least trace-coefficient denominator | 360 | 288 |

The Sp(5) formal highest-weight reconstruction used denominator576. After
expansion and combining equal trace terms, the least denominator is288.
This reduces fractions; it does not rescale the representative.

All defect coefficients were compared exactly to rational combinations of
structural finite-rank trace identities. The rational duals annihilate Q0
of every one of the188 complete source identities, and evaluate to1 on
the constructed defect. Thus each finite class is nonzero and cannot
come from stable cohomology. This argument alone proves fortuity; it does
not require assuming that a modular cohomology dimension lifts to Q.

The stable acyclicity calculations additionally give

```
11627 - 6940 - 4687 = 0,
37085 - 22495 - 14590 = 0.
```

Every identity is generated by antisymmetrizing ten fundamental indices
in dimension eight, or twelve indices in dimension ten, followed by
explicit highest-weight projection, Taylor-profile extraction, or actual
Q0. Source independence and exact character deficits certify completeness.
The final two target identities are not needed for positive existence.
The enumeration addendum now completes these two specified symplectic sectors, with exact dimension1.

## Exact SO(15) matrices

Rows are the ordered trace probes specified in `scalar/README.md`. The
first columns are the complete source relation images; the last column
is the specified scalar Pfaffian defect.

For (7,7,3):

```
      -57960       -7560       272160
   -12983040    -1693440     20321280
    32336640     2419200    -52254720

det = 4236912357801984000
dual = (1/181440, -1/40642560, 0)
dual * matrix = (0,0,1)
```

For (9,5,3):

```
       86940        49140        11340       -374220
    17388000      9828000      2268000     -24948000
     1179360       997920      -272160             0
   -19051200     -4536000     -1088640      35925120

det = 13523262578496405504000000
dual = (-1/249480, 1/49896000, 0, 0)
dual * matrix = (0,0,0,1)
```

Exact source/terminal character dimensions are2/3 and3/4. These small
matrices therefore prove global nonboundary, not merely independence
inside a sampled incomplete relation family.

The actual cochains are fully specified by the original Pfaffian marker
operators and the explicit flavour-equivariant odd-trace primitive. The
coefficient oracle in `scalar/minimal_scalar_certificates.py` computes any
requested coefficient. For example, in the original native convention:

```
(7,7,3): [AA.AA.AA.Aa.BB.BB.BB.CC] = 45360,
(9,5,3): [AA.AA.AA.AA.Aa.BB.BB.CC] = -62370.
```

## Strict rank support of the two SO defects

At SO(13), the same specified defects have exact Q0 preimages within the
finite trace-identity ideal. The complete terminal frames have dimensions
45 and90. Their preimages use44 and88 actual Q0-image rows; their
additional old-product coordinates are exactly0 and(0,0).

All50 and100 rational coefficient equations pass. The underlying integer
frames are reconstructed from sixteen primes using the a priori bound
`abs(entry)<2^248<2^256`, with CRT modulus greater than2^257. A seventeenth
prime is an additional check. The standalone rebuild verifies193800
residue entries in34 prime-case replays.

Relation ideals decrease under block embeddings. Consequently these
fixed representatives have zero **fortuitous component** at all odd
orthogonal dimensions at most13. They have nonzero fortuitous components
at15. At odd dimensions exactly17, their nonzero degree-seventeen defects
are no longer identities, so those same representatives are not closed.
Their fortuitous support among odd orthogonal groups is exactly SO(15).
This does not claim that their total smaller-rank cohomology classes vanish.

The birth mechanism is concrete: a defect remains a finite trace identity
while a smaller-rank identity that supplied its Q0 preimage ceases to be
an identity. This can create a new fortuitous component as rank increases.
It does not exclude a brane-fusion interpretation with additional scalar
blocks.

## Algorithmic conclusion

The finite/stable comparison is controlled by the connecting quotient

```
(finite trace identities intersect stable Q0 image)
/ Q0(finite source trace identities).
```

This keeps the universal lift obstruction, rather than choosing a metric
orthogonal complement to gravitons. A closed identity in the next degree
need not be in the stable Q0 image: this intersection must actually be
checked. In the two symplectic sectors the large stable spaces have
11627 and37085 centre dimensions, while the complete source relation
spaces both have188 dimensions. Sparse trace-support blocks supply the
stable-image test; exact reconstruction and rational duals finish the proof.
These dimension reductions are not wall-clock speedup benchmarks.

The scalar calculation goes further. `structure/UNIVERSAL_TWO_PROBE_DUAL.md`
derives a closed dual for every all-odd highest weight with smallest row
at least3, in odd orthogonal defining dimension `sum(shape)-2`. Only the
coefficients of the one- and three-cubic-trace probes are needed. The
source-completeness proof uses the first relation tensor functor, and the
coefficient identity has a finite exact interpolation certificate with a
20-by-20 parameter evaluation determinant110592. Direct coefficient
checks on19 highest weights and50 trace counts agree exactly.

This proves existence for that specified minimal scalar family; it does
not enumerate all fortuitous classes. The SO tower and its fusion
interpretation were already anticipated. The new methodological result
is the explicit parameter-uniform separating functional and its proof.

## Historical scope

Novelty is relative to the pinned data and current branches actually
inspected, documented in `novelty/`. The current manuscript's Sp(4)/Sp(5)
census stopped atj24. Its threshold table and omega construction covered
Sp(2)/Sp(3); no all-rank symplectic omega assertion was found there.
The two new symplectic results can belong to one extending family;
two unrelated physical families or full ring primitiveness are not claimed.

The original SO(13) j30 candidates and the SO(15) (7,5,5) control are
deliberately excluded from the new-class count. Their nonvanishing was
already implicit in earlier product certificates or smaller-rank
persistence. The two selected SO(15) weights pass the stronger exact
smaller-rank boundary comparison.

Historical audit snapshots retain their earlier provisional wording.
The completed mathematical receipts in `audit/` and the two selected
rank-birth proofs in `scalar/` supersede that earlier pending status.
No claim is made to have surveyed unpublished work.

## Reproduce

Requirements: Python3, python-flint, NumPy, SciPy; a C++17 compiler for
rebuilding the scalar CRT banks. The exact48 imported upstream Python
modules are included with source pins and SHA256 hashes, so no external
checkout or original workspace paths are required.

`search/` also contains the actual discovery implementations and their
source hashes. Large modular-search checkpoints are reproducible and are
not needed to verify the frozen rational proof data.

```
python replay.py dependencies
python replay.py scalar
python replay.py theory
python replay.py symplectic --rank 4
python replay.py symplectic --rank 5
python replay.py identities --rank 4
python replay.py identities --rank 5 --seconds 300
```

The Sp(5) identity replay saves a checkpoint; repeat it after a time cap
until `replay_outputs/sp5_identities.json` reports `passed: true`.
Only repeated deterministic graph expansions are memoized. The `.pkl.xz`
files decompress byte-for-byte to the original frozen exact data, with
round-trip hashes in `COMPRESSED_DATA.json`. There are no compiled
executables in this delivery.

`python replay.py all` runs the dependency, scalar, closed-form, and both
rational-witness checks. The longer full identity reconstructions remain
explicit commands above. Existing completed receipts are in `audit/`.

## Source pins

- Main upstream: https://github.com/lbotao399-stack/quantum-lifts-results
- Core trace algorithms: `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`.
- Native runtime snapshot: `082be7070d207e2cdecfd0a91b03d5e014fd303d`.
- Primary Fortuity inventory: https://github.com/davidsd/fortuity-ai/tree/4256f207986029540060fb44dd494746828f65ce
- Existing mixed-lift branch audited: `78d170a26afc2a16cae43f123df86d5880be6196`.
- Current manuscript audited: `All-loop-holomorphic-kernels@bba7a4857710c1baef433341f1428b63ba7642f0`.

The intended publication destination is the independent branch
`research/fortuitous-dual-discovery-20261123` in
`https://github.com/lbotao399-stack/quantum-lifts-results`.
Automatic approval review rejected commit creation despite the renewed
general upload authorization. It requires explicit confirmation of this
public destination and the complete code, data, certificates, included
upstream files, enumeration addendum, and stability proof.
No Research Discussion announcement or reviewer request was made.
