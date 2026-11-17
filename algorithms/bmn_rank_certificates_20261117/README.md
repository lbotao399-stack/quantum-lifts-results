# Direct highest-weight BMN rank certificates

This release contains the November 17, 2026 `Sp(3)=USp(6)` and `G2` Cartan-product calculation at level `j=2(E+J_L)=30`: 41 product sectors per group, grouped into 12 highest-weight columns. The field prime is 1000003.

## Reproduce

Use Python 3.13 and numerical packages built for that same interpreter:

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python runtime/work/sp3-native.nosync/replay_direct_index_evidence_v1.py --evidence evidence --full-minors
.venv/bin/python verify_index_tables.py
```

Run these commands from this directory. `evidence/replay.json` is regenerated. The saved payloads use Python pickle and are intended for this repository's hash-verified data only.

The full replay verifies 209 payload hashes, reconstructs all 121237 entries of 62 nonempty minors with the independent dense matrix trace evaluator, recomputes their nonzero determinants, checks 7168 additional probe coordinates, checks exact closure of 185 distinct generators, checks 141 selected modular highest-weight rows, and verifies the rational `G2` torus realization. The second command checks source hashes, dimension sums, parity signs, and every index-bound entry. It does not rerun the character-integral calculation that supplied the archived exact dimensions.

## Calculation

1. Build products of exactly `Q0`-closed even-length trace generators.
2. Impose both flavor-raising equations on the formal product coefficients before evaluating. Compute their simultaneous kernel modulo the prime.
3. Lift a modular kernel row to integer coefficients and apply the existing rational extremal projector. Its denominators are invertible at this prime. The result is exactly highest weight and `Q0`-closed, with the same reduction modulo the prime.
4. Evaluate the projected products on Cartan matrices. A nonzero minor gives a rational lower bound on independent cohomology classes, since every `Q0` boundary restricts to zero.
5. Subtract this lower bound from the exact even-minus-odd-length index. A negative answer forces odd-length cohomology. A nonnegative answer leaves cancelling classes unresolved.

The index subtraction removes only the exhibited even-length closed subspace. The result is an upper bound on the index of the quotient by that subspace; it is not a complete calculation of the Cartan kernel. In particular, do not apply signed Weyl differences to sampled weight-space rank lower bounds: subtraction does not preserve a one-sided bound. The archived direct highest-weight minors avoid this problem.

## Results

Each entry below is `(exact index, Cartan-product rank lower bound, difference)`.

| Highest weight | Sp(3) | G2 |
|---|---|---|
| (11,2,2) | (77,77,0) | (14,14,0) |
| (10,3,2) | (161,161,0) | (24,24,0) |
| (9,4,2) | (211,211,0) | (25,25,0) |
| (9,3,3) | (111,110,1) | (12,11,1) |
| (8,5,2) | (192,192,0) | (20,20,0) |
| (8,4,3) | (183,181,2) | (19,17,2) |
| (7,6,2) | (113,113,0) | (12,12,0) |
| (7,5,3) | (163,161,2) | (21,18,3) |
| (7,4,4) | (87,86,1) | (13,4,9) |
| (6,6,3) | (67,67,0) | (11,8,3) |
| (6,5,4) | (75,74,1) | (15,5,10) |
| (5,5,5) | (12,12,0) | (4,0,4) |

No negative difference occurs. This phase certifies no new fortuitous generator or quantum lift. The ongoing undeformed N=4 campaign ledger remains: 5 audited `d1` highest-weight pairings, 3 audited fortuitous-generator families, and 0 certified `d2`, `d3`, or nonperturbative lifts. The family count does not assert global ring minimality or absolute literature novelty.

## Files and provenance

- `evidence/`: 41 shared product families, 82 group-sector results, 82 probe/minor payloads, 2 index gates, 2 exact-index inputs, and a fresh replay record.
- `evidence/SHA256SUMS.json`: hashes of the 209 scientific input payloads; generated replay output is excluded.
- `SOURCE_MANIFEST.json`: hashes of the bundled Python runtime.
- `runtime/`: the production direct-product algorithm and the independent replay code, retaining relative imports.

The vendored `fortuity` modules originate from [davidsd/fortuity-ai](https://github.com/davidsd/fortuity-ai), in particular its graded trace words, BMN generators, finite-field routines, and exterior-algebra matrix evaluator. The original BMN vendor provenance pins upstream commit `4256f207986029540060fb44dd494746828f65ce`. Original file headers are preserved. The exact bundled file hashes are recorded in `SOURCE_MANIFEST.json`; no upstream license file was present in the imported snapshot, and this release does not assign it a new license. The direct highest-weight product kernel, group-specific evaluation, rational certification gates, and replay packaging are the additions here.

## Nonabelian G2 initial rank evidence

The `g2_initial_rank_evidence/` directory adds the initial singlet complex at highest weight `(5,5,5)`, with six saved `Q0` rank minors and the original evaluation probes. Reproduce its audit with:

```sh
.venv/bin/python runtime/work/sp3-native.nosync/replay_g2_initial_rank_evidence_v1.py --evidence g2_initial_rank_evidence
```

This audit recomputes all six determinants, the exact character dimensions, rational highest-weight projections and symbolic `Q0` nilpotence for 18 selected source recipes, and 1056 dense-matrix probe coordinates. The remaining minor entries are retained production evaluations; this audit does not independently reconstruct every entry.

| Spin label | Exact dimension | Outgoing rank lower bound | Incoming rank lower bound | Cohomology upper bound |
|---:|---:|---:|---:|---:|
| 0 | 85 | 0 | 85 | 0 |
| 1 | 416 | 85 | 331 | 0 |
| 6 | 1387 | 1066 | 321 | 0 |
| 7 | 393 | 321 | 72 | 0 |
| 8 | 82 | 72 | 10 | 0 |
| 9 | 10 | 10 | 0 | 0 |

The spin label is the first coordinate of the archived BMN sector `(spin,5,5,5)`. Subtract the two rank lower bounds from the exact dimension. These six sectors are empty. Positive remainders in the four other sectors remain upper bounds, with no new lift or generator counted. Degree continuations must be merged in common evaluation frames before updating a full differential rank.

The production G2 algorithm is `runtime/work/sp3-native.nosync/g2_bmn_complex_v1.py`; its native trace and `Q0` evaluators are included as C++ source. The independent replay above needs no compiled native evaluator. To build the production backend on the compute host, from `runtime/work/sp3-native.nosync/`:

```sh
c++ -O3 -std=c++17 -shared -fPIC -DTRACE_FRAME_MATRIX_SIZE=7 -DTRACE_FRAME_PRIME=1000003 bmn_derivation_rank_v1.cpp -o bmn_derivation_rank_v1_d7_p1000003.dylib
```

The `.dylib` filename is the driver's fixed library name on both macOS and Linux. Production workers require a scheduler allocation; source completeness is asserted only at exact character-dimension saturation.
