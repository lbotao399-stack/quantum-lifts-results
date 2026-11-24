# Public claim ledger and finite coverage inputs

**The fixed baseline semantic/status audit is complete.** Its preserved
receipt records 244 canonical statements and attributed evidence packets,
all 255 distinct public texts, all 764 immutable document/issue occurrences,
zero unresolved semantic metadata and zero missing references. Its eight
deliberate corruptions are rejected. The current release gate additionally
requires the finite B05, F02/F03 and D03 publication addendum described in
`RELEASE_SCHEMA.md`; its authoritative status is `STRICT_GATE_RECEIPT.json`
and the generated `CURRENT_SUMMARY.md`. Historical proof limitations and open
physical questions retain their stated status. This audit does not re-prove
every historical result or award itself a roadmap point.

## Files

| File | Meaning |
|---|---|
| `canonical.py` | Curated semantic claims, exact native charges, proof status and supersession decisions |
| `claims.json` | Generated canonical records, including the five N4 records imported from the public nonperturbative claim index |
| `CURRENT_SUMMARY.md` | Generated current summary, exact counts and unresolved metadata |
| `source-audit.json` | Every ingested document occurrence, textual duplicate groups and semantic review state |
| `advertisement-spans.jsonl.gz` | Exact public paragraphs and their source locations; these are review units, not a claim count |
| `inventory/public_sources.json.gz` | Public source snapshots at immutable commits, plus the frozen Issue 1 text |
| `inventory/manifest.json` | Reachable source graph, link resolution, exclusions and source hashes |
| `coverage-blocks.json` | Coverage index binding the core and named-frontier manifests |
| `coverage/coverage_manifest.json` | Full finite core, counts, exact row hashes and independent alphabet check |
| `coverage/native_slots.jsonl.gz` | Every native core charge/length row |
| `coverage/physical_blocks.jsonl.gz` | Every complete physical core charge block and its native rows |
| `native-frontiers.json` | Every finite compulsory frontier input and all admissible neighbouring/correction spaces |
| `advertised-native-inputs.json` | Every explicit finite native tuple and exact claim occurrence, complete neighbours/corrections and all quantified charge domains |
| `deduplication.json` | Repeated assertion uses, evidence packets, supersession and explicit non-additive mathematical relationships |
| `reviews/` | Complete document reviews by the ledger, nonperturbative, F4, method and all-loop manuscript reviewers |
| `release-audit.json`, `release-sources.json.gz` | Immutable post-baseline release texts and complete semantic reviews |
| `BASELINE_GATE_RECEIPT.json` | Preserved passing baseline audit and eight deliberate-corruption controls |
| `STRICT_GATE_RECEIPT.json`, `CHECK_RECEIPT.json` | Current strict gate, release requirements and all nine deliberate-corruption controls |

The inventory includes only the public baseline and its reachable immutable
documentation links, together with Issue 1. Quarantined or isolated material
is neither read nor used. An identical document text at two commits shares a
textual review unit, but both immutable occurrences remain recorded: their
relative evidence links may resolve to different code. Identical prose does
not imply identical computational inputs.

The fixed review queue has 121 documents in the ledger lane, 79 in the
nonperturbative lane, 43 in the F4 lane, 11 in the method lane and the one
complete 5470-line all-loop manuscript. The complete diffs of historical
summary/discussion variants preserve their separate immutable contexts.
The accepted seven-letter SU3 certificate has the additional immutable pin
`eac7db42a6e27f4900826fe7b42cc75365698dc3`.

## Native and physical charges

Every finite native input stores the ordered five entries
`theta1, theta2, theta3, n_plus, n_minus` and its `letter_length`. It also stores
twice each physical charge, as exact integers, using the repository's
unchanged convention. The checker recomputes this translation.

The omega package originally orders the two spacetime Taylor degrees first.
Its tuple `(3,3,5,0,0)` is recorded alongside the flavour-first tuple
`(5,0,0,3,3)`. Both have length nine and physical level 28. No tuple is
silently reinterpreted in another package's order.

The finite core includes all eleven groups through level 36, all admissible
cochain lengths, and every relevant classical, quantum and correction
neighbour. The named-frontier extension additionally includes the fixed
orthogonal targets through level 84, the exceptional targets, both required
G2 level-36 source lengths, and the specified SU3 sources and targets.
All homogeneous earlier correction cochains remain present. The 22 compulsory
inputs contain 511 adjacent/correction slots and 148 kinematically allowed
maps, counting each named frontier separately. At the baseline, the additional
advertised input index contains 548 tuples, 1111 exact occurrences, 14488
adjacent/correction slots and 4317 maps; new release inputs are added explicitly
and the generated summary carries the updated counts. These counts describe
coordinates; they do not count scientific results. A mixed evidence packet's
group list is never multiplied into an invented group/charge census.

### Why the native list is finite and exhaustive as a kinematic list

Each letter can carry each of the three odd Taylor coordinates at most once.
Therefore each flavour degree is at most the letter count. Every relative
letter carries at least one Taylor coordinate, so the letter count is at
most the sum of the five native degrees. The vacuum is treated separately.

These conditions are also sufficient for an alphabet monomial. If the total
flavour degree is at least the requested letter count, place the flavour
occurrences successively around the letter positions; each flavour occupies
distinct positions because its degree is at most that count, and together
the occurrences cover all positions. If the total flavour degree is smaller,
put each occurrence in a separate position and place one dotted Taylor
coordinate in each remaining position. The upper-bound inequality guarantees
enough dotted coordinates. Distribute the remaining dotted degrees among the
positions. This constructs a nonempty allowed letter in every position.

This is an alphabet-existence argument. Gauge invariance, finite-colour
identities, statistics, highest-weight multiplicities and cohomology can make
many listed spaces vanish. Their dimensions remain uncomputed unless a
separate certificate establishes them. The 1,376,452 core native group slots
and 675,081 core physical group blocks are not counts of BPS states. They do
not define a discovery or lift percentage.

## Semantic records and evidence states

Each canonical record contains the group, charges, source, target, order,
field, normalization, scope, certificate, supersession state and explicit
limitations. A source or target is marked not applicable only when the record
does not assert a map. A uniform theorem instead records its quantified
domain and exact native parameter formulas.

The checker uses a canonical proposition key, group and coefficient field to
find incompatible active outcomes. Distinct source normalizations receive
their actual scope: the two SO11 chains are not two target discoveries. A
historical finite-field result cannot silently become a rational theorem.
An open physical matrix does not become zero because one classical profile
or selected topology vanishes.

`published_pending_independent_gate_review` reports a public proof claim
without asserting that its separate scientific acceptance gate has closed.
`attributed_public_claim_without_integrated_replay` identifies a published
assertion whose complete native packet has not been integrated. Neither
status is a new proof. The full F4 sixth-power incoming second-map zero is
separately accepted by the root audit: its complete published 3017-minor
lemma is an explicit prerequisite, while the new closed-factor proof was
checked directly. The all-page formal assertion has its own record.

## Commands

From the repository root:

```sh
python3 research-roadmap/claim-ledger/inventory.py
python3 research-roadmap/claim-ledger/coverage.py
python3 research-roadmap/claim-ledger/ledger_review.py
python3 research-roadmap/claim-ledger/build.py
python3 research-roadmap/claim-ledger/check.py --negative-controls --require-complete --receipt research-roadmap/claim-ledger/STRICT_GATE_RECEIPT.json
```

The ordinary check exits successfully when the ledger is internally
consistent. The strict check additionally requires every fixed semantic
completion condition. It exits with status 2 for incomplete semantic coverage,
and status 1 for structural corruption. The preserved baseline run exits
with status 0. The current strict receipt also requires every post-baseline
release review and publication pin.

The negative controls operate on copies in memory. They introduce a
contradictory active cohomology dimension, a missing immutable certificate,
a nonexistent document anchor, a duplicate claim identifier, a wrong
native-to-physical charge, a missing earlier correction slot, an omitted
semantic document review, an omitted advertised native occurrence and an
omitted required release. Each must be rejected. No original evidence is
modified by a negative control.

No script in this directory edits milestone weights, completion tests or
the roadmap score.
