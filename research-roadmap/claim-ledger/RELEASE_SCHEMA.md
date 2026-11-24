# Finite post-baseline release review

The 255-text baseline queue is frozen. Three finite mathematical release
packets supplement it: `b05`, `physical`, and `f4`. The root reviewer owns
`releases/*.json`; the builder and checker authenticate their exact immutable
text bytes. Generated ledger, navigation and score reports are excluded from
recursive self-review.

Required source roots:

| Release identifier | Complete text root |
|---|---|
| `b05` | `algorithms/su3_generator_20261124` |
| `physical` | `algorithms/physical_block_20261124` |
| `f4` | `algorithms/f4_complete_source_20261124` |

Every `.md` and `.tex` file below the corresponding immutable root must occur
exactly once in its packet. Packaging and publication notes can be classified
as navigation, with an explicit scope reason. This finite set is checked
against the Git tree; a missing mathematical text fails the check.

```json
{
  "schema_version": 1,
  "release_id": "physical",
  "commit": "immutable 40-character Git commit",
  "reviewer": "root acceptance reviewer",
  "documents": [
    {
      "path": "algorithms/physical_block_20261124/README.md",
      "sha256": "SHA-256 of exact UTF-8 Git blob bytes",
      "classification": "in_scope",
      "review_complete": true,
      "claim_ids": ["the exact canonical identifiers"],
      "covered_assertions": [
        {
          "location": "literal heading or exact source location",
          "claim_ids": ["the exact canonical identifiers"],
          "description": "Complete assertion, exact normalization, acceptance scope and remaining limits."
        }
      ],
      "scope_exclusions": [],
      "metadata_gaps": [],
      "notes": "Review and independent-check scope."
    }
  ],
  "additional_claims": [],
  "claim_updates": []
}
```

`additional_claims` uses every required field from `claims.json`, including
exact charges, source, target, coefficient field, normalization, complete
scope, immutable certificate references and supersession state. F02 and F03
remain separate claims with their actually accepted coupling orders and
physical channels. Unaccepted later physical-survival work is excluded.

Existing B05 and D03 canonical identifiers can receive this limited update:

```json
{
  "id": "f4-power6-full-second-zero-report",
  "append_certificate": [
    {
      "kind": "git",
      "commit": "immutable 40-character Git commit",
      "path": "algorithms/f4_complete_source_20261124/D03_PROOF.md",
      "url": "immutable public source URL"
    }
  ],
  "proof_status": "accepted_published_theorem_with_independent_closed_factor_audit",
  "append_limitations": ["The published full minor remains an attributed prerequisite; no new full-minor reconstruction is claimed."],
  "acceptance_note": "Exact accepted scope and independent decisive check."
}
```

Only `id`, `append_certificate`, `proof_status`, `append_limitations` and
`acceptance_note` are permitted in an update. A new scientific assertion
requires a complete additional record. Publication or packaging alone does
not change mathematical acceptance.

The builder generates `release-audit.json` and `release-sources.json.gz`, then
joins the release records into the current ledger and summary. The checker
requires all three releases, verifies the exact tree text list, every text
hash, semantic assertion mapping and immutable claim reference. Its ninth
negative control removes a required post-baseline release and requires
rejection. The baseline receipt is preserved in `BASELINE_GATE_RECEIPT.json`.
