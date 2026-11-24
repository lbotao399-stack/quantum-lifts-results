# Semantic review exchange

The finite document queue is `review-queue.json`. A document is one distinct
SHA-256 of public text; `occurrences` preserves every immutable occurrence.
The source bytes remain in `inventory/public_sources.json.gz`.

Review outputs use the following JSON structure:

```json
{
  "schema_version": 1,
  "reviewer": "agent or reviewer identifier",
  "documents": [
    {
      "sha256": "full source text hash from the queue",
      "classification": "in_scope | mixed | out_of_scope | navigation",
      "review_complete": true,
      "claim_ids": ["existing or new canonical identifiers"],
      "covered_assertions": [
        {
          "location": "heading, source anchor, exact row or short literal locator",
          "claim_ids": ["canonical identifier"],
          "description": "The exact assertion or group of claims governed by one stated theorem."
        }
      ],
      "scope_exclusions": [
        {
          "location": "heading or exact source scope",
          "reason": "The actual reason this material is not an undeformed-N4 result."
        }
      ],
      "metadata_gaps": [],
      "notes": "Any limitation of this semantic review."
    }
  ],
  "additional_claims": []
}
```

`additional_claims` uses the complete schema of `claims.json`: identifier,
theory, group, title, canonical proposition key/value, all charges, source,
target, order, coefficient field, normalization, exact scope, immutable
certificate, supersession state, proof status, superseded identifiers,
limitations, unresolved metadata and advertised locations. Native input
references can point to a detailed immutable proof or input table that
actually states those data. A generic catalogue link is insufficient.

Classify the complete document. A paragraph archive is not a semantic review.
An in-scope document can support several already-recorded assertions without
creating duplicate scientific results. Explain excluded method instructions,
literature background, other theories and historical scheduling details.
Do not call a document complete while a result-bearing section remains
unreviewed. Identical text shares textual classification; the immutable
evidence links at each occurrence remain checked separately.

Use historical attributed status for old claims whose proof payload is not
available. Record the exact stated field and scope; an absent fresh replay
does not by itself make the semantic record incomplete. Actual unstated
charges, an ambiguous active normalization, or an unlinked advertised result
remain explicit metadata gaps. Do not upgrade an advertised claim to a
scientifically accepted theorem during this bookkeeping review.

The nonperturbative review owner writes `reviews/nonperturbative.json` and
handles every queue item with `review_lane` equal to `nonperturbative`.
The ledger owner handles every other item. Reviewers must not read isolated
or quarantined paths; the queue contains only permitted public documents.
