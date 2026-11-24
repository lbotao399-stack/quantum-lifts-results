# Exact claim-ledger transport

The large reviewed JSON and gzip artifacts are stored losslessly in bounded
XZ chunks. `ARTIFACTS.json` records every original file's byte count, offset
and SHA-256, the whole uncompressed concatenation and archive hashes, and
every chunk hash. All mathematical and semantic receipt hashes remain
unchanged.

From the repository root, Python 3 alone restores the original paths:

```sh
python3 -B research-roadmap/claim-ledger/assemble_artifacts.py
```

The command authenticates every byte before writing. It rejects a corrupted
existing decoded artifact instead of silently replacing it. Decoded files are
ignored by Git; code, review schemas, current human-readable summary, release
reviews and strict acceptance receipts remain directly readable in the tree.
Transient inventory checkpoints and reviewer transfer files are unnecessary
for the completed ledger and are not published.

The complete semantic/reference audit is:

```sh
python3 -B research-roadmap/claim-ledger/check.py --negative-controls --require-complete --receipt /tmp/claim-ledger-check.json
```

This final command needs the public Git objects at every immutable reference
in the ledger. Its published strict receipt was produced with those objects
present and all references checked. A shallow checkout does not automatically
contain discussion-linked commits. The bounded continuous-integration tier
authenticates this complete-audit receipt and every named artifact; it does
not label that operation a new semantic review or a new Git-reference audit.
