# Discovery implementations

These are the actual stable-complex, identity-generation, relation-quotient,
and rational-reconstruction implementations used in the research run.
Only repository loader prefixes were changed to locate the included
`upstream/` snapshot. `SOURCE_PROVENANCE.json` records both original and
delivered hashes. No mathematical algorithm or normalization was changed.

`jets/pfaffian_relations.py` expands antisymmetrized symplectic index
contractions with Grassmann signs. `paired_relation_basis.py` and
`path_relation_completion.py` construct the complete small highest-weight
relation spaces. `relation_quotient.py` and `sp5_relation_quotient.py` reduce
target identities against the actual universal Q0 image in sparse blocks,
extract rational preimages, and check closure coefficient by coefficient.

The stable basis/map builders are `jets/formal_frontier.py` for Sp(4) and
`sp5_universal/stable_complex.py` for Sp(5). Their checkpoints can be
regenerated from the included upstream Python code. Large search
checkpoints and intermediate modular matrices are not part of the proof
delivery; the exact finite identities, operators, and separating duals are.

For a complete independent proof replay, use the top-level `replay.py`.
These research search scripts operate on their own working-directory
checkpoints and are supplied as implementations rather than a polished
one-command campaign interface. The frozen recipes fix every identity
used by the published certificates, independently of random search order.
