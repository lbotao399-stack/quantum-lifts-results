# Source map

Use [the claim index](../claims.json) to determine a statement's current status. Use these sources to recover its conventions and calculation. The [source inventory](inventory.json) stores the origin, byte count, and SHA-256 for each pinned file. Historical uses of “pending”, “passed”, or “exact” refer to that source's own date and scope.

| Question | Open | What to take from it |
|---|---|---|
| Where did the program start? | [Claude's technical ledger](claude/ledger.tex) and [manuscript](claude/article.tex) | The original definition of the operators, proposed partners, and conjectures. For later completions or restrictions use the corresponding claim ID. |
| What defines the differential and quantum product? | [All-loop manuscript](all-loop/main.tex), [theorem guide](all-loop/theorems.md), [product-kernel specification](all-loop/product_kernels/README.md) | The common all-loop source revision is `3ef33773984a525d2520a4f1ffaf60f3b29ee04f`. The specification relates the saved tables through seven loops to the interaction and ordered contractions. |
| Where do the coefficients and signs enter code? | [Differential engine](engines/qengine.py), [product engine](engines/exact_star.py), [kernel implementation](all-loop/product_kernels/code/product_kernels.py) | Preserve the original field order, trace-dual convention, and loop factors. [Section 4](../research/04-complex-and-products.md) gives the sign change from the stored to the article's $C_4$. |
| What did the later calculations change? | [Research ledger](notes/claim_status.md), [duality audit](notes/duality_nonpert/research_status.md), [homotopy notes](notes/homotopy), [SQCD notes](notes/sqcd_duality_survey) | The targeted notes distinguish completed rational identities, effective models, and uncomputed local insertions. [The research chapter map](../research/README.md) links these questions to their full English derivations. |
| Where are the historical algorithms? | [Pure-theory products](engines/originals/rational_c4), [SU(4) search](engines/originals/su4_cartan), [instanton insertions](engines/originals/instanton/contact_matching), [SQCD scripts](engines/originals/sqcd_duality_survey) | Consult as algorithms. Some read sibling files or perform calculations and writes on import. The [portable replay catalog](../checks/catalog.json) lists the checks with an isolated entry point. |
| Which statements were taken from the literature? | [References by section](LITERATURE.md) and the two-adjoint source note below | Use the cited primary work for the physical input; use the local derivation for what follows from it. The citations are not an exhaustive novelty search. |

<a id="casimir-bracket"></a>

## Two trace invariants: the missing product-rule term

The [trace calculation and explicit primitive](notes/universal_casimir/trace_bracket.md) use the one-loop differential $Q_1$ on two traces of powers of the undifferentiated adjoint field $b$. Acting on their product produces a cross term. With powers $m,n\ge2$, the calculation expresses that term as $4mn$ times a tree-level variation. The proof differentiates two commuting powers of the same matrix; it applies before any finite-rank trace reduction. The [component check](engines/originals/universal_casimir/verify_symbol.py) inspects selected ranks and powers. Generality comes from the trace identity, not from those finite tests. Claim: `casimir-bracket`.

<a id="two-adjoint"></a>

## Two adjoints: what the duality evidence actually tests

- [Intriligator and Nardoni, arXiv:1604.04294](https://arxiv.org/abs/1604.04294) study flows of theories with two adjoint matrices. They describe the proposed truncations of the $D_{\mathrm{even}}$ and $E_7$ chiral rings as conjectural and give evidence against them. A proposed finite spectrum on the magnetic side does not supply an electric local primitive.
- [Bryan, Rajaraman and Shirman, arXiv:2512.04067](https://arxiv.org/abs/2512.04067) test nilpotent directions in the $A_k$ and $D_{k+2}$ models and report a failure of the proposed $D_{k+2}$ duality. Their analysis does not test $E_7$.

These are literature inputs to `two-adjoint-duality`. Neither paper is used here as evidence for a nonzero microscopic instanton differential. For the separate one-adjoint calculation and its specific unresolved meson terms, see [the cubic SU(5) comparison](../research/20-cubic-adjoint-su5.md).
