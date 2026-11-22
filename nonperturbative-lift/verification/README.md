# What was verified

Read the `scope` of each entry in [the calculation report](replays.json). A `passed` entry is about the specified equation or algorithm. `proved` in [the claim index](../claims.json) may instead rely on an analytic proof, a historical exact certificate, or an implication from a cited input; it is not a count of freshly rerun theorems.

| Kind of evidence | Where it lives | What it establishes |
|---|---|---|
| Project replay | [Replay report](replays.json) · [commands and scopes](../checks/catalog.json) | The listed scripts ran in isolated temporary directories with assertions enabled. The report records the environment; each entry gives a script digest, output, and exit status. The exact interpreter and package versions are recorded in the run report. |
| Historical exact calculation or independent readback | [Expanded evidence](../evidence) · [certificate ZIP](../artifacts/nonperturbative-certificates.zip) · [certificate table](../research/26-certificates.md) | Exact coefficients, nonvanishing detectors, primitives, rank arguments, and previously saved residuals in their stated sectors. Opening a saved zero-residual report is different from reconstructing the polynomial equation today. |
| Hash and structural integrity | [Integrity report](integrity.json) · `python3 -B tools/project.py verify` · [source manifest](../sources/inventory.json) · [archive manifest](../artifacts/inventory.json) | Rehashes pinned originals, the ZIP and each member; compares the expanded records with the archived bytes. Verifies the digest of each split chapter, its statistics, claim dependencies, continuation targets, and local Markdown links. It neither computes a differential nor performs a new peer review. |
| Literature input | [Primary references](../sources/LITERATURE.md) · [two-adjoint source note](../sources/README.md#two-adjoint) | The physical theorem, superpotential, or duality premise is attributed to its source. Local algebra tests consequences under the specified hypotheses. |

The [massless source](../checks/massless_source.py) builds the full tree matrix and four signed detector rows. The [scalar check](../checks/scalar_saturation.py) eliminates the inverse weight, expands the Pfaffian, and verifies the local second-order primitive. The [translation check](../checks/translation_jets.py) acts on exterior jet monomials and recomputes conformal descendant coefficients. The proofs giving all ranks or orders are in the corresponding chapters.

**Example: the SU(2) nilpotency certificate.** The [boundary replay](../checks/su2_boundary.py) reconstructs the product $C_0\star C_2$ and the differential of the stored primitive from the bosonic and ghost letters. The three coefficients in loop order have product sizes `6, 75, 135`, primitive sizes `24, 221, 105`, and residual sizes `0, 0, 0`. Higher coefficients vanish by bosonic field count. The separate [hidden-square primitive record](../evidence/rational_c4/full_completion_summary.json) and [nonzero detector record](../evidence/rational_c4/nonzero_summary.json) are historical certificates. The [nilpotency proof](../research/04-complex-and-products.md#np-5c85ee5e8a2c3293-1) uses all three ingredients and associativity; the fresh boundary computation checks one of the three.

From the project root:

```sh
python3 -B tools/project.py replay --list
python3 -B tools/project.py replay su2-boundary
python3 -B tools/project.py verify
```

The essential command, reading a claim, needs only Python: `python3 -B tools/project.py show sqcd-microscopic-image`. For a binary witness, run `python3 -B tools/project.py extract rational_c4/C0_C2_boundary_exact.pkl /tmp/np-lift-boundary.pkl` with a **new** destination. The tool checks the ZIP and member digest before writing. The package setup for symbolic replays is in the [project entry point](../README.md). Standalone checks never require the old reader or original checkout. Historical continuation scripts are different: some execute on import and require expanded polynomials or sector-specific matrix columns.

<a id="external-witnesses"></a>

## Published polynomial data

The original certificate ZIP contains 78 payloads. Its historical manifest lists four further coefficient files. All four are available individually. The [publication manifest](../artifacts/publication.json) records their original paths, compressed and expanded byte counts, and both SHA-256 digests.

| Polynomial | Expanded bytes | Calculation record |
|---|---:|---|
| [SU(2) third correction](https://github.com/lbotao399-stack/quantum-lifts-results/releases/download/nonperturbative-lift-v1/su2-third-correction.pkl.gz) | 686635425 | [934381 terms and original SHA-256](../evidence/rational_c4/m6_lift/compact_first/lift/level3_exact_summary.json) · [independent calculation](../evidence/rational_c4/m6_lift/compact_first/lift/level3_independent_readback.json) |
| [SU(3) product](../artifacts/witnesses/su3-c0-product.pkl.gz) | 21971105 | [278026 product terms](../evidence/su3_rational/c0_product_Q_certificate.json) · [rank argument](../evidence/su3_rational/proofs.md) |
| [SU(4) degree-four first correction](../artifacts/witnesses/su4-k4-first-lift.pkl.gz) | 12005798 | [7712 leading terms, 81401 primitive terms](../evidence/su4_cartan/k4_first/first_lift.json) · [independent calculation](../evidence/su4_cartan/k4_first/independent_readback.json) |
| [SU(4) degree-five first correction](../artifacts/witnesses/su4-k5-first-lift.pkl.gz) | 59680646 | [19628 leading terms, 439142 primitive terms](../evidence/su4_cartan/k5_first/first_lift.json) · [independent calculation](../evidence/su4_cartan/k5_first/independent_readback.json) |

The recovered SU(2) bytes reproduce the SHA-256 recorded before publication. For the other three files, decoding the coefficient dictionaries reproduced the term counts in the table. Their first published SHA-256 digests identify those files in subsequent downloads. The differential residuals in the linked records come from the original calculations.

The [historical calculation ZIP](https://github.com/lbotao399-stack/quantum-lifts-results/releases/download/nonperturbative-lift-v1/historical-calculations.zip) contains the early engine and checks, continuation programs, coefficient data, research notes, and argument diagrams. The archive includes a copy of its [file index](../artifacts/historical-files.json), which lists each research file by path, size, and SHA-256. Numerical solver caches, installed environments, and third-party paper copies sit outside this collection. The archived manuscript predates the chapter updates. Its source reproduced the historical manuscript digest, and the public copy redirects 27 reader links to their corresponding theorem or derivation anchors in this edition. The 17 registered replays use the files in the Git checkout.

To verify a complete download, run these commands from the project directory. The two downloads can also be saved by following the links above.

```sh
release_dir=$(mktemp -d)
gh release download nonperturbative-lift-v1 --repo lbotao399-stack/quantum-lifts-results \
  --pattern 'su2-third-correction.pkl.gz' --pattern 'historical-calculations.zip' --dir "$release_dir"
python3 -B tools/project.py verify --release-dir "$release_dir"
```

The verifier streams each gzip file to check its expanded coefficients byte for byte. It also checks the historical ZIP and every indexed member. Running `verify` without the release directory checks the three polynomials stored in Git and reports the number checked. To recover an individual pickle, use `gzip -dc input.pkl.gz > output.pkl` in a fresh directory.
