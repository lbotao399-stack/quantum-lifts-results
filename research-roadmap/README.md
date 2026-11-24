# Quantum-lifts research programme: 100-point roadmap

**Current stage: 23/100 points, with 19 of 60 milestones met.**

Assessment date: **2026-11-24, Asia/Shanghai**. Evidence is frozen at public commit
[`2f268a0fb47a7b17119cda86a95f2fac574648cf`](https://github.com/lbotao399-stack/quantum-lifts-results/commit/2f268a0fb47a7b17119cda86a95f2fac574648cf)
and the immutable public inputs linked from it. The user requested and authorized this roadmap and its publication.

The current programme has exact finite-sector fortuity certificates, uniform first-order results, independently checked orthogonal third-order lifts, and substantial all-order **formal** survival theorems. Autonomous complete discovery, general higher-order coverage, microscopic physical instanton maps, and a physical BPS-chaos dataset remain unfinished. The score counts completed deliverables under the explicit rubric below; it does not measure the fraction of an infinite classification problem already solved.

## Read the roadmap

| Document | Purpose |
|---|---|
| [SCORECARD.md](SCORECARD.md) | Exact score, all 60 milestones, and execution queue |
| [MILESTONES.md](MILESTONES.md) | Every completion test, dependency, current finding, next calculation, and required artifact |
| [CURRENT_STAGE.md](CURRENT_STAGE.md) | What the current results establish; structural bottlenecks and the first attacks |
| [SCOPE.md](SCOPE.md) | The finite coverage contract, uniform-family goals, and physical/statistical boundaries |
| [EVIDENCE.md](EVIDENCE.md) | Immutable sources and their precise limits |
| [milestones.json](milestones.json) | Machine-readable milestone and score source |
| [scope.json](scope.json) | Machine-readable coverage specification; the expanded block manifest is still to be generated |
| [evidence.json](evidence.json) | Machine-readable provenance register |
| [AGENTS.md](AGENTS.md) | How subsequent research selects, tests, closes and revises a milestone |
| [HISTORY.md](HISTORY.md) | Version history and the rules for changing this programme |

## Point allocation

| Area | Current | Maximum |
|---|---:|---:|
| Scientific reliability | 6 | 8 |
| Fortuitous discovery and classification | 4 | 16 |
| Transformative algorithms and autonomous execution | 2 | 18 |
| Perturbative higher differentials | 8 | 20 |
| Physical nonperturbative lifts | 2 | 14 |
| Physical explanation and dual descriptions | 1 | 14 |
| BPS chaos | 0 | 10 |
| **Total** | **23** | **100** |

The weights express research priorities, not a law of nature or an estimate of labour. Thirty-four points target discovery and algorithms because those improvements determine whether the rest of the programme can scale. The physical interpretation and chaos tracks require physical operators and observables, not just larger cochain matrices. Dependencies carry no additional points.

## Scoring discipline

- A milestone earns its entire integer weight only when every acceptance condition holds. Partial and open milestones earn zero. The initial score audits published certificates, with their stated replay limits; it does not claim every large computation was rerun during this assessment.
- A new document, a proposed framework, a selected ansatz, an index deficit, a modular lower bound, a fast kernel, or a formal existence theorem cannot receive credit for a stronger deliverable it does not establish.
- The same certificate may support distinct capabilities, but repeating a target, rescaling a source or multiplying by a known class creates no extra result count. Targets, multiplets, map ranks and independent generator families remain separate units.
- Complete negative results count when the milestone asks for a decision. They do not become positive lifts or positive chaos findings. An inconclusive statistical result cannot close a scaling or mechanism gate.
- Deformed, Abelian and lower-supersymmetry examples are controls. They earn no positive undeformed-N4 lift credit. The historical Konishi control earns one normalization point, with its long-multiplet scope explicit.
- No point is awarded for writing this roadmap. In particular, A07 remains partial because a complete repository-wide claim and coverage ledger is still absent.

## Verify and update

From the repository root, with Python 3 and no third-party packages:

```sh
python3 research-roadmap/score.py
python3 research-roadmap/score.py --ready
```

After a justified change to `milestones.json` or `evidence.json`:

```sh
python3 research-roadmap/score.py --write
python3 research-roadmap/score.py
```

The checker rejects invalid totals, duplicate identifiers, missing evidence, cyclic dependencies, met milestones with unmet prerequisites, and stale generated documents. It validates the ledger; the milestone's mathematical or physical certificate remains the substantive proof.

## Immediate research direction

The first scientific target is **D03: decide the full (F_4) sixth-power second differential**. C03 and C04 supply the reusable complete input and verification machinery; an exact analytic full-space argument may also resolve D03 without building a dense matrix. B05 supplies a separate concrete route to the missing (SU(3)) generator and its physical channel. F02 and F03 are the necessary physical work before an instanton or chaos claim can be tested.

The roadmap is revisable. Changing a question, scope, weight or completion test requires a versioned reason and an old/new score comparison. A newly discovered obstruction can change the research question; it cannot silently lower an existing completion test.
