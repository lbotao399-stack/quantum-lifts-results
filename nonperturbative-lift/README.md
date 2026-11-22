# non perturbative lift

This project asks which supersymmetric local operators survive perturbation theory but become a supersymmetry variation after instantons are included. A surviving operator modulo such variations is a **cohomology class**. The distinction that matters here is *which differential acts on which operators, over which coefficients*. The [comparison by physical effect](research/23-operator-data.md) separates changes to classes, products, vacuum equations, and interactions.

**Start with one claim.** Search the [claim index](claims.json); follow its `derivation` to the exact passage, then `evidence` to the proof or witness. Use the [results chapter](research/02-results.md) for the formulas, the [review](research/03-review.md) for first-use definitions, and [open problems](OPEN_PROBLEMS.md) for the equations that remain. No reader, server, or external checkout is needed to navigate the project.

| Research question | Start here | What is established; what remains |
|---|---|---|
| What survives in pure supersymmetric Yang–Mills theory? | `su2-nilpotency` · [products](research/04-complex-and-products.md) · [finite-rank classes](research/05-leading-modules.md) · [SU(3) lifts](research/06-su3-lifts.md) | The perturbative $SU(2)$ class $[C_2]$ has nonzero square and zero cube without inverting the loop parameter. $SU(3)$ has complete representatives for $C_3,C_4$. The $SU(4)$ degree-six nonabelian primitive is open. |
| What would turn a pure-theory survivor into an instanton pair? | `secondary-instanton` · [operator identity](research/09-secondary-pairing.md) · [physical matching](research/13-connection-matching.md) | The algebraic implication has explicit hypotheses. The instanton connection contribution has been integrated; the complete $SU(2)$ local coefficient has not. |
| Which SQCD operators disappear? | [Full determinant ideal](research/16-determinant-lift.md#determinant-full-torsion) · [same confined theory in Pfaffian coordinates](research/18-sp-and-so.md#pfaffian-full-torsion) · [fermionic targets](research/18-sp-and-so.md#pfaffian-jet-images) | The determinant generates the entire scalar kernel in magnetic quark coordinates. The Pfaffian coordinates give fifteen quadratic generators, and their ungauged effective model gives 105 independent first fermionic images at six flavors. |
| Does the microscopic magnetic source survive? | [All-power massless proof](research/16-determinant-lift.md#massless-all-powers) · [necessary mass pole](research/16-determinant-lift.md#mass-pole-ward) | Sixteen independent current components span a free loop-parameter module. The incoming class detects a [failure of joint mass and loop specialization](research/16-determinant-lift.md#mass-loop-base-change). The proposed massive microscopic image requires a specific inverse-mass residue on an exact family of tree vacua. Its local instanton matrix element is open. |
| What decides the instanton order? | [Saturation and invariant factors](research/23-operator-data.md#complete-tests) · [circle and second-page example](research/19-higher-f-terms-and-monopoles.md#multiplier-pages) | In a multiplier model the parameter acts by the defining polynomial in its Jacobian algebra. A superpotential linear in the parameter can have a second-page scalar target. |
| Which derivative operators are protected? | [N=2 bulk translation](research/21-seiberg-witten.md#bulk-translation) · [N=4 descendant count](research/22-n4-sym.md#n4-descendant-detector) | The complete basic bilinear source excludes instanton pairings of two towers in the stated N=2 bulk complex. N=4 half-BPS total translations have an independent nondegenerate two-point detector. |
| When is an ordinary superpotential enough? | `ads-global-unit` · [nonpolynomial potentials](research/17-nonpolynomial-superpotentials.md) · [symplectic and orthogonal duals](research/18-sp-and-so.md) · [higher F-terms and monopoles](research/19-higher-f-terms-and-monopoles.md) | Eliminating fields can produce a pole or several branches. The calculations specify the allowed denominators and the equations before elimination. Interactions with additional fermions and magnetic-flux insertions have separate tests. |
| What do other dualities and extended supersymmetry constrain? | `su5-kss-matching` · [cubic adjoint](research/20-cubic-adjoint-su5.md) · [two-adjoint sources](sources/README.md#two-adjoint) · [Seiberg–Witten periods](research/21-seiberg-witten.md) · [N=4 SYM](research/22-n4-sym.md) | Ordinary vacua do not detect nilpotent operators. Exact diagonal vacua do detect the specified pure-scalar classes in $\mathcal N=4$ SYM. No instanton pairing of a complete nonscalar perturbative survivor is certified here. |
| How were the identities obtained? | `forest-normalization` · `period-hessian` · [residues](research/07-twisted-periods.md) · [exact arithmetic](research/08-exact-arithmetic.md) · [source map](sources/README.md) | Product kernels fix coefficients and signs. Rational primitives, rank bounds, and evaluations detect different things; the records state which one proves each claim. |

The [polynomial inventory and verification instructions](verification/README.md#external-witnesses) link all four expanded polynomials omitted from the original certificate ZIP. The [versioned release](https://github.com/lbotao399-stack/quantum-lifts-results/releases/tag/nonperturbative-lift-v1) also contains the earlier calculations and argument diagrams, searchable by path in the [historical file index](artifacts/historical-files.json).

**Try it from the project root:**

```sh
python3 -B tools/project.py claims --status open
python3 -B tools/project.py claims determinant
python3 -B tools/project.py show su2-nilpotency --json
python3 -B tools/project.py show sqcd-microscopic-image --next
python3 -B tools/project.py verify
python3 -B tools/project.py replay su2-boundary
```

These commands use the Python standard library. For a script invoked by absolute path the working directory can be anywhere. To run the full [replay catalog](checks/catalog.json), use a separate environment; the recorded run used CPython 3.13.2 and the packages pinned in [requirements.txt](requirements.txt):

```sh
python3 -m venv /tmp/np-lift-venv
/tmp/np-lift-venv/bin/python -m pip install -r requirements.txt
/tmp/np-lift-venv/bin/python -B tools/project.py replay --list
/tmp/np-lift-venv/bin/python -B tools/project.py replay all --report /tmp/np-lift-replays.json
```

[Agent instructions](AGENTS.md) · [Chapter index](research/README.md) · [Verification and downloads](verification/README.md) · [Certificate inventory](artifacts/inventory.json) · [Integrity report](verification/integrity.json) · [Replay results](verification/replays.json). All paths in the claim index are relative to this directory; `archive:` denotes a member of the certificate ZIP. The index records `proved`, `literature-input`, `conditional`, `finite-evidence`, `conjecture`, and `open` separately. A replay records only what its `scope` says it tested.
