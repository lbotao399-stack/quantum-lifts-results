# F4 sixth-power proof audit and seventh-power source inventory

The complete incoming second differential to the sixth glueball is zero over
the rationals in the relative holomorphic jet complex. This packet audits the
published full-space theorem and supplies a direct closed-factor proof and
independent literal-tensor checks. It does **not** claim a fresh reconstruction
of the complete order-3017 minor.

The published exact bounds are

```text
3720 - 703 = 3017
3016 classical boundaries + 1 closed product = 3017 independent cycles
```

The complete minor is an attributed input at public commit
`c778cd7c3ceb72ad8082d791b2aee5f0df082164`, with residues `449` and `60`
at primes `1009` and `1013`. The rational outgoing rank `703` is the published
independent result at `3b30b8ff242977debd90581c12f4a767d7875c6c`.
See [the complete dependency argument](D03_PROOF.md),
[method and scope audit](METHOD_AUDIT.md), and
[immutable proof inputs](proof_dependencies.json).
The [public reproduction guide](PUBLICATION.md) gives the immutable input
download command and the exact permitted-file manifest.

## Fresh permitted calculations

| Calculation | Complete checked result |
|---|---|
| Formal product coefficients | 172 equations, 96 columns, rational rank 95, one primitive vector with 87 nonzero entries |
| Separate marked-factor closure | All 22 classical equations vanish, for the 45 active pairs of spectator factors |
| Entire three-slot invariant bilinear complex | 60 columns, classical rank 20, closed dimension 40; exact first-quantum row identity on all 60 columns |
| Native colour normalization | All 2704 entries of the adjoint Casimir equal 18 times the identity |
| Literal Albert extension | All 150850 entries of the 3017-by-50 bank rebuilt in bounded integer arithmetic |
| Outgoing prerequisite | Both saved order-703 determinants, exact characters and literal Cartan Gram complement checked; no native graph rebuilder invoked |
| Additional product Ward control | The complete rational order-8 determinant is `523032987223594133/125439321093750` |

The direct bilinear calculation proves that every closed marked factor is
strictly closed under the bare perturbative differential. The first-jet
spectator factors are strictly closed by flavour counting. The matched quantum
product chain identity then completes their actual polynomial product. This
argument uses no Taylor or colour subspace projector and does not require the
broader projected-product construction in the historical discussion.

The older local compiler used colour-subspace maps. That route was stopped and
its code and receipts removed from this publication directory. Its outputs are
not counted as permitted fresh verification. No default command below runs it.
The distinction is recorded in [stopped-route provenance](stopped_route_provenance.json).

## Replay

Python, NumPy, SciPy, python-flint and a C++17 compiler are required. Keep BLAS
single threaded and use at most two workers. The input materializer reads only
the 98 public numerical inputs and two additional public theorem/provenance
files listed in the immutable manifests.

```sh
python3 -B publication.py
python3 -B survival_inputs.py --download --destination /tmp/f4-public
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export F4_PUBLIC_SOURCE=/tmp/f4-public/c778cd7c3ceb72ad8082d791b2aee5f0df082164/algorithms/dual_singlets_20261118/f4/j36_source/focused_powers
g++ -O3 -std=c++17 -shared -fPIC albert_tensor_runner.cpp -o /tmp/f4-albert.so
python3 -B verify_closed_product.py --inputs "$F4_PUBLIC_SOURCE" --out /tmp/f4-closed-product.json
python3 -B verify_closed_bilinears.py --inputs /tmp/f4-public --out /tmp/f4-closed-bilinears.json
python3 -B albert_tensor_runner.py --inputs "$F4_PUBLIC_SOURCE" --library /tmp/f4-albert.so --points 12 32 --out /tmp/f4-controls.json
python3 -B rebuild_albert_extension.py --inputs "$F4_PUBLIC_SOURCE" --library /tmp/f4-albert.so --workers 2 --out /tmp/f4-literal-extension
python3 -B dependency_replay.py outgoing --inputs /tmp/f4-public --out /tmp/f4-outgoing.json
python3 -B verify_packet.py --inputs /tmp/f4-public --out /tmp/f4-packet-check.json
```

The extension command checks literal invariant contractions, not an action on a
selected subspace. It makes no assertion about the determinant of the complete
order-3017 matrix. The retained literal extension is
`literal_albert_extension_1009.u16.gz`, with ordered points, shape and hashes in
`literal_albert_extension_1009.json`.

## Exact native inventories and correction freedom

Charges below use two dotted degrees followed by three flavour degrees.

| Target | Incoming second-source charges | Sole allowed Taylor type | Lengths | Dimensions |
|---|---|---|---|---|
| Sixth power | `(4,4,2,2,2)` | `(4,4,3,3)` | `12,13,14` | `10007,3720,704` |
| Seventh power | `(5,5,2,2,2)` | `(5,5,3,3)` | `14,15,16` | `217297,62866,8553` |

All ten physical spin/flavour singlet types are listed in
`second_source_inventory.json`; their respective full singlet dimensions are
`228269,70238,10160` and `2558557,663400,79101`. The full Taylor Ward identity
excludes maps from the other types. No selected polynomial ansatz replaces that
identity.

| Target | First-correction charges | Lengths | Complete singlet dimensions | All homogeneous first corrections |
|---|---|---|---|---:|
| Sixth power | `(5,5,1,1,1)` | `11,12,13` | `13570,4109,633` | `4109-633=3476` |
| Seventh power | `(6,6,1,1,1)` | `13,14,15` | `93168,23895,2959` | `23895-2959=20936` |

The character code is byte-for-byte pinned to the independent public character
package. Rebuild either inventory with:

```sh
python3 -B character_inventory.py --out /tmp/f4-second-source.json
python3 -B character_inventory.py --order 1 --out /tmp/f4-first-corrections.json
```

## Seventh power and the remaining native-input gate

The seventh-power first-page Ward identity was reconstructed completely over
the rationals: 2704 scalar pairs, 1378 symmetric pairs, 1117 differentiated
scalar tests and 18948 mixed-ghost components. The first-correction map has
full terminal rank 2959, with the order-68 determinant residues
`909,962,176568` at primes `1009,1013,1000003`. These are prerequisites for its
second-map calculation, which remains open. See [D04 status](D04_STATUS.md).

```sh
python3 -B dependency_replay.py seventh-first-page --inputs /tmp/f4-public --out /tmp/f4-seventh-first-page.json
python3 -B dependency_replay.py seventh-intermediate --inputs /tmp/f4-public --out /tmp/f4-seventh-intermediate.json
python3 -B dependency_replay.py seventh-colour --inputs /tmp/f4-public --out /tmp/f4-seventh-colour.json
python3 -B seventh_ward_inventory.py --inputs "$F4_PUBLIC_SOURCE" --out /tmp/f4-seventh-ward-inventory.json
```

C04 remains open. Neither the published full-kernel bound nor the literal
extension supplies a complete native basis of all 10007 predecessors and its
full action kernel. No performance milestone is claimed from pointwise timings
or the reused checkpoint rows. No physical instanton statement follows from
these formal perturbative calculations.
