# Sp(6): exact classical fortuity and a reusable inverse-Q calculation

This addendum extends the four-certificate package published at commit
`2dc1dec8eee723c44e9da10ac71b7648d5c20fe1`. It supplies an explicit
Sp(6)=USp(12) highest-weight cochain in paper sector
`(11,0,0;3,3), Y=15`, with `j=40` and native weight `(3,3,11,0,0)`.

**The fortuitous multiplicity in this sector is exactly one over Q.**
The full finite-sector cohomology dimension is not inferred. This is
a new certified member relative to the fixed preceding calculation;
independence from the existing symplectic rank family is not asserted.
All statements concern classical `Q0`. Loop corrections have not been
computed for this cochain.

## Exact result

| Quantity | Value |
|---|---:|
| Highest-weight operator trace terms | 152833 |
| Least trace-coefficient denominator | 3 |
| Universal `Q0`-defect trace terms | 98882 |
| Complete relation dimensions at lengths 14, 15, 16 | 67, 188, 183 |
| Rank of `Q0` on length-15 identities | 121 |
| Rank of `Q0` on length-16 identities | 45 |
| Stable abelianization rank on length-16 identities | 24 |
| Combined `Q0` and abelianization rank | 61 |
| Exact combined minor | `-2^109 * 3^21 * 5^7` |
| Nonzero coefficients of the connected-trace dual | 111 |
| Source-identity images annihilated by the dual | 188 |
| Dual pairing with the supplied defect | 1 |
| Fortuitous dimension | `(183 - 61) - 121 = 1` |

The cochain satisfies all three native highest-weight equations.
Its `Q0` variation equals the supplied integer defect exactly, and the
defect is an explicitly recorded combination of genuine finite-rank
Pfaffian identities. The rational dual annihilates `Q0` of every
identity in the complete source frame and evaluates to one on the
defect. It therefore proves non-exactness and exclusion from the
stable cohomology image without an inner product or graviton
orthogonality prescription.

The 61-by-61 integer minor supplies the upper bound. The 121 independent
relation boundaries and the explicit nonboundary defect supply the
matching lower bound. The exact character dimensions are:

| Field length | Universal highest weights | Sp(6) highest weights | Identities |
|---|---:|---:|---:|
| 14 | 107126 | 107059 | 67 |
| 15 | 100776 | 100588 | 188 |
| 16 | 47346 | 47163 | 183 |

## What changed in the algorithm

1. **Cycle expansion.** The new Pfaffian implementation groups endpoint
   matchings into trace cycles and handles repeated scalar-pair paths
   by multiplicity. Its integer coefficients were compared directly
   with the original endpoint implementation, including three Sp(6)
   inputs with all `13!! = 135135` matchings.
2. **Defect-specific inversion.** The program first identifies the one
   nonboundary direction in the small relation complex, then solves
   only for a universal preimage of that defect. Its 1961 trace-weight
   blocks all have exact rational solutions. A complete universal
   cochain census was not used to construct this witness.
3. **Reusable scalar spectators.** After stripping `Q0`-closed
   pure-scalar trace factors, 1374 of the Sp(6) inverse tasks match
   solved Sp(5) tasks exactly, with their recorded integer factors.
   Only 587 tasks are new in this comparison. These are counts of
   linear inverse problems, not independent physical families.

The largest inverse block used 17718 highest-weight source expressions
and an image of rank 7627. Every reconstructed rational answer was
checked with the native differential. A failed source ansatz is never
used as a non-existence certificate.

## General exclusion and the adjacent rank

[FIRST_LAYER_EXCLUSION.md](FIRST_LAYER_EXCLUSION.md) proves that stable
abelianization is injective on the first symplectic trace-identity
layer, at field length `2r+2`. Consequently finite-rank `Q0` cohomology
coincides with universal cohomology through length `2r+1`.

For the present length-15 class, the adjacent Sp(7) sector is therefore
outside fortuity. The independent explicit check has four complete
length-16 identities, abelianization rank four and determinant eight.
It proves fortuitous dimension zero there. The Sp(6) class cannot lift
to a Sp(7) class under the trace covering.

The cited invariant-theory input is the symplectic second fundamental
theorem. This exclusion argument is not a claim to have discovered a
new general invariant-theory theorem.

## Files

- `data/sp6_fortuitous_operator.json.gz.part-000` and `.part-001`:
  the two byte segments of the original verified gzip, containing every
  rational coefficient of the highest-weight operator in native trace
  notation. `data/sp6_operator_parts.json` records both segment hashes
  and the original whole-file hash. This storage split avoids a stalled
  single large GitHub upload; it does not change the mathematical data.
- `data/sp6_defect.json.gz`: its integer universal `Q0` variation.
- `certificates/sp6_fortuitous_certificate.json`: final exact results,
  hashes, character counts and scope.
- `certificates/portable_verification.json`: successful exact replay of
  the frozen cochain and separating certificate using the packaged
  pinned upstream source.
- `certificates/sp6_small_complex.json.gz`: obstruction-screen
  snapshot, including the complete defect coefficients, exact dual
  and integer-minor probes. Its preimage-pending flags describe that
  earlier stage; the final certificate supplies the completed witness.
- `certificates/sp6_rational_coordinates.json.gz`: exact small-matrix
  coordinates, including all first-image inclusions.
- `certificates/sp6_support_solver_result.json.gz`: all 1961 block
  verification records.
- `CYCLE_EXPANSION.md`: cycle-sign derivation, validations, and the
  unsuccessful one-/two-probe simplification at Sp(5).
- `discovery_scripts/`: additional experiments and coefficient-comparison
  programs. The portable entry points below use the parent package's
  pinned upstream source; no original scratch checkout is needed.

The construction receipts also retain hashes of historical pickle
caches. Those hashes record the discovery run. Verification checks
exact polynomials and equations, independent of Python pickle's
object-sharing layout.

## Reproduce

Keep this folder as `sp6_extension/` inside the complete parent package.
Use the parent package's documented Python dependencies.

The verifier reads and hashes the two operator parts automatically.
To obtain the original single gzip file separately:

```bash
python sp6_extension/assemble_operator.py --output /tmp/sp6_fortuitous_operator.json.gz
```

Its SHA256 is
`7b1c026ab478dea9f20e51c3f3c14086df52ab547049200015a62dcef9c6a7a5`,
identical to the file used in the completed exact mathematical replay.

```bash
python sp6_extension/run.py reproduce --workdir /tmp/fortuity-sp6
```

This rebuilds the 438 identities, recomputes the small-complex proof,
finds the universal preimage, projects the cochain to highest weight,
verifies the frozen operator, and checks the adjacent-rank exclusion.
Long stages checkpoint and resume. The default source cap for a single
inverse block is 30000, larger than the 17718 expressions used by the
largest block in this run.

To verify the supplied cochain using an existing reconstructed frame:

```bash
python sp6_extension/run.py verify --workdir /tmp/fortuity-sp6
```

If the frames are absent, this command first reconstructs them. To
repeat the original-versus-cycle comparison:

```bash
python sp6_extension/run.py benchmark --workdir /tmp/fortuity-sp6
```

The Sp(5) recovery and exact spectator-reuse comparison can be repeated
after the Sp(6) small-complex stage:

```bash
python sp6_extension/run.py recover-sp5 --workdir /tmp/fortuity-sp6
python sp6_extension/run.py reuse --workdir /tmp/fortuity-sp6
```

`--base` can point to the parent package explicitly when this folder is
stored elsewhere. The `workdir` contains generated checkpoints, not
changes to the native field or differential conventions.
