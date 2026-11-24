# Public files and portable reproduction

`publication_manifest.json` is the exact allowlist for this packet. It records
the SHA-256 and byte count of every permitted file below this directory, plus
all 100 pinned public input blobs. The manifest itself is excluded from its
own file digest list. `publication.py` checks the allowlist and rejects extra
files apart from Python bytecode caches. The excluded-route provenance is a
textual audit record; none of the stopped-route source files is included.

## Fetch and authenticate

Run these commands from this directory in a Python 3.12 environment. The
retained receipts used Python 3.12.14 and the exact package versions below.

```sh
python3 -m pip install -r requirements.txt
python3 -B publication.py
python3 -B survival_inputs.py --download --destination /tmp/f4-public
```

The downloader reads only the immutable paths listed in
`survival_input_manifest.json` and `proof_dependencies.json`; it checks the
byte count and SHA-256 before saving each blob. Existing authenticated blobs
are reused. It does not execute downloaded code.

An existing public Git checkout can supply the same input bytes:

```sh
python3 -B survival_inputs.py --repository /path/to/quantum-lifts-results --destination /tmp/f4-public
```

To authenticate an already materialized input tree without fetching:

```sh
python3 -B survival_inputs.py --destination /tmp/f4-public
```

All three modes finish with 100 verified public files: 98 numerical or code
inputs and two additional theorem/provenance files, totaling 2724659 bytes.

## Decisive small rational replay

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 -B verify_packet.py --inputs /tmp/f4-public --out /tmp/f4-packet-check.json
```

This command recomputes the complete 172-by-96 rational product equations and
the complete 20-by-60 bilinear classical action, including the quantum row
identity and all 2704 Casimir entries. It authenticates every saved literal
extension entry and its source-code hashes. It also checks the historical
complete-minor arithmetic. It does not regenerate the historical 3002-column
boundary block or the full order-3017 determinant.

The outgoing-rank prerequisite has its own permitted replay:

```sh
python3 -B dependency_replay.py outgoing --inputs /tmp/f4-public --out /tmp/f4-outgoing.json
```

This reconstructs the exact character dimensions and literal Cartan Gram
forms and computes both saved order-703 determinants. It invokes no native
column rebuilder. The complete-rank theorem at the published immutable pin
remains an explicitly attributed prerequisite.

## Full literal extension replay

The following regenerates all 3017-by-50 literal tensor entries. It uses at
most two workers and writes resumable checkpoints outside this packet.

```sh
export F4_PUBLIC_SOURCE=/tmp/f4-public/c778cd7c3ceb72ad8082d791b2aee5f0df082164/algorithms/dual_singlets_20261118/f4/j36_source/focused_powers
g++ -O3 -std=c++17 -shared -fPIC albert_tensor_runner.cpp -o /tmp/f4-albert.so
python3 -B albert_tensor_runner.py --inputs "$F4_PUBLIC_SOURCE" --library /tmp/f4-albert.so --points 12 32 --out /tmp/f4-controls.json
python3 -B rebuild_albert_extension.py --inputs "$F4_PUBLIC_SOURCE" --library /tmp/f4-albert.so --workers 2 --out /tmp/f4-literal-extension
```

The decompressed little-endian unsigned-16 bank must have SHA-256
`9a762aea923e782cd4e7f14680d088b5d673a3ee1a1fad2081e2a4e030f5ff5a`.
The deterministic gzip must have SHA-256
`611f1be43d43b1d16310efb07e3c60d573baec6000f9fbf006db1dacc4d5516c`.
These 150850 entries comprise 49 literal boundary columns and the closed
product. Their verification scope is not the old 3002-column block.

## Additional inventories

The remaining scripts and receipts document the complete character counts,
earlier correction freedoms and seventh-power prerequisites listed in
`README.md` and `D04_STATUS.md`. They establish no complete seventh-power
second-map decision. All input ranks that were only read from a published
certificate are identified as such in their receipts.

The publication manifest is generated only after the files are final:

```sh
python3 -B publication.py --write
python3 -B publication.py
```

No physical instanton matrix element or untwisted protection conclusion is
inferred from this formal perturbative packet.
