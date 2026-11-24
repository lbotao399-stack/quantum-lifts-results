# Complete physical certificates: portable reproduction

The mathematical F02 manifest is unchanged. Nine large files are transported
as 49 lossless chunks; `ARTIFACTS.json` records every chunk and original-file
SHA-256. No coefficient or matrix entry is omitted. The assembler rejects a
corrupted chunk and a mismatching pre-existing assembled file.

From a fresh repository checkout:

```bash
python3 -m pip install -r research-roadmap/verification/requirements.txt
python3 algorithms/physical_block_20261124/assemble_artifacts.py
python3 algorithms/physical_block_20261124/freeze_manifest.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algorithms/physical_block_20261124/reproduce.py --full
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algorithms/physical_block_20261124/f03/channel.py
```

The recorded execution command in the frozen manifest retains the original
machine's Python dependency path as provenance. It is not required when the
declared dependencies are installed normally. The full F02 replay regenerates
execution timings and receipts; use the manifest check before regeneration
to verify the originally published packet.

The two `physical_f02_root_*.json` files record additional independent review.
The separate `f03/MANIFEST.json` and `f03/REVIEW.md` fix the first quantum-order
physical comparison and its narrower scope. No F03 result is part of the
frozen F02 mathematical manifest.
