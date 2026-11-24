# Exact certificate transport

The three largest compressed certificate files are stored as bounded binary
chunks because the available GitHub transfer did not complete for a large
single blob. This changes storage only. Their original bytes, file names,
SHA-256 values, mathematical coefficients and all verification commands remain
the same.

After a fresh checkout, run this once before the commands in README.md:

```sh
python3 -B algorithms/su3_generator_20261124/assemble_artifacts.py
```

`ARTIFACTS.json` records every part, its exact byte count and SHA-256, and the
same information for the assembled original. The standard-library assembler
checks all parts and every entire file. It refuses an already present file
with incorrect bytes. The public verification runner performs this step
automatically. Reconstituted originals are ignored by Git so a verification
job leaves its checkout clean.

The original `MANIFEST.json` continues to validate the original 34-file
mathematical package and six repository dependencies. No cochain, coordinate,
relation, matrix entry or homogeneous direction is removed by this packaging.
