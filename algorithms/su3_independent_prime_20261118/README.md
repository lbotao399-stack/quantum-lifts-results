# SU(3) independent-prime calculation and lossless map packing

This is the prime-1000033 continuation of the [SU(3) second-correction calculation](../su3_second_correction_20261118/README.md) in undeformed N=4 SYM. The first calculation, modulo 1000003, completed both map batches and all 37 second-equation elimination stages. They report consistency with rank 50147 for 77639 variables. Particular-solution residual replay and rational certification remain pending. No new lift is certified.

## Reused mathematical inputs

The same rational trace recipes are reused: 28311 first-correction directions, 49328 second-correction directions and 63983 second-target directions. The original basis evidence was obtained at prime 1000003; reusing the recipes does not assert that every evaluation frame remains injective at the second prime. Both correction equations retain the complete first-correction homogeneous freedom.

`prime_changes.json` records the original and modified source hashes. The production prime is changed only in the first-map, second-map and first-solver modules. Dependent modules import that value before defining their defaults. The first RHS reader additionally accepts the existing hash-verified rational `Q1.json.gz` archive. The original source recipes and rational Q1/Q2 coefficients are unchanged.

| Production map | Source rows | Chunks |
|---|---:|---:|
| First-equation Q0 | 28311 | 451 |
| Second-equation Q0 | 49328 | 804 |
| Second-equation Q1 | 28311 | 892 |
| Total | 105950 | 2147 |

There are also 13 first-equation and 37 second-equation RHS tasks. `tasks.json` partitions all 2197 tasks into 24 disjoint shards. Every map row occurs exactly once. The uncompressed matrix payload is 6957312128 bytes. Each production shard uses four worker processes with eight native evaluation threads per process.

## Packing small checkpoint files

The original A6 solver was observed waiting in `cl_sync_io_wait` during matrix reads. `compact_su3_maps_v1.py` addresses repeated small-file reads:

1. Verify each original manifest, source-recipe hash, recipe interval, evaluation frame and matrix hash.
2. Require complete, nonoverlapping row coverage for every source-target coupling.
3. Copy the unchanged residues into one contiguous matrix per source-target coupling.
4. Preserve the source ordering, global column offsets and augmented RHS convention used by the existing assemblers.
5. Retain a hashed inventory of all input payloads.

Packing adds no rank reduction and discards no homogeneous parameter. It requires additional disk space for the packed maps. A production speedup has not yet been measured.

The controls compare every entry of both assembled equations before and after packing at primes 1000003 and 1000033. Each check compares 126 coordinates, packs 24 blocks into 12, and rejects missing rows and corrupted payloads. The coupled-solver control also includes the false obstruction produced by freezing a particular first correction.

## Reproduce controls and run

The runtime was checked with Python 3.13 locally and Python 3.12.11 on Linux. The cluster used NumPy 2.4.6, SciPy 1.15.3 and python-flint 0.8.0. A C++17 compiler is required for the native map libraries.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python preflight.py
.venv/bin/python validate_compact_su3_maps_v1.py \
  --runtime runtime/work/sp3-native.nosync \
  --output controls.json
```

The cluster preflight completed in 86.02145722601563 seconds. It includes Q0 controls, Q1 checks across all 17 profiles, the coupled-system controls, and an independent dense check of 8940 coordinates for the JSON RHS reader and its sign. The saved Q1 control checks 4732 coordinates at both primes. Validation receipts are in `validation/`.

For production, launch `distributed_maps.py` as a scheduler array with task IDs 0 through 23, 32 CPUs per task, and `SLURM_CPUS_PER_TASK=32`. Limit BLAS libraries to one thread per worker to avoid oversubscription. Give each task a one-hour cap; the launcher stops starting new work near that limit and preserves completed checkpoints. Only after every shard has completed successfully, launch `solve_compact.py` with 64 CPUs. It verifies the packing control again, packs all three maps, solves the first equation, then couples the second equation with all first-order freedom retained. Incomplete shards prevent the dependent solve.

The running cluster schedule and progress are recorded in the [numerical snapshot](../../data/N4_SU3_independent_prime_20261118.json). At its timestamp all 24 shards held 32 CPUs each. A second-prime solution is not by itself a rational certificate. Rational reconstruction must also check compatible pivot/free-variable choices and the resulting exact equations.

## Counts and provenance

Phase increments: d1 pairings 0; audited fortuitous-generator families 0. Cumulative counts: d1 5; generator families 3; d2, d3 and nonperturbative lifts 0, 0, 0. These are supercharge-obstruction orders. The family count does not assert global ring minimality or absolute literature novelty.

The graded trace and matrix-evaluation modules originate from [davidsd/fortuity-ai](https://github.com/davidsd/fortuity-ai), with upstream pin `4256f207986029540060fb44dd494746828f65ce`. Imported headers are retained. The additions here are the independent-prime runtime, complete second-basis recipes, task partition, and lossless matrix packing. Credentials, account paths, installed wheels and native binaries are excluded.
