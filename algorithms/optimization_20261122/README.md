# Exact trace enumeration without permutation deduplication

This contribution accelerates two existing trace-enumeration calls. It keeps
every canonical trace coordinate, the ordering of the coordinates, the
transpose signs, and all published source coefficients unchanged. It does not
change a differential, add a finite-colour relation, or establish a new lift.

The baseline is the research-source commit
[`0b2e1a6724361fb2666089bfbe0c1a36a23d48a1`](https://github.com/lbotao399-stack/quantum-lifts-results/tree/0b2e1a6724361fb2666089bfbe0c1a36a23d48a1),
which is linked from the November 22 research discussion. Its complete source
tree must be checked out separately; the discussion-only default branch does
not contain all the linked production scripts.

## Measured result

These are observed wall times from fresh Python processes on the same machine.
The full observations, source hashes, and exact output hashes are in
[`benchmark.json`](benchmark.json). Peak memory includes the same baseline
module imports in both processes and is measured before output hashing.

| Labelled degree | Exact coordinates | Existing enumeration, seconds | New enumeration, seconds |
|---:|---:|---:|---:|
| 6 | 130 | 0.005463856999995187 | 0.0004159520000257544 |
| 7 | 822 | 0.05369108499996855 | 0.0024481539999214874 |
| 8 | 6,202 | 0.5403879990000178 | 0.012962569000023905 |
| 9 | 52,552 | 5.842376918000014 | 0.12026659200000722 |

At degree nine, the measured speed ratio is `48.57855220508504`; peak resident
memory changes from **272,444 KiB to 64,136 KiB**. The old implementation visits
2,903,040 permutations before deduplication; the new implementation emits
52,552 final words directly.

The new degree-ten run emits **499,194** coordinates in
`1.343708242000048` seconds, with peak resident memory **159,328 KiB**. The old
degree-ten function was not run. Its source would visit 43,545,600 permutations;
that is a combinatorial work count, not a measured timing or memory comparison.

The active F4 degree-fourteen single-trace source has content `(4,4,4,2)`.
Its old catalogue visits **3,153,150** words; the replacement constructs the
same **112,950** canonical `int32` codes. Three measured catalogue speed ratios
are `14.546945913188251`, `20.03028379034609`, and `13.437204471422687`.

**The complete archived source-row replay does not demonstrate a reliable speedup.** Reconstructing
all **474** recipes through the original `single(seed)` caller reproduces all
**8,283,801** sparse coefficients. With catalogue caches already populated, the
three paired speed ratios are `1.1493730027683784`, `1.279476352088898`, and
`0.9712092149390688`; the final run is slightly slower. This excludes subsequent elimination, finite-colour rank
jobs, and construction of the quantum differentials. There is no claim of a
fortyfold speedup of the complete lift search.

## What changes

### Labelled multi-traces

The existing function
[`mixed_colour_relations.trace_words`](https://github.com/lbotao399-stack/quantum-lifts-results/blob/0b2e1a6724361fb2666089bfbe0c1a36a23d48a1/algorithms/dual_singlets_20261118/f4/j36_source/mixed_colour_relations.py)
loops over trace-length partitions, then all permutations of the distinct
labels, and deduplicates the canonically rotated and reversed traces.

[`dihedral_traces.trace_words`](dihedral_traces.py) instead emits the final
canonical words directly:

1. The least unused label begins the first trace.
2. Choose the other labels in that trace.
3. For a trace of length at least three, choose its second and last labels in
   increasing order. Permute only the labels between them.
4. Recurse on the unused labels and sort the final output list.

The least label fixes cyclic rotation. Comparing the second and last labels
fixes reversal because all labels are distinct. The least remaining label
fixes the order of the next trace. Consequently each old basis coordinate is
generated exactly once, and the final sorting preserves all old column IDs.

The basis consists of ordinary even Lie-algebra colour letters. Reversing a
length-ℓ trace still has sign `(-1)**ℓ`; the unchanged caller applies that sign
when specializing or permuting a word. The replacement does not reinterpret
odd traces as unsigned invariants.

The independent exact count chooses the trace containing the smallest label.
For its length ℓ, there are `binomial(n-1,ℓ-1)` choices of other labels, one
orientation when ℓ=2, and `(ℓ-1)!/2` orientations when ℓ≥3. Multiply by the
count for the unused labels and sum over ℓ. With empty count 1 and degree-one
count 0, this gives

the exact counts `1, 0, 1, 1, 6, 22, 130, 822, 6202, 52552, 499194` for
degrees zero through ten, in order.

### The actual F4 source catalogue

[`native_catalogue.py`](native_catalogue.py) and
[`dihedral_catalogue.cpp`](dihedral_catalogue.cpp) implement a fixed-content
necklace recursion followed by reversal filtering. They emit the same
base-four codes used by `single_trace_flag.catalogue`. The independent Burnside
count at content `(4,4,4,2)` is

```
(3153780 + 14 * 630) / 28 = 112950.
```

The optional native canonicalizer performs exactly the same fourteen cyclic
rotations in both orientations as the existing `canonical` function. It keeps
the array shape and integer dtype, including scalar and empty-array inputs.
The native catalogue supports one through four even colours and total length
2 through 31. Odd-length reversal-invariant words are omitted because their
trace equals its negative. Arbitrary superletter parities are outside this
API; existing supertrace enumeration is unchanged.

Necklace and bracelet enumeration are established algorithms. See S. Karim,
J. Sawada, Z. Alamgir, and S. M. Husnine,
[Generating bracelets with fixed content](https://www.socs.uoguelph.ca/~sawada/papers/fix-brace.pdf),
*Theoretical Computer Science* **475** (2013), 103–112,
[doi:10.1016/j.tcs.2012.11.024](https://doi.org/10.1016/j.tcs.2012.11.024).
This is a checked integration into these trace coordinates, not a new claim in
enumeration theory. The implementation does not assert the constant-amortized
complexity of a more sophisticated bracelet generator.

### Repeated-content multi-traces

The optional pure-Python `trace_words_fixed_content` produces the same lists as
the existing labelled enumeration followed by signed identification. On the
existing degree-eight two-colour cases `(5,3)` and `(4,4)`, direct enumeration
is faster even when the old labelled list is already cached. For `(2,2,2,2)`,
it is slower than specializing an already-cached labelled list. This optional
function is therefore **not installed as a general replacement**. The measured
slower case is retained in the report.

## Exact verification

The benchmark verifies:

- Ordered-list equality at degrees 6, 7, and 8, plus matching complete output
  hashes between separate original and optimized processes at degree 9.
- Every coordinate and sign of the adjacent-transposition actions at degrees
  6, 7, and 8: 48,996 entries in total.
- Every trace evaluation at degrees 6, 7, and 8 in the native **26-dimensional
  F4 representation**, at each of the two published primes 1,000,003 and
  1,000,033: 14,308 entries in total. The representation constructor checks its
  invariant metric and Lie brackets before these evaluations.
- The entire degree-fourteen catalogue, including order and `int32` dtype.
- Every sparse coefficient in all 474 archived F4 source recipes, rather than
  a selected set of probes.

The independent root audit in `independent_verify.py` additionally checks the
canonical representatives with a separately written canonicalizer, exact
count recurrences, signed Coxeter relations, repeated-content cases, and
Burnside counts. The separate native API audit is in
`../alternative_method_20261122/audit_native_catalogue.py` and its JSON receipt;
it checks edge cases and invokes the C++ backend with undefined-behaviour
sanitization.

These checks establish equality of the implemented combinatorial coordinates
and calculations. They do not establish a new finite-colour rank, a new
cohomology class, or a nonzero higher differential.

## Reproduce

Use Python with the baseline dependencies (`numpy`, `scipy`, `python-flint`)
and a C++17 compiler. From a repository containing this contribution:

```bash
git fetch origin 0b2e1a6724361fb2666089bfbe0c1a36a23d48a1
git worktree add --detach ../quantum-lifts-pinned 0b2e1a6724361fb2666089bfbe0c1a36a23d48a1
python algorithms/optimization_20261122/benchmark.py --baseline ../quantum-lifts-pinned --out benchmark-replay.json
```

Each labelled benchmark runs in a fresh subprocess with a 45-second timeout.
Native compilation uses `c++ -O3 -std=c++17 -shared -fPIC`; its cached binary is
named by the source hash and atomically installed. Compilation or loading time
is reported separately and is not included in the native function timings.

`benchmark.optimized_calls(trace_module, source_module)` is the explicit opt-in
adapter. It replaces only the three relevant imported function references,
then restores them on exit. The benchmark runs the actual archived source
recipes through that adapter. No existing repository source is overwritten.
This is a scoped adapter, not a universal import hook. A caller that previously
used `from mixed_colour_relations import trace_words`, such as
`degree_eight_colours` or `quartic_arm_complex`, retains its copied reference;
opt that caller in with an explicit replacement of its own `trace_words` alias.
The benchmark refuses to run with assertions disabled.

## Status

Algorithm improvement: reproduced on the stated functions and complete input
sets. New physical lifts: **0**. New fortuitous generator families: **0**.
