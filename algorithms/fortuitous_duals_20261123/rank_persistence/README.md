# Rank persistence of the Sp(3) omega block

This directory records classical highest-weight calculations at native
weight `(3,3,5,0,0)`, length 9 and `j=28`. The first two weight entries are
dotted/Taylor counts; the last three are flavour weights. `Sp(N)` means
`USp(2N)`, with the original defining trace and no rank-dependent rescaling.

## Exact dimensions and ranks

Here C8, C9 and C10 are the adjacent highest-weight trace-cochain spaces in
this one grading. The differential increases length by one.

| Group | dim C8 | dim C9 | dim C10 | Incoming rank | Outgoing rank | dim H9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sp(1) | 64 | 25 | 3 | 22 | 3 | 0 |
| Sp(2) | 1527 | 1255 | 448 | 850 | 405 | 0 |
| Sp(3) | 2559 | 2689 | 1281 | 1607 | 1081 | 1 |
| Sp(4) | 2624 | 2870 | 1451 | 1672 | 1198 | 0 |
| Sp(5) | 2624 | 2870 | 1455 | 1672 | 1198 | 0 |
| Formal stable cover | 2624 | 2870 | — | 1672 | 1198 | 0 |

The Sp(3) row uses the companion complete rational relation space, its exact
differential rank, and the rational omega and finite-map certificates.
Its class is the previously known omega; this is not a claim
of discovering a new BPS state. The lower and higher rank calculations test
the entire indicated highest-weight block.

The displayed identities are

```
25   - 22   - 3    = 0;
1255 - 850  - 405  = 0;
2689 - 1607 - 1081 = 1;
2870 - 1672 - 1198 = 0.
```

All dimensions are exact character computations. The incoming and outgoing
minors are evaluated modulo 1000003 and have nonzero determinants:

| Group | Incoming determinant | Outgoing determinant |
| --- | ---: | ---: |
| Sp(1) | 592760 | 363806 |
| Sp(2) | 827918 | 155065 |
| Sp(4) | inherited stable incoming rank | 435941 |

Nonzero minors give rational rank lower bounds. Nilpotence and the exact
middle-space dimensions force equality in the zero-cohomology rows; the
Sp(3) case requires one additional argument. The exact relation space has
dimension 181 and its differential has rank 116. Thus its closed part has
dimension `181 - 116 = 65`. Stable cohomology vanishes, so those 65 closed
relations are exactly the intersection of the stable incoming image with
the finite-rank relation space. The finite incoming rank is therefore
`1672 - 65 = 1607` over Q. The modular boundary functional pairs 1 with
omega, so appending omega increases this already exact boundary rank.
The exact rational omega is closed, giving outgoing rank at most
`2689 - 1607 - 1 = 1081`. The nonzero 1081-minor gives the opposite bound.
Consequently the Sp(3) cohomology dimension is exactly 1 over Q.

## Why the zero result holds for every N >= 4

The source-degree trace specialization to Sp(4) is onto by defining-trace
generation. Both its domain and codomain have dimension 2870, so it is
injective. Thus there are no source-degree trace identities at Sp(4).

For N >= 4 choose a symplectic direct sum of an eight-dimensional space and
its complement, and put every matrix equal to `X_Sp4` on the first summand
and zero on the second. Traces restrict without a scalar factor, including
the graded trace signs. Every Sp(N) identity would therefore restrict to an
Sp(4) identity. Consequently the source-degree relation kernels satisfy

```
I_N is contained in I_4 = 0, for N >= 4.
```

The 1672-dimensional stable incoming image therefore remains independent in
every Sp(N) source space. Every finite predecessor is trace-generated, so
the incoming rank is exactly 1672.

The nonzero Sp(4) outgoing 1198-minor uses actual matrix evaluations. The
same evaluations embed into Sp(N) by the symplectic direct sum above. Hence
every N >= 4 has outgoing rank at least 1198. The complete stable outgoing
rank is 1198, which is also an upper bound after specialization. Thus

```
dim C9(Sp(N)) = 2870;
incoming rank = 1672;
outgoing rank = 1198;
dim H9(Sp(N)) = 0, for every N >= 4.
```

**No target-degree injectivity is assumed.** Indeed C10 has dimension 1451
at Sp(4) and 1455 at Sp(5). Source injectivity and the saturated outgoing
minor are the sufficient ingredients.

Since the known Sp(3) omega is a nonzero cohomology class, it has no closed
lift to any Sp(N) with N >= 4: the restriction of a boundary is a boundary,
whereas every closed class in those larger-rank blocks is a boundary.
Its restrictions to Sp(1) and Sp(2) are boundaries. They need not be zero
cochains: the native evaluations saved here are nonzero at those ranks.

This argument concerns Q0 and this fixed highest-weight block only. It does
not establish rank persistence for other charges, products or quantum maps.

## Exact chosen-lift control

`exact_rank_control.py` loads the original `omega_rational.json.gz` from the
companion directory and differentiates its 4264 trace terms. The resulting
stable derivative has 4502 terms. At the exact integer frame with seed
159772 its six Sp(4) coordinates are

```
561727/864, 7235557/864, -70747/216,
1017713/288, -1877243/864, 428041/144.
```

The corresponding Sp(1), Sp(2) and Sp(3) control coordinates vanish.
These are exact sample values. The nonzero Sp(4) control concerns this chosen
lift; the complete minor and rank proof above establish the stronger
statement excluding every closed lift.

## Reproduce

The standalone saved-certificate replay requires only NumPy and python-flint:

```bash
python3 verify_rank_certificates.py
```

It checks hashes, recomputes all saved determinants, and verifies the exact
dimension/rank arithmetic. It does not regenerate native matrix entries.

For native generation, first run the sibling `positive_sector` pinned-input
fetch and stable-trace reconstruction. The helper `rank_runtime.py` imports
that directory's `probe_stable.py`; it contains no machine-specific source
paths. If the companion is elsewhere, set `FORTUITY_POSITIVE_SECTOR` to its
absolute directory. Then run:

```bash
python3 dimensions.py
python3 lower_rank_scan.py
python3 sp4_outgoing_rank.py
python3 exact_rank_control.py
python3 portable_smoke.py
```

`dimensions.py` has a main guard and accepts `--sizes`, `--lengths`, and
`--output`. Importing it does not start a computation or write output.
The omega loader uses the portable rational JSON; `omega_exact.pkl` is not
required. The stable-row cache is generated by the companion pipeline and
is not part of this deliverable.

Original native runs have already passed. The final portability check only
reruns saved-minor arithmetic and the inexpensive exact chosen-lift control;
it does not repeat the complete positive-sector certification.
`portable_smoke.py` copies the production scripts to a temporary directory,
checks that importing them writes no JSON and prints nothing, and reproduces
the exact rational control from its portable omega JSON.

## Inputs and files

The algorithms and original omega inputs are pinned to
`9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9`; the shared runtime is pinned to
`082be7070d207e2cdecfd0a91b03d5e014fd303d` in
`lbotao399-stack/quantum-lifts-results`.

`MANIFEST.json` selects only scripts, JSON receipts and three compact minor
archives. Logs, compiled libraries, pickle caches and copied upstream code
are excluded. The complete native sample identifiers and selected row and
column indices remain in the original receipts.
