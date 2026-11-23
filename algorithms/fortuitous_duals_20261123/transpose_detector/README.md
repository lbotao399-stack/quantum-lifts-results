# Sparse source scoring by transposing Q0

The complete rational Sp(3) highest-weight block at native weight
`(3,3,5,0,0)`, length nine, has a forty-term target coefficient functional.
`build_transpose.py` applies the transpose of the actual pinned classical
differential to that functional. It enumerates all possible adjacent
inverse splittings in the forty supported cyclic trace words and asks the
original Q0 implementation for every resulting sign and coefficient.

The 187 candidate source words give 186 nonzero integer coefficients, with
maximum absolute value 522. Of those words, 106 occur in the original
4264-term rational omega. The score is exactly one, with no cochain
normalization change. A new candidate can be scored without constructing
its full differential.

## Complete verification

- All 2870 formal highest-weight source basis rows agree with direct Q0 scoring.
- All 2624 formal predecessor images have score zero.
- All 181 certified rational finite-colour source identities have score zero.
- The original omega has score one by both methods.

The finite-colour descent statement is for this complete highest-weight
block. An arbitrary unprojected expression must first pass a highest-weight
test or be projected. Closure is also a separate prerequisite. A coefficient
functional is not a Hilbert-space inner product or a physical operator.

## Use

The standalone scorer uses only Python's standard library:

```bash
python3 score_candidate.py
python3 score_candidate.py /path/to/candidate.json.gz
```

Input uses the same `weight`, `length`, and `[rational coefficient, ordered
native traces]` entries as `../positive_sector/omega_rational.json.gz`.
The default command scores that unchanged source and prints `1`.

For full regeneration, run the sibling positive-sector reproduction first,
then `python3 build_transpose.py`. It consumes the final forty-word catalogue
in `../enumeration/sp3_intrinsic_dual.json`, not an exploratory search file.

## Timing scope

`transpose_result.json` records seven paired trials after an explicit warmup
of both methods. Trial order alternates. Every timed result is checked to
equal the same rational value one. Each comparison includes either (a) a
full Q0 expansion plus target scoring or (b) sparse source scoring, on an
already loaded cochain. It excludes parsing, detector construction, closure
certification and the census. No end-to-end enumeration speedup is inferred.

The exact certificate, rather than these machine-dependent timings, is the
reproducible mathematical result.
