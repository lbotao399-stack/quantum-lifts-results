# Native polynomial Bianchi experiment

**Result:** retaining all positive-order mark corrections enlarges the allowed
leading mark spaces on seven of the thirteen archived commuting curves tested.
It supplies no additional independent local evaluation constraint. The two
degree-seven source bounds remain **1 and 2**. This is not a new quantum lift,
not a global vanishing result, and not a demonstrated search speedup.

The baseline is commit
`0b2e1a6724361fb2666089bfbe0c1a36a23d48a1` of
`lbotao399-stack/quantum-lifts-results`, directory
`algorithms/dual_singlets_20261118/f4/j36_source/first_page`.

## Alternative construction

The archived `nilpotent_curved_arcs.mark_columns` constructs a polynomial mark
by choosing a particular preimage at each successive order. That is a valid
ansatz for producing exact families, but it sets homogeneous higher mark
corrections to zero.

`native_bianchi.py` instead keeps every native coefficient of each positive
power of the curve parameter. For four Taylor axes and the 52-dimensional
adjoint, each power has `6 * 52 = 312` two-form coefficients. At degree three
there are **936** independent higher-mark unknowns. The leading mark belongs
to the unchanged archived rational frame: 70 coordinates for `(8,19,25,46)`
and 66 for `(8,19,27,47)`.

For each exact scalar curve, the code constructs **all coefficients of the
original equation**

$$[X(t),Y(t)]=0.$$

Here the bracket includes the alternating Taylor-axis product, as in the
baseline's local Bianchi equation. The curve is held fixed. Exact rational
column elimination orders the higher-mark unknowns before the leading marks.
Each free leading coordinate then supplies an explicit polynomial mark, with
all required positive-order coefficients. Purely higher-order solutions have
zero leading evaluation and need not be emitted.

Every emitted mark is substituted back into the native bracket and every
polynomial coefficient is checked to vanish. Every scalar-commutator
coefficient is independently checked as well. The construction retains all
native colour and Taylor-axis equations; it introduces no reduced equation
space. The old allowed leading space is checked to lie in the new one by an
exact rational rank test. The new evaluations are also checked to annihilate
the archived incoming matrix before adjoining them to the existing constraints.

## Executed results

| Root tuple; short axis is 19 | Distinct archived curves | Old allowed leading dimensions | Simultaneous native dimensions | New independent evaluations | Final invariant bound |
|---|---:|---|---|---:|---:|
| `(8,19,25,46)` | 3 | `48,44,44` | `48,48,48` | 0 | 1 |
| `(8,19,27,47)` | 10 | `38,37,45,45,38,37,45,45,44,45` | `45,45,45,45,45,45,45,45,45,45` | 0 | 2 |

For the first curve of the second tuple, the native matrix has 933 nonempty
equation rows, 936 higher unknowns and 66 leading coordinates; its rational
rank is 587. Its allowed leading dimension increases from 38 to 45. The second
curve has 953 rows and rank 579; its allowed leading dimension increases from
37 to 45. These are gains in the ansatz's completeness, but the additional
functionals are already spanned by the previously combined curve tests.

A separate extension of the original polynomial-family search to maximum
degree five and twelve coefficient samples also left both bounds unchanged.
The exact original equations are verified by the baseline for every accepted
family. Increasing polynomial degree alone did not improve this search.

The native simultaneous solve was slower than the fixed-preimage ansatz on
each tested curve. The JSON files record the measured times. No speed
improvement is claimed.

## Reproduce

Requires Python 3 and `python-flint` (tested with version 0.9.0). The baseline
uses NumPy and SymPy elsewhere in its imported colour construction.

Set `BASELINE` to the exact pinned `first_page` directory, then run from this
directory:

```bash
python test_native_bianchi.py "$BASELINE" --site 25 --degree 3 --limit 12 --out /tmp/native_marks25.json
python test_native_bianchi.py "$BASELINE" --site 27 --degree 3 --limit 12 --out /tmp/native_marks27.json
python extend_commuting_arcs.py "$BASELINE" --site 25 --degree 5 --samples 12 --out /tmp/site25_degree5.json
python extend_commuting_arcs.py "$BASELINE" --site 27 --degree 5 --samples 12 --out /tmp/site27_degree5.json
```

Compare mathematical counts, ranks and bounds with the saved JSON; elapsed
times will differ. `new_native_evaluations` is empty in both native test files.
The scripts do not change the baseline. The archived finite-stabilizer test is
reused only to compare the same local source bound after adding the native
equations. No global source or quantum differential is computed here.

## Independent native catalogue audit

`audit_native_catalogue.py` and `native_catalogue_audit.json` separately audit
the compiled trace catalogue in the adjacent `optimization_20261122` directory.
They are not evidence for the local Bianchi experiment. The audit checks
twenty dtype/shape cases against independently written tuple rotations and
reversals, eleven edge contents including the maximal supported length 31,
and invalid Python input handling. A standalone undefined-behavior-sanitized
C++ executable checks invalid native counts, valid boundary lengths and 10,000
random canonicalizations. The receipt pins both native and wrapper hashes.
Full degree-fourteen coordinate equality and performance are verified in the
optimizer's own report.

```bash
python audit_native_catalogue.py ../optimization_20261122 --out /tmp/native_catalogue_audit.json
```
