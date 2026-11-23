# Exact cycle expansion of the native Pfaffian recipes

The endpoint implementation enumerates every signed perfect matching
of the `2r+2` marked endpoints. At Sp(6) this is `13!! = 135135`
matchings per recipe.

The new implementation groups those same terms by their disjoint
trace cycles. It does not replace the Pfaffian relation by a different
normalization.

## Combinatorial reduction

- Every connected component of the union of a marked matching and the
  fixed paths is a trace cycle.
- Anchor each cycle at its first labelled path, then order and orient
  its other paths. Reversing a path gives the exact native reversal
  sign from `path_relation_completion.canonical_paths`.
- Each undirected cycle is obtained twice from directed orientations.
  The implementation divides by two only after summing all anchored
  orders. Dividing each individual directed contribution is generally
  invalid; its coefficients need not be even.
- Moving an odd path past preceding odd paths contributes the native
  Grassmann sign. Repeated even paths can be grouped by multiplicity;
  repeated identical odd paths vanish.
- The endpoint-first convention for a path with internal vertices has
  one descending fixed edge. Its sign is included once for each such
  path, independently of the number of its internal vertices.
- Arbitrary fixed matching graphs are reduced to paths and closed
  cycles. A matching that closes each open path fixes the overall
  orientation and Grassmann sign exactly.

`fast_paths.py` implements the recurrence with integer arithmetic.
Every final division by two is asserted to be exact. Optional caches
memoize unchanged native canonicalization functions.

## Validation

- 80 generated mixed boson/fermion path recipes at ranks 0–3 agree
  coefficient by coefficient with the endpoint expansion. Rank zero
  is used only as a graph-algebra sign test, not as a gauge theory.
- 47 actual Sp(4) and 39 actual Sp(5) frontier graph/path recipes agree
  coefficient by coefficient.
- All 438 Sp(5) identity recipes, including `Q0` components and the two
  target-completion identities, reproduce their stored integer
  polynomials exactly.
- Three actual Sp(6) paired recipes were compared against all 135135
  endpoint matchings with cold native canonicalization caches. Each
  polynomial agrees exactly; the number of resulting trace terms is
  respectively 135, 740 and 9756.

The recorded wall times are workload-specific. They establish an
improvement for the tested inputs, not a uniform speedup factor or a
new complexity theorem. The individual inputs, hashes and timings are
in `benchmark_sp6.json`; the lower-rank checks are in
`cycle_engine_validation.json` and `frontier_validation_sp4.json` /
`frontier_validation_sp5.json`.

## Universal-preimage solver

The second change is independent of the Pfaffian expansion. Given a
candidate defect, `inverse_support.py` finds possible `Q0` parents by
merging adjacent native letters. The differential preserves the
multiset of the individual trace weights, so these are separate linear
systems. In blocks with one derivative-carrying trace, the solver also
uses the native highest-weight equations.

Modular elimination proposes a preimage. Rational reconstruction and
a fresh native `Q0` computation certify it over the rationals.
Failure of a restricted source set is never treated as proof of
non-existence. The solver enlarges the source set from the residual
support, subject to its documented resource cap.

On the previously certified Sp(5) defect, all 958 trace-weight blocks
were solved and the combined 43283-term preimage was checked exactly.
The first run solved 957 blocks; the last block exceeded the initial
5000-source cap. Raising that cap to 12000 admitted its 8535
highest-weight source expressions and completed the exact check.
These data concern recovery of an existing class, not an additional
independent physical family.

For Sp(6), the candidate defect has 98882 trace terms and decomposes
into 1961 trace-weight blocks. Every block was solved over the
rationals. The largest block used 17718 highest-weight source
expressions obtained from 66308 adjacent-letter parents; its selected
image rank was 7627. This is a preimage calculation for that defect,
not a complete census of the single-trace complex.

After removing pure-scalar trace factors, 1374 of the 1961 Sp(6)
inverse tasks agree exactly, up to their recorded integer factors,
with already solved Sp(5) tasks. This leaves 587 new tasks. The
comparison concerns linear inverse problems; it does not identify or
count independent physical families. `spectator_reuse.json` records
every exact match.

An attempted simplification to one- or two-coefficient duals failed
for the Sp(5) candidate. All 66310 active trace-coefficient columns
were checked, giving 6427 proportional-column signatures. None
supports a one- or two-probe detector with zero boundary pairings and
nonzero defect pairing. The successful connected-trace dual instead
has 111 nonzero coefficient probes. This negative result prevents
extrapolating the two-probe scalar SO witnesses to these symplectic
sectors.
