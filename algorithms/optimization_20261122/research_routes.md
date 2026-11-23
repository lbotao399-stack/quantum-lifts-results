# Research routes examined in this run

The working baseline is commit
`0b2e1a6724361fb2666089bfbe0c1a36a23d48a1`, selected from the recent entries in
`RESEARCH_DISCUSSION.md`. The default-branch README and STATUS are older than
that discussion. The baseline already contains exact dual obstruction tests,
trace-relation constructions, and higher-page augmented-rank calculations.

## Demonstrated improvement: generate the required trace orbits directly

The original multilinear routine visits every permutation for every partition
of the trace lengths, then removes duplicate cyclic and reversal orbits. The
new routine chooses the trace containing the least unused label, fixes its
cyclic starting point, fixes its reversal orientation, and recurses on the
remaining labels. It emits the same ordered coordinate list without those
redundant visits. A separate fixed-content necklace recursion handles repeated
ordinary even colour letters. The reversal test includes the original sign
and rejects precisely the self-reversing odd-length traces.

This adapts established combinatorial ideas to this repository; it is not a
claim of a new necklace-enumeration theorem. The primary reference is Karim,
Sawada, Alamgir and Husnine, [Generating bracelets with fixed
content](https://doi.org/10.1016/j.tcs.2012.11.024), Theoretical Computer Science
475 (2013), 103–112 ([author
manuscript](https://www.socs.uoguelph.ca/~sawada/papers/fix-brace.pdf)). The
constant-amortized-time theorem for their optimized algorithm is not asserted
for the implementation here.

The native adapter also builds the existing degree-fourteen F4 catalogue
directly and preserves all source-coordinate indices. The recorded benchmark
distinguishes enumeration, catalogue construction, and the complete archived
source-recipe replay. The latter has a smaller measured gain. No total
finite-colour solver acceleration is inferred from the enumeration timings.

## Exact constraint solving: executed, with no stronger source bound

The separate [native Bianchi experiment](../alternative_method_20261122/README.md)
replaces successive fixed-preimage choices by a simultaneous solve retaining
all 936 positive-order native mark coefficients through degree three. Every
coefficient of the scalar commuting equations and the mark Bianchi equation
is checked over the rationals. This enlarges the allowed leading-mark space
on seven of thirteen archived curves. All additional evaluation constraints
are already in the previous combined span. The two local bounds remain one
and two, and the global F4 rank question remains unresolved.

This was a concrete test of an exact constraint-based alternative, not just
a proposed use of a bootstrap. The result also explains why merely keeping
more homogeneous corrections does not settle this particular source.

The supplement's equation `eq-E-rank` already supplies a higher-page
augmented-rank formula, and
`nonperturbative-lift/research/23-operator-data.md` already discusses formal
Smith valuations. Repackaging those constructions without a measured
implementation gain would repeat existing work.

## Giant gravitons: useful relation organization, insufficient differential data

Ji Hoon Lee's [Trace relations and open string
vacua](https://arxiv.org/abs/2312.00242) connects finite-rank trace relations
and their successive relations to a resolution, with an explicit giant-brane
interpretation in the half-BPS setting. This motivates organizing the search
by relation generators and their syzygies. It does not by itself determine
the interacting quantum differential in the full relative jet complex.

Liu and Rajappa's [Finite N indices and the giant graviton
expansion](https://arxiv.org/abs/2212.05408) also demonstrates nonuniqueness of
giant-graviton expansions once multiple giants contribute. An index fixes an
alternating count; it does not supply individual source representatives,
corrections, or lower-page certificates. No index coefficient was promoted
to a lift claim in this run. The repository already uses relation generators,
product bounds and local syzygies extensively.

## Exceptional tensors: a concrete remaining enlargement

Gandhi, Savage and Zainoulline's [Diagrammatics for
F4](https://arxiv.org/abs/2107.12464) proves fullness of a diagrammatic tensor
functor built from the Albert representation. It supports a systematic
spanning search with the native 26-dimensional representation, its metric
and cubic invariant. Fullness does not supply a complete faithful local
rewriting system in arbitrary degree. Native exact evaluations and an
independent character dimension are still needed for a completeness claim.

The pinned repository already finds 2877 degree-eight invariants from defining
traces against a total of 2891, and already supplies the missing fourteen
with Albert contractions. Thus simply proposing Albert tensors is not a new
result. The recent discussion instead identifies transport through the
26-dimensional covariant across the ghost brackets as a remaining network
construction. This run did not establish its corrected quantum pairing.

The new trace enumerator applies to closed trace factors and trace spectators.
An open matrix word between two Albert covariants has no cyclic symmetry;
it must not be fed to the closed-trace canonicalizer. Its reversal exchanges
the endpoints and requires the appropriate signs. This is a material boundary
on reusing the optimization in a larger tensor-network calculation.

## Claims made

The positive result is a measured, exactly checked trace-enumeration and
catalogue improvement. The Bianchi result is a recorded negative experiment
with a more complete fixed-curve ansatz. New physical lifts: zero. New
fortuitous generator families: zero.
