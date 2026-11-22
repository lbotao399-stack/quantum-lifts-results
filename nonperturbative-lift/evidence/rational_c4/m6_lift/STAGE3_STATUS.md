# M6 compact-chain third primitive

Scope: the preserved SU(2) rational B5 leading representative, weights (9,9), and its perturbative lift. This is not a completed expanded M6 lift or a physical instanton calculation.

## Exact chain

|Component|Terms|Status|
|---|---:|---|
|Leading|423|Exact; equivalent to the 885-term source modulo a certified Q0 boundary|
|First|31263|Exact full-Q0 equation; saved readback passed|
|Second|264322|Exact full-Q0 equation and third Cartan constraint; saved readback passed|
|Third|934381|Exact full-Q0 equation and independent saved readback passed|

The third target has 636836 terms. Its full-Q0 closure and zero Cartan restriction passed. The actual-chain d1 certificate compares two 575497-term rational polynomials with zero residual. The cohomology argument only uses Cartan injectivity on the 12-dimensional d1 kernel in this sector; it does not assert that Cartan detects all H(Q0).

## Third expanded witness

- File: `compact_first/lift/component3_exact.pkl`.
- SHA256: `7af3a5300cdb9f679090e064e3a89cb3e36b3238c808b5a0d8b6b351268b8cee`.
- Rational matrix rank 3532; 3527 nonzero solved coefficients.
- `compact_first/lift/level3_exact_summary.json`: all 8 disjoint output buckets have exact integer residual 0; full-Q0 gate 661.0563111305237 seconds.
- `compact_first/lift/level3_full_residual.pkl`: empty rational polynomial.
- `compact_first/lift/level3_independent_readback.json`: all 4 fresh-process output buckets have exact residual 0; all 934381 monomials have grading `(Nb,Nc,n1,n2)=(2,5,9,9)`; full-Q0 gate 350.42749094963074 seconds.
- The independently reopened target SHA256 is `fccc388a396e5672809a2e565b954bdb4408a13ebf468efb530b6588cde57eeb`.
- `compact_first/lift/level3_readback_run.json`: returncode 0; peak RSS 1626464 KiB, below the authorized single-worker guard 2400000 KiB.

Partition into eight output buckets changes memory use only: every output monomial remains a complete dictionary key. All bare-ghost terms and original ordered differential signs are included. Ordinary canonical Fraction coefficients are retained; equal immutable denominator integers are shared only to reduce storage.

## Still unexpanded

- Q1 and Q2 of the third component.
- The fourth target and its joint Q0/Cartan-Q1 primitive; the relevant joint rank is 2597+6=2603.
- The last target and its pure-c primitive; the relevant Q0 rank is 418.
- The order 4 and 5 full perturbative closure equations and the assembled full representative.

The source's stored rational operators obey `Q_n(physical)=(64*pi^2)^(-n)Q_n(stored)`. The rational B5 leading normalization has not been identified with the article's complete M6 secondary-operation normalization. Physical nonperturbative identification remains separate.

The current phase is complete and stopped. Controller 92187 and child 93681 exited successfully and released their leases. No subsequent images or solver were started.
