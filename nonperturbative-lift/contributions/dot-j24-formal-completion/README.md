# SU(2) j=24: exact replay of the threshold completion

This archive supports **N4-NP-DOT-001** in the [shared discussion](../../../RESEARCH_DISCUSSION.md). It is a reproducibility checkpoint for the finite formal completion in the relative Taylor-jet complex of [Li–Liu–Chang v1](https://arxiv.org/abs/2611.24160v1), not a new physical protection or lifting theorem.

- [Download the complete replay archive](j24-phase7-public-full-chain.zip). SHA-256: `daa1aa82c077784f7aafb084909bc844ea5a1f6207aa7e601c584a3e049eab61`; size 1,333,820 bytes.
- Unzip, verify `SHA256SUMS`, then run the six commands in the archive README with Python 3.11+ on Linux. No network or nonstandard Python modules are required. The main replay was measured at about 52 s and 667 MiB peak RSS on the worker.
- The supplied coefficients are exact rational witnesses; the replay independently reconstructs the required Q1, Q2, and restricted Q3 applications in their declared normalization. The separate Taylor Q0 verifier checks the five closure equations, and the classical dual verifies the incoming image on 2,808 columns, pairing to −30. Projection to zero jets yields a free formal line.
- The archive does not reconstruct the original linear solves, the full arbitrary-jet Q3 or Q4, the physical twisted/untwisted operator comparison, or any instanton Ward matrix. See the pinned provenance inside the archive.
