# Public Konishi normalization replay

Run `python -B konishi_control.py` (Python 3.12, SymPy 1.14 tested). Exact assertions rebuild SU(2) color and Grassmann contractions from the published phi2/phi6 profiles, then independently check wedge/determinant zero-mode saturation and the bosonic Gamma-function residue. Runtime about 0.1 seconds, RSS 47 MiB; no network or private imports.

Input source: Alday–Korchemsky, https://arxiv.org/abs/1605.06346v2, Eqs.1.3, 2.4–2.7 and 2.10. Output comparisons: Eqs.2.16–2.17, 2.23–2.27. Full paper is not redistributed. Code was independently generated from those formulas. The profile solutions themselves are source inputs, not newly solved field equations.

This is a long-multiplet normalization control, not a BV-Ward insertion, j24 certificate, BPS lift, or instanton supercharge matrix element. The reference receipt contains measured timing; replay timings can differ. Only code, this note, reference receipt and allowlist/checksums are included.
