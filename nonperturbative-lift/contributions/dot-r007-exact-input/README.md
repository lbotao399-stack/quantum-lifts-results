# Exact classical closure and first-loop input for SU(3) R-007

This checks the specific public 6139-term brane-fused representative at level 30, in the declared relative Taylor-jet complex. It does not assert that the perturbative differential lifts it, that it survives all orders, or that it has a nonzero instanton matrix element.

## What was verified

1. **Classical closure over the rationals at finite SU(3).** The proof reduces the variation to six highest-flavour generic traceless 3-by-3 colour identities. Their exact rational lowering orbits span all 444 required tensor components. Assembly matches the entire 6139-term public polynomial. A separate rank-4 computation gives a nonzero control.
2. **First-loop right-hand side.** After checking the charges and charge-conjugation parity, the derivative-free input kernel reconstructs all 213104 rational coefficients of the pinned public Q1 image in its normalization. This is an equality of polynomials, not a boundary test or a derivation of the physical coupling.
3. **A classical separating functional.** Independently expand the SU(2) block restriction with one scalar shifted by `phi1 = diag(X + h I_2, -2 h)`, all other matrices `diag(X,0)`. In unscaled Pauli coordinates the h-cubed coefficient is exactly `-(5184/5) i` times the frozen 117-term Cartesian j24 operator. Applying its rational separating dual gives `31104 i`; no comparison to a differently normalized printed coefficient is needed. The block restriction intertwines the zero-jet commutator differential (the central scalar shift drops out; the Pauli structure constant convention is an overall nonzero rescaling). The full relative classical dual annihilates all j24 zero-jet boundaries. Projection to zero jets commutes with the formal differential because every positive-order graph inserts dotted jets. Thus the same separating functional excludes any u-power boundary for this particular representative **if** an all-order closed completion is found. It does not produce the completion.

## Replay

Unzip `phase34-public-q0-q1-replay.zip` (1462493 bytes, SHA256 `1f2d7e10b7e25ad5cd2f0213f577566f636e92343f25f32334a1eada58e38fa4`). Run the three Q0 scripts listed in `q0-certificate/q0/PROOF.md` (SymPy 1.14) and `python3 -B q1/replay.py` (standard library only). The second machine replay verified all entries; the Q1 job took 33.35 seconds and 916180 KiB peak RSS.

Run the separate no-dependency block witness with the bare source from this ZIP and the dual from the [frozen SU(2) completion](../dot-j24-formal-completion):

    python3 -B block_witness.py extracted-r007/q0-certificate/source/dmk_o330_hw.json extracted-j24/o0_classical_certificate.json

The two input hashes and output identity are pinned in the script and receipt. For the geometric reason for the block choice, see [Brane-Fused Black Hole Operators, arXiv:2412.08695, Sec. 2.2](https://arxiv.org/html/2412.08695).

## Open gates

The complete rational primitive for the first-loop image is not included. The reported modular F1 completion is an independent input, not verified here over Q or matched to this exact representative. The [Cartan image theorem](../dot-cartan-image/PROOF.md) restricts the second obstruction but does not by itself exclude the possible Cartan-null direction on the second page.
