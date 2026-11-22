# Reproduce the finite checks

- Fetch O32_C_even.json.gz at the pinned repository link in ALL_ORDER_PROOF.md, using the original gzip bytes. The SHA-256 is checked in both replayers.
- Unzip the adjacent dot-r007-exact-input audit ZIP and set R007_EXTRACT to its extraction root (containing q0-certificate). The frozen dependencies are hash checked before import. The traceless projector requires SymPy.
- Run:

      python3 -B factorized_replay.py /path/to/O32_C_even.json.gz
      python3 -B characters32.py
      python3 -B replay_small_ranks32.py /path/to/pinned/taylor_hooks/native
      python3 -B block_witness32.py /path/to/O32_C_even.json.gz /path/to/j24/o0_classical_certificate.json

  The exact native source is commit 692864c36c9586533b0d459a44edfa6661beabea, directory algorithms/dual_singlets_20261118/taylor_hooks/native, also used by dot-r007-all-order. A C++17 compiler is needed only for that command. The j24 dual is in the adjacent dot-j24-formal-completion ZIP. Its complete validation is described there.
- The 47 factor records encode Eq. (2.26) in the listed order; the ungraded negative control must have 716 monomials and the graded syzygy zero.
- The prefactor nonzero test follows directly by the explicit diagonal substitution in ALL_ORDER_PROOF.md; it does not require sampling, a full rank projection at j=32, or computing any quantum-product coefficient.
