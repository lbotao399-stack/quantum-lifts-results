# Replay the characteristic-zero existence proof

Read ALL_ORDER_PROOF.md for the argument and conventions. No big R-007 correction matrix is reconstructed.

1. Unzip the frozen exact R-007 input packet in the adjacent dot-r007-exact-input contribution. Point DMK_VENDOR at its extracted q0-certificate/vendor directory. The scripts hash-check dmk.py and words.py before importing them. SymPy 1.14 is needed by the frozen traceless projector.
2. Run:

       python3 -B free_polynomial.py
       python3 -B exact_small_characters.py
       python3 -B cartans_replay.py

   The last regenerates cartan-300-minor.json from exact data and compares the assembled prefactors with the original dMK traces. It checks the original selected determinant at both stated primes.

3. Run the small rank replay with the native source from the pinned threshold certificate commit 692864c36c9586533b0d459a44edfa6661beabea:

       python3 -B replay_small_ranks.py /path/to/algorithms/dual_singlets_20261118/taylor_hooks/native

   This builds a temporary C++17 engine, checks the exact source file hashes, runs the scalar evaluator, and replays the two-prime rank lower bounds. The relevant joint rank is outgoing rank plus the extra Cartan rank on the kernel of the *same* evaluated outgoing map, i.e. a lower rank of the true stacked rational map. Exact integer character dimensions from step 2 provide the upper bounds. The trace polynomials are invariant/HW by their gauge-trace, flavour determinant and dotted bracket construction. Probe selection is a method of finding nonzero residue minors, not a hypothesis that rank does not fall at special points.

4. The non-torsion proof additionally uses the block-separation replay in the adjacent dot-r007-exact-input contribution.

The proof is for the normalized formal holomorphic Taylor-jet complex. These scripts do not print the full repaired cochains or an instanton matrix element.
