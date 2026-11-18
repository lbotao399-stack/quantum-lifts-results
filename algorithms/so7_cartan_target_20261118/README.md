# An incoming rational first-order pairing at level 30 in $SO(7)$ SYM

At prime $1{,}000{,}003$, two square minors have the following determinants. The entries are reductions of rational projected trace polynomials. The projection denominators are invertible at this prime. The evaluation points lift to integer antisymmetric matrices by lifting the independent upper-triangular entries and enforcing antisymmetry.

| Map | Order of the minor | Determinant modulo $1{,}000{,}003$ |
|---|---:|---:|
| Classical supercharge $Q_0$ from length 11 to length 12, at weight $(3,3,2,2,2)$ | 2,307 | 488,963 |
| Length-10 classical images with the $SF$ first-order pairing appended | 12,526 | 942,245 |

The weight lists the two holomorphic derivative degrees followed by the three flavour degrees. The length-11 space has rational dimension 14,832. The first minor and $Q_0^2=0$ give

$$
\operatorname{rank}_{\mathbb Q}(Q_0:\text{length }10\longrightarrow\text{length }11)
\leq 14{,}832-2{,}307=12{,}525 .
$$

Here $S=\operatorname{Tr}(\lambda_+\lambda_-)$ is the gaugino bilinear and $F$ is the known source whose first obstruction is $S^4$. The target $SF$ spans the classical cohomology at weight $(4,4,1,1,1)$ and length 10. Its dimension is $1{,}739-1{,}443-295=1$. In the 1,870-coordinate target evaluation, the 1,443 saved boundary rows have modular rank 1,443. The saved target functional $\ell$ satisfies $\ell Q_0=0$ on these boundaries and $\ell(SF)=1$. The rational boundary-rank upper bound is $4{,}752-3{,}309=1{,}443$. [Target frame and weights](check_incoming.py) · [Target polynomials](replay.py)

The modular functional lifts to a rational functional with the same reduction. Adjoin the evaluated $SF$ to the 1,443 evaluated boundary rows and choose a 1,444-column minor invertible modulo the prime. Lift the free coordinates of the saved functional to integers. Solving the remaining square system over the rationals enforces the 1,443 zeros and the one, and preserves the reduction of every coordinate. The rational boundary-rank bound then makes this functional vanish on every classical boundary.

The final column of the second minor is $\ell Q_1$, with $Q_1$ the first quantum supercharge. Its other 12,525 columns evaluate the classical differential. The nonzero modular determinant gives a nonzero rational determinant for the lifted rows. Together the two minors imply

$$
\operatorname{rank}_{\mathbb Q}Q_0=12{,}525,\qquad
\operatorname{rank}_{\mathbb Q}(Q_0,\ell Q_1)=12{,}526 .
$$

The signed cofactors along the pairing column give coefficients of an explicit rational linear combination of the selected 12,526 source words. Laplace expansion makes its evaluations against the 12,525 classical columns zero. Those columns already detect the full rational classical image, so the combination is $Q_0$-closed. Its pairing is the nonzero rational determinant. The supercharge identity $Q_0Q_1+Q_1Q_0=0$ makes its quantum image closed and makes the first obstruction of every boundary vanish. Hence the source defines a nonzero classical class, and its $d_1$ image spans the line of $SF$. The target is absent on every page from the second onward.

The initial 54,042 source rows span 32,287 of the 42,991 source cochains at this prime. The proof needs the indicated minor and the bound on the entire classical image. It does not require these sampled source words to span the domain. Classification of the source modulo previously known fortuitous products and descendants remains open.

## Replay

The sources for the orthogonal matrix evaluation, extremal projection and $Q_1$ are in `native/`. Use the vendored full-letter runtime and the NumPy, SciPy and python-flint requirements in the independent-prime algorithm. From this directory run

```sh
python3 build_incoming_native.py
python3 verify_pairing.py
python3 rank_successor.py --runtime ../su3_independent_prime_20261118/runtime/work/sp3-native.nosync --native native --replay successor_minor.json.gz --output /tmp/so7-successor-replay
python3 -c 'import numpy as np; z=np.load("target_frame.npz"); np.save("/tmp/so7-functional.npy",z["functional"])'
python3 incoming_minor.py --replay incoming_minor.json.gz --runtime ../su3_independent_prime_20261118/runtime/work/sp3-native.nosync --native native --functional /tmp/so7-functional.npy --output /tmp/so7-incoming-replay --workers 16
```

The last replay forms the selected rows directly from their trace words, recomputes all 212 pivot-block determinants and multiplies them in their stored order. [Successor recipes](successor_minor.json.gz) · [Augmented recipes](incoming_minor.json.gz) · [Extremal projections and evaluated rows](incoming_rows.py)

A8 array `5628225` completed all 12 shards and independently rebuilt the entire augmented minor from its trace recipes. All 212 panel determinants agreed; the replay contains the SHA-256 of the 12,526-by-12,526 residue matrix. [Full replay receipt](incoming_replay.json) · [Sharded replay driver](replay_strips.py) · [Preceding and adjacent classical trace recipes](classical_recipes).
