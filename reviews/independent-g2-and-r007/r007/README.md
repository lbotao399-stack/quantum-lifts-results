# R-007 radial target and one complete boundary projection

Dot constructed the all-order formal R-007 source and the 5,080-term outgoing candidate $T$. The defining trace, divided Taylor jets and the free singleton oscillator metric give

$$
\|T\|^2=5171748864,\qquad \|T|_{\mathfrak h}\|^2=7987200 .
$$

Here $\mathfrak h$ is the Cartan subalgebra of $\mathfrak{su}(3)$. The first norm contracts the three quadratic traces and one symmetrized cubic trace of every term. The cubic tensor has squared norm $20/3$ on all eight colours and $2/3$ on the Cartan. The free one-letter metric includes the factorials of both dotted oscillators and of the undotted oscillator. [Dot's classical target and recipe](https://github.com/lbotao399-stack/quantum-lifts-results/blob/04fe44e3a06f3fc589f72fe7bfd8568da5afc9a3/nonperturbative-lift/contributions/dot-r007-t2-primary/PROOF.md) · [Dot's characteristic-zero all-order result for the source](https://github.com/lbotao399-stack/quantum-lifts-results/blob/04fe44e3a06f3fc589f72fe7bfd8568da5afc9a3/nonperturbative-lift/contributions/dot-r007-all-order/ALL_ORDER_PROOF.md)

For one complete predecessor profile the matrix-generator contraction gives the value $912$ on $Q_0^\dagger T$. The classical image is therefore present in this bare state. The profile consists of five copies of the first undifferentiated scalar, one copy of the second scalar, and the two divided jets encoded as $(0,2,(1,2))$ and $(2,0,(2))$ in Dot's zero-based flavour alphabet. Enumerating all cyclic products gives 175 formal words. Diagonalize the repeated scalar. The expression is degree five in its two independent diagonal entries and is multilinear in the three other letters. Evaluate the first diagonal entry at one and the second at $0,1,2,3,4,5$. For the other three entries use all 56 ordered torus-neutral triples of the eight traceless generators. These 336 coordinates determine the polynomial: homogeneity recovers the other diagonal ratios, the torus forbids the omitted triples, and regular semisimple matrices are dense. The evaluation matrix has rational rank 44. Traces generate simultaneous-conjugation invariants of matrices; there is only one odd slot in this profile.

Project $Q_0^\dagger T$ onto this profile using that distinct one-letter occupations are orthogonal. On its 44-dimensional space solve the exact normal equations for a primitive $u$. The matrix of the projected $Q_0^\dagger Q_0$ has rank 44. The saved solution uses 16 words and its boundary has 219 terms. In all 336 interpolation coordinates the residual adjoint vanishes. A separate trace pairing gives

$$
\|Q_0u\|^2=\langle Q_0u,T\rangle=\frac{401271808}{589},
\qquad
\|T-Q_0u\|^2=\frac{3045758809088}{589}.
$$

The last number is the norm after projection on this profile. Other profiles can reduce the harmonic norm further. The nonzero Cartan test survives every classical boundary subtraction. The [one-instanton propagator calculation](PROPAGATORS.md) works at separated bare insertions.

~~~sh
python3 replay_radial_checks.py
~~~

The seven replays reconstruct their receipts in a temporary directory. Use CPython 3.13 with NumPy and python-flint. The standalone polarization replay uses only the Python standard library.
