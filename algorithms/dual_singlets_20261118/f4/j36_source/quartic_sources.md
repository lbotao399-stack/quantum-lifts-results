# Classical elimination of 80 candidate sources for the sixth glueball power

In the full relative holomorphic jet complex of $F_4$ at $j=36$, take the 20 leading recipes in Dot's archived `project-bank.json` and the four scalar dressings in the [tree calculation](dressed_tree.py). Among the 80 parameter directions, the physical classical differential $Q_0$ has rank 61. All 19 directions in its kernel are classical boundaries. Their image in classical cohomology has dimension zero. This leaves the incoming $d_2$ on the glueball power $S^6$ open for sources outside this family. [Inputs and two minors](dressed-tree.json) · [Complete replay](quartic-source-test.json)

### Physical closure from ordinary colour polynomials

The 20 tree images have two distinguished dotted slots. Their ordinary colour tensors are symmetric under exchange of those slots. Identifying the two slots with an even variable $T$ gives polynomials in $X_1,X_2,X_3,T$ of degrees $(4,4,4,2)$ after dressing. The three scalar columns $X_i$ transform by the fourth power of the determinant. The verifier splits the second jet in each of the 20 seeds, checks the two dotted assignments of the last mixed jet, and checks all six scalar raising and lowering operations separately on the seeds and dressings.

Use the archived six odd transfers $P$ from a scalar to a dotted gaugino. On the tree image, $P^6$ is top coefficient extraction after the substitutions

$$
 X_i\longmapsto X_i+\eta_{i+}\bar\lambda_{\dot+}+\eta_{i-}\bar\lambda_{\dot-},
 \qquad
 T\longmapsto \xi_+\bar\lambda_{\dot+}+\xi_-\bar\lambda_{\dot-}.
$$

Here the six $\eta$ and two $\xi$ are exterior parameters. The two $\xi$ assignments restore both original symmetric colour slots, with a common nonzero factor determined by the exterior ordering. The six transfers anticommute with $Q_0$ on all five first jets and all twelve second jets. Direct substitution gives zero in each of the 102 ordered anticommutators.

Write $K=\sum_{\substack{i=1,2,3\\\text{scalar flavours}}}(T\partial_{X_i})(X_i\partial_T)$ for the three scalar-column transfers, with each directional derivative acting on the named column. On degree two in $T$, $K$ is six times the identity plus the sum of all interchanges of a marked slot and a scalar slot. Top coefficient extraction kills $KF$ for every input polynomial $F$. After substitution the last transfer $T\partial_{X_i}$ is differentiation in $\eta_{i+},\eta_{i-}$ with coefficients $\xi_+,\xi_-$. An exterior derivative has zero top coefficient. Thus $(K-6)(K-14)F$ has exactly 84 times the physical top component of $F$.

The Pieri rule applied to the scalar rectangle $(4,4,4)$ and the symmetric pair of marked slots gives just $(6,4,4),(5,4,4,1),(4,4,4,2)$, once each in the colour tensor functor. To evaluate $K$, take the leading three-row minor of the scalar columns and the leading four-row minor after appending $T$. For $a=0,1,2$, their powers $4-a,a$, multiplied by the first entry of $T$ to power $2-a$, are highest-weight vectors of these three shapes. In $X_i\partial_T$ the derivative of the four-row minor vanishes. In the reverse derivative the cofactor identity for the three-row minor gives eigenvalue $(2-a)(7-a)$, respectively $14,6,0$.

Top coefficient extraction commutes with all linear changes of colour coordinates. On the last summand it is nonzero. For the squared three-row and four-row minors, set the three scalar columns to the first three unit vectors. In the replay's exterior order the coefficient of all eight parameters and all eight gaugino components in four colour directions is $50400$. Irreducibility gives injectivity on that summand over the rationals. This argument works on the 52 colour coordinates of $\mathfrak f_4$ before taking gauge invariants.

The physical tree map and the projected polynomial map therefore have the same parameter kernel. The determinants of order 61 are $975516\bmod 1000003$ and $444735\bmod 1000033$. The 19 independent rational null vectors have distinct free coordinates. Three projected images vanish in formal traces. The other sixteen reduce to zero using polarizations of

$$
12\operatorname{Tr}_{26}(X^4)-\big(\operatorname{Tr}_{26}(X^2)\big)^2=0
$$

on Lie-word arguments, with spectator traces. The identity follows from its 35 Cartan coefficients and density of regular semisimple elements. The reduction never uses a trace with more than eight first-jet slots. Substitution of the exterior variables preserves each polarized identity. [Quartic identity in the native representation](README.md#trace-identity-and-primitive) · [All 19 reductions](dressed_tree.py)

### A classical primitive for the finite-colour relation

Lift a closed candidate to formal orthogonal traces and call its classical derivative $R$. Let $I$ be the ideal of finite-colour identities. We show $R\in Q_0 I$ for each of the nineteen directions.

1. In a polarized quartic with a nontrivial Lie-word argument, choose an innermost nonzero graded commutator of two first jets. It is $Q_0$ of the corresponding second jet, with the divided-jet factor for repeated gauginos. Replace that bracket by its primitive inside the nested graded brackets. Every other first jet is closed. The resulting polarized quartic is still an identity because its arguments are Lie-algebra valued. Its differential is the original term, including the spectator factors.
2. For a quartic on four single letters use the map from formal first-jet traces to graded-commuting traces. It has a multiplicative section given by graded symmetrization inside each even trace. Odd-length traces vanish after graded commutation and orthogonal reversal. The quartic on single letters is fixed by projection followed by this section. In each spectator subtract its symmetrized section. Every subtracted difference is the differential of a sum of adjacent-interchange primitives, including the reversal primitive on an odd trace. Multiplication by the quartic puts these primitives in $I$.
3. The terms left over equal the section of the commuting projection of $R$. That projection vanishes because $R$ is a classical derivative of a one-excess cochain. This constructs a finite identity with derivative $R$. Subtracting it from the lifted candidate gives a closed formal trace cochain with thirteen letters and fourteen Taylor slots. Every trace still has at most eight slots.

### Odd-length stable cycles at one excess

Polarize the Taylor slots of a single trace so that they are distinct. A block of labels specifies one divided jet. The classical map splits a block into two ordered nonempty blocks with its native sign. At excess zero all blocks are singletons. At excess one there is one double block. The incoming columns have one triple block or two doubles. The verifier enumerates the signed cyclic and reversal basis and checks $Q_0^2=0$ over the integers on every incoming column.

| Total slots of the trace | Incoming columns | One-excess columns | Terminal columns | Incoming rank modulo 3 | Outgoing rank modulo 3 |
|---:|---:|---:|---:|---:|---:|
| 2 | 0 | 0 | 1 | 0 | 0 |
| 3 | 0 | 3 | 1 | 0 | 1 |
| 4 | 7 | 6 | 3 | 4 | 2 |
| 5 | 25 | 30 | 12 | 14 | 12 |
| 6 | 195 | 180 | 60 | 121 | 59 |
| 7 | 1680 | 1260 | 360 | 894 | 360 |
| 8 | 15960 | 10080 | 2520 | 7561 | 2519 |

Thus terminal cohomology at odd trace length and one-excess cohomology at odd trace length both vanish over the rationals throughout the required range. A one-letter trace vanishes identically. The multilinear calculation takes distinct odd Taylor directions. Evaluation on the actual Taylor superspace is tensoring the multilinear permutation module with the signed tensor power of that space and taking symmetric-group coinvariants, with the sign twist for the initial odd labels. This is exact over the rationals by averaging, so repeated directions and the two even Taylor coordinates preserve the vanishing.

The classical differential preserves the slots in each trace. The rational tensor-product formula for cohomology applies before averaging permutations of equal traces. At total excess one exactly one factor has excess one. Thirteen letters force at least one factor to have odd length. Its cohomology vanishes by the table. The corrected formal cycle is therefore a classical boundary. Specializing its primitive to $F_4$ gives a primitive of the original candidate, since the subtracted identity specializes to zero. This is a rank-and-existence proof for all nineteen primitives. The previously published 6417-term polynomial remains the explicit primitive for Dot's first recipe. [Signed rank calculation](stable_trace_ranks.py)

With NumPy, SciPy and python-flint, run in this directory

~~~sh
python3 quartic_source_test.py /path/to/F4-j36-d2-verification.zip --out /tmp/f4-quartic.json
cmp /tmp/f4-quartic.json quartic-source-test.json
~~~
