# A lower-page certificate for the eighth glueball power in SO(13)

In the relative holomorphic local-operator complex of undeformed four-dimensional $\mathcal N=4$ SYM, write $S=\operatorname{Tr}_{13}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$ using the defining trace. Let $E_1=H(Q_0)$ and $E_{r+1}=H(E_r,d_r)$. The calculation below gives $[S^8]\ne0$ on $E_3$. With Dot's level-42 source $x_0$ and Dot's choice of coupling and trace order, the quantum ring then gives, at level $j=2(E+J_L)=48$,

$$
d_3([x_0][S])=[d_3x_0][S]=-921600[S^8]\ne0.
$$

The page product multiplies leading terms. A cochain representative through order two is the quantum product of Dot's corrected source with $S$, truncated through that order. All bare $Q_r$ annihilate $S$. Its product with a pure gaugino power has no positive-order term by flavour-degree counting. Dot's uniform pure-subring threshold already predicts survival at this rank. The new calculation supplies a separate finite-colour certificate and the displayed product consequence. The two source normalizations in the linked replays are kept in their respective coupling conventions. Identification with physical local operators assumes the quantum comparison in Appendix B. [Dot's source and historical threshold](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md) · [Quantum product and comparison](https://arxiv.org/abs/2611.24160)

### The two earlier differentials

At pure gaugino weight the stable tree quotient is one-dimensional. Inside a trace the classical second jets impose graded commutation, repeated odd letters vanish, and an odd-length orthogonal trace vanishes after sorting. The balanced finite-colour relation dimensions at dotted and flavour degrees $(8,8;0,0,0)$ and $(7,7;1,1,1)$ are both zero. The three graded inputs of a bare $Q_2$ would need six flavour derivatives. A source with one second jet and otherwise first jets supplies at most four. The possible incoming second pairing therefore uses only the twice-iterated one-loop action and the classical homotopy.

Let $\lambda$ extract the pure glueball coordinate. The pinned stable duals obey $\lambda Q_1=L_1Q_0$ and $-L_1Q_1=L_2Q_0$. The replay tests the identities on all 33 and 1266 possibly contributing columns in its inverse-word enumeration. The first dual descends because the finite relation space at its weight is zero. At $(6,6;2,2,2)$ the flavour-singlet relation dimension is 301 at balanced dotted weight, and 81 after also taking the dotted singlet. Let $C_i$ be the stable distinct-plane coordinates in the zero-based order of the pinned program. The corrected second dual is

$$
L_2+\frac{9}{1280}C_1+\frac1{80}C_2 .
$$

The $C_i$ vanish on classical boundaries. To prove that the corrected row descends, consider the specialization kernel $I$ from formal orthogonal traces to thirteen colours. The native $Q_1$ commutes with specialization, including its colour-dependent empty traces. If $X$ is a predecessor relation, $Q_1X$ belongs to the preceding first-jet relation space, which has dimension zero. Thus $L_2Q_0X=-L_1Q_1X=0$. It suffices to span the relation space modulo $Q_0I$.

### Spanning the quotient by internal commutators

The orthogonal second fundamental theorem generates every relation by antisymmetrizing fourteen left endpoints against fourteen right endpoints and pairing the remaining defining indices. Removing the antisymmetrized contractions gives open paths and closed cycles. A nonsingle open path of length $k$ costs $k-1$ remaining pairings, a closed cycle costs $k$. There are at most $18-14=4$ pairings. A determinant with an adjacent graded commutator is $Q_0$ of the same determinant with the corresponding second jet inserted. This is a boundary inside $I$. The code checks all twelve nontrivial first-jet commutator rules, including each repeated odd species. Internal sorting leaves at most one of each dotted species per path. A sorted closed cycle must have even length. For each dotted species, the repeated identical odd singleton paths vanish within each of the three symmetric auxiliary endpoint types. [Orthogonal invariant tensor theorem](https://arxiv.org/abs/1102.3221) · [Path construction](https://github.com/lbotao399-stack/quantum-lifts-results/blob/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/orthogonal_glueball_relation_span.py)

Write $m,p,s$ for the minus gaugino, plus gaugino, and a scalar slot inside a path. Independently enumerating all multisets of cyclic or open words before sorting gives 76 possible raw types. Repeated odd letters or odd cycles remove 66. The ten remaining types give the following eight sorted multisets. Parentheses denote the one closed spectator.

| Open nonsingle paths and closed spectators | Remaining plus, minus singletons |
|---|---|
| $mp,mp,mp$ | $3,3$ |
| $mp,mp,mp,mp$ | $2,2$ |
| $mp,mp,mp,ms$ | $3,2$ |
| $mp,mp,mp,ps$ | $2,3$ |
| $mp,mp,mp,ss$ | $3,3$ |
| $mp,mp,mps$ | $3,3$ |
| $mp,mp,ms,ps$ | $3,3$ |
| $mp,mp,(mp)$ | $3,3$ |

An even-length sorted open word is symmetric in colour after graded reordering. Its auxiliary skew form has one component. An odd-length word has the three symmetric auxiliary endpoint components. The three singleton matrices for each scalar flavour can be chosen through a general linear coupling of flavour to this auxiliary three-space. A generic such coupling is invertible. In the flavour projection it can be absorbed by a change of the three marker columns, acting by its determinant squared and transforming the directions of scalar slots on nonsingle paths. Polynomiality extends the resulting span to singular couplings. Polarization recovers coefficients with different singleton endpoints. Odd singleton choices give the usual exterior basis for each dotted species.

For the six scalar slots the invariant flavour projection is the degree-six functional $\det(\partial)^2/144$ on the entries of a three-by-three flavour matrix. Indeed the normalized Haar average on $SU(3)$ is invariant on both sides; the bi-invariant component of degree six is the line spanned by the squared determinant, on which this functional takes value one. On the 90 balanced ordered assignments the integer projector numerator has rank five and square equal to 144 times itself. At most two scalar slots lie on nonsingle paths, and the code tries all three coordinate directions at each slot. Hence the 56 polarizations cover the flavour-singlet, dotted-weight-zero quotient relevant to the target. Central $-I$ identifies the adjoint invariants of $O(13)$ and $SO(13)$.

The ranks before and after appending $L_2$ are both three. The minor on rows $0,2,28$ and columns $0,1,2$ is $-200385994162176000$. Re-expanding ordered trace coefficients and changing basis reconstructs every entry of all 56 rows in 257 stable coordinates, including the closed-spectator row and all vanishing physical-plane coordinates. The same functional fits 280 additional raw ordered-path samples, nine of which have nonzero uncorrected scores. At $SO(11)$, $S^6$, the same implementation recovers the previously independently tested functional $(-1/32,-3/16,-1/4)$. [Earlier full ordered-word control](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36)

### Replay

Use Python 3.12 or later with NumPy, SciPy, and python-flint. In a clone of the [public trace runtime](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7), check out commit `ee9e06494a926fb88baeb6fd875cee36762f63a7`. From this certificate directory run the following commands with `REPO` set to the absolute path of that checkout.

```sh
python3 second_path_quotient.py --repo "$REPO" --colours 13 --power 8 --workers 3 --out replay13.json
cmp replay13.json so13-second.json
gzip -dc power-13-8.json.gz > raw-rows.json
python3 audit_paths.py --repo "$REPO" --case replay13.json --workers 3 --uncompressed raw-rows.json --out audit.json
cmp audit.json so13-audit.json
python3 second_path_quotient.py --repo "$REPO" --colours 11 --power 6 --workers 2 --out replay11.json
cmp replay11.json so11-second.json
```

Adding `--rebuild-uncompressed` to the audit command also regenerates all 280 sampled Pfaffian rows from their stored endpoint specifications.
