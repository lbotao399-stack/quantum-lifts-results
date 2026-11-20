# The $F_4$ fifth-power channels

The terminal fifth gaugino power has level $j=2(E+J_L)=30$. Its incoming spin and flavour singlets have doubled charges $(1,1,10,10,10)$. At derivative bidegree $(t,t)$ their three flavour degrees are $5-t$ and their length is $14-t$. The first, second and third perturbative sources for the pure fifth power have $t=4,3,2$, respectively. The following are full invariant cochain dimensions before taking classical cohomology.

| $t$ | Length $13-t$ | Source, length $14-t$ | Target, length $15-t$ |
|---:|---:|---:|---:|
| 5 | 31 | 20 | 12 |
| 4 | 2053 | 745 | 144 |
| 3 | 18533 | 6862 | 1278 |
| 2 | 46192 | 17499 | 3092 |
| 1 | 31394 | 12353 | 2333 |
| 0 | 4295 | 1862 | 428 |

The integral character algorithm multiplies adjoint Adams operations by reflecting shifted highest weights into the dominant chamber. Symmetric and exterior invariant dimensions through degree 12 agree with the independent generator-degree formulas. Nine $B_4$ cochain calculations agree with the orthogonal symmetric-group calculation, including each separate fermion degree. [Character algorithm](../f4_singlet_characters.py) · [All cochain degrees and checks](characters.json)

A8 job `5670599` used all 52 independent derivations of the Albert algebra, acting on its 27-dimensional space of Hermitian octonionic matrices. The scalar identity line contributes zero to each nonempty product. Each generator satisfies the derivation identity on all $27^2$ basis pairs and preserves the trace metric. [Colour matrices](../f4_colour.py) · [Albert-algebra realization](https://arxiv.org/abs/math/0105155)

At $t=4$, evaluating the first classical map on these trace expressions gives a 134-dimensional image at the recorded prime. The full source has dimension $666+79=745$, so its rational classical kernel has dimension at most $745-134=611$. In formal orthogonal traces the outgoing map has rank 358 over the rationals. The direct sum of the finite-colour source evaluation and the formal classical map gives a minor of order 948. It follows that at least $948-358=590$ independent rational closed directions are restrictions of formal classical cycles. The kernel modulo these cycles has dimension at most $611-590=21$. Counting possible lifts after removing these directions uses one-loop protection of universal multi-graviton classes. The quotient estimate itself uses only classical maps. [Projection and determinants](ranks.json) · [A8 receipt](a8_5670599.json) · [First-order protection](https://arxiv.org/abs/2610.22332)

From the repository root, using the formal source plan from the [$SO(9)$ calculation](../rank_four/README.md), run

```sh
python3 algorithms/dual_singlets_20261118/f4_singlet_characters.py --out /tmp/f4-characters.json
python3 algorithms/dual_singlets_20261118/f4_formal_cycles.py --plan /tmp/so9-s5/plan.json.gz --work /tmp/f4-singlets --workers 32
```

## The classical target and its pure-gaugino predecessor

Write $S=\operatorname{Tr}_{27}(\lambda_+\lambda_-)$ for the bilinear in the two gaugino components in the Albert representation. Its scalar identity line contributes zero. At derivative bidegree $(5,5)$, with all three flavour degrees zero, the singlet dimensions and classical ranks over the rationals are

| Number of letters | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Cochains | 1 | 2 | 18 | 62 | 99 | 74 | 31 | 20 | 12 |
| Incoming rank | 0 | 1 | 1 | 17 | 45 | 52 | 22 | 9 | 11 |
| Classical cohomology | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 1 |

Every differential meets a length at which its rank lower bound and that of the adjacent differential sum to the character dimension. The identity $Q_0^2=0$ forces equality at both maps. At length ten, integer evaluations of all 45 formal predecessor images have rank 11. Appending $S^5$ raises the rank to 12. The predecessor evaluations have full cochain rank 20. The saved rational functional vanishes on all 45 images and takes value one on $S^5$. The 96 target evaluations also agree with independent exterior Pfaffians. At length nine the entire classical singlet source has dimension $20-9-11=0$. [Integer target and functional](pure.json.gz) · [Character dimensions and minor certificates](pure_chain.json)

## Exceptional colour contractions in the first-order source

At derivative bidegree $(4,4)$ two open words in the 27-dimensional representation can be paired through either the 52 derivations or the 26 traceless Jordan multiplication operators. The latter span an invariant subspace because a derivation acts on multiplication by $x$ through multiplication by its action on $x$. For each subspace the contraction uses the inverse trace Gram matrix. Reversing an open word includes the transpose sign of its insertion, negative for a derivation and positive for Jordan multiplication. Each sampled parent is a spin and flavour singlet. Its classical variation was formed before evaluating any colours.

| Insertions | Rank detected in the length-eleven image | Lower bound on formal cycles plus nine-letter boundaries | Upper bound on the remaining length-ten kernel |
|---|---:|---:|---:|
| 52 derivations | 134 | 590 | 21 |
| 26 Jordan multiplications together with the ordinary trace rows | 139 | 599 | 7 |

The outgoing determinant of order 139 is 345573 modulo 1000003. The determinant of order 957 for the source evaluation augmented by the formal differential is 118669 at the same prime. The formal differential has rational rank 358. Each new nine-letter boundary is entered with zero formal coordinate. The rational dimension of the resulting closed span is at least $957-358=599$. Hence its codimension in the complete classical kernel is at most $745-139-599=7$. Direct source evaluations attain rank 738 of 745. These bounds use classical maps. A quantum pairing requires evaluating the physical first differential on the surviving quotient. [Outgoing minor](jordan_outgoing.json) · [Incoming minor](jordan_incoming.json) · [A8 job 5671196](a8_5671196.json) · [Adjoint outgoing test](adjoint_outgoing.json) · [Adjoint incoming test](adjoint_incoming.json)

From the repository root, using the source plan and ranks above, run

```sh
python3 algorithms/dual_singlets_20261118/f4_pure_target.py --out /tmp/f4-pure.json.gz --workers 4
python3 algorithms/dual_singlets_20261118/f4_pure_chain.py --pure /tmp/f4-pure.json.gz --work /tmp/f4-pure-chain
python3 algorithms/dual_singlets_20261118/f4_joined_boundaries.py --plan /tmp/so9-s5/plan.json.gz --previous /tmp/f4-singlets/ranks.json --work /tmp/f4-jordan --channel jordan --count 450 --workers 48
python3 algorithms/dual_singlets_20261118/f4_joined_incoming.py --plan /tmp/so9-s5/plan.json.gz --previous /tmp/f4-singlets/ranks.json --outgoing /tmp/f4-jordan/jordan_outgoing.json --work /tmp/f4-jordan --channel jordan --count 480 --source-count 300 --workers 48
```

## A necessary and sufficient classical equation for the first differential into $S^5$

At derivative bidegree $(4,4)$ the 666 seven-fermion directions have three separate scalar flavours. The first-order pair rule requires all three flavours at two sites, hence vanishes. In the other 79 directions the two distinguished letters are a scalar and a matter fermion. The remaining eight letters are gauginos. Projecting formal cycles and nine-letter boundaries to this component gives rank at least $436-358=78$ over the rationals. The determinant of order 436 is 710453 modulo 1000003. The full nine-fermion evaluation reaches rank 79. [Projection](pair_support.json) · [A8 job 5671321](a8_5671321.json)

To evaluate the coefficient of $[S^5]$ after acting on an ordinary trace word, insert a colour generator and its trace dual through the two commutators. In pure terminal cohomology the two gauginos anticommute, and each squares to zero. Cutting at the two insertions leaves at most two open strings of gauginos, each of length at most two. A nonvanishing term therefore has at least three spectator bilinears. After removing three, the four-letter classical singlet has dimension $2-1=1$. Two commuting colour generators give top exterior coefficient 4608 on $S^2$. Division by this value extracts the normalized coefficient. The pure fifth power is independently nonzero.

Of the 1,904 ordinary trace singlets, only rows 717, 727 and 1348 have nonzero first-action coefficients, respectively $-18,1,-1/2$. On every row this coefficient equals three times formal classical image coordinate 11 plus one half of each of coordinates 47 and 236. All coordinates are zero-based. All 1,904 residuals vanish over $\mathbb Q$. The first action therefore vanishes on the restrictions of all formal cycles. It vanishes on classical boundaries in cohomology by $Q_1Q_0=-Q_0Q_1$. [Coefficients and rational row identity](first_functional.json) · [A8 job 5671727](a8_5671727.json)

For the three-term cochain in the saved canonical trace order, $Q_1$ is $18 S^5$ and $Q_0$ is $6S^4\operatorname{Tr}_{27}(\phi_1\phi_2\phi_3)$. The first-action functional is consequently nonzero on the 79-dimensional component. The preceding 78-dimensional projected closed span is its entire kernel. Subtracting the corresponding closed representatives from any putative preimage and rescaling reduces existence of a first-order pairing to the single equation

$$
Q_0 X=-6S^4\operatorname{Tr}_{27}(\phi_1\phi_2\phi_3),
\qquad X\text{ in the seven-fermion singlet of dimension }666 .
$$

Conversely a solution added to the saved cochain has zero classical action and first action $18[S^5]$. A8 job 5671407 reached seven-fermion source rank 659 and classical image rank 138. Adjoining the prescribed right-hand side gives rank 139 at the same prime. The unrestricted rational equation is open. [Cochain and integer actions](scalar_source.json) · [Minors](scalar_correction.json) · [Slurm receipt](a8_5671407.json)

```sh
cd algorithms/dual_singlets_20261118
python3 f4_first_source.py --out /tmp/f4-scalar.json
python3 f4_pair_support.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --jordan f4/jordan_incoming.json --out /tmp/f4-pair
python3 f4_scalar_correction.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --work /tmp/f4-correction --count 420
python3 f4_first_functional.py --plan /tmp/so9-s5/plan.json.gz --previous f4/ranks.json --work /tmp/f4-first
```
