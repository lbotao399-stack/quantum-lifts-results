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
