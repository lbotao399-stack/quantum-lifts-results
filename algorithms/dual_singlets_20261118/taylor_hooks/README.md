# Classical cohomology of a single Taylor representation

The two even Taylor directions are the two spacetime derivatives. Retain one odd flavour direction. A partition with rows $(s,s,1^c)$, where $s\ge1$ and $1^c$ means $c$ rows of one box, already occurs on this smaller space. Its multiplicity can be recovered there because the classical differential $Q_0$ is natural in the Taylor directions: polarizing the jets commutes with splitting one jet into two. At total degree $D$ the multilinear cochain space has a permutation action on its $D$ Taylor slots. The Schur projection has denominator dividing $D!$. These computations use primes larger than $D$.

Let $x_1,x_2$ record the even degrees and $y$ the odd degree. For rows $(u,v,1^k)$ with $v\ge1$, the hook character is

$$
s_{(u-1,v-1)}(x_1,x_2)y^k(x_1+y)(x_2+y).
$$

Working in descending powers of $y$, the inverse of $(1+x_1/y)(1+x_2/y)$ has coefficient $(-1)^i h_i(x_1,x_2)$ at $y^{-i}$, where $h_i$ is the character of the $i$th symmetric power. To extract the even rectangle $(s-1,s-1)$, the Pieri rule removes a horizontal strip of $i$ boxes from its lower row. Thus at letter length $Y$ the multiplicity equals

$$
\sum_{\substack{0\le i<s\\\text{horizontal strips}}}(-1)^i\,
\dim H(Q_0)_{Y,\,(c+2+i,0,0;\,s-1,s-1-i)} .
$$

The tuple gives the three flavour degrees followed by the two even degrees. A one-row diagram has at most one odd box and contributes nothing in this expansion at $y^{c+2}$. Expanding instead about $y=0$ gives the independent sum with terms $(-1)^i\dim H(Q_0)_{Y,(c-i,0,0;s+i,s)}$, for $0\le i\le c$. The same identities apply to the cochains and to each classical image. The calculation takes every term in one field before extracting the Schur component.

| Group and target | $Y$ | Cochains | Incoming rank | Outgoing rank | Classical cohomology at 2147483647 |
|---|---:|---:|---:|---:|---:|
| $SU(3)$, $(5,5,1^3)$ | 9 | 7268 | 4267 | 2983 | 18 |
| $USp(6)$, $(4,4,1^3)$ | 8 | 1316 | 538 | 763 | 15 |
| $USp(6)$, $(4,4,1^3)$ | 9 | 1198 | 763 | 434 | 1 |
| $USp(6)$, $(4,4,1^3)$ | 10 | 491 | 434 | 57 | 0 |
| $SO(7)$, $(4,4,1^3)$ | 8 | 1317 | 538 | 764 | 15 |
| $SO(7)$, $(4,4,1^3)$ | 9 | 1201 | 764 | 437 | 0 |
| $SO(7)$, $(4,4,1^3)$ | 10 | 495 | 437 | 58 | 0 |

The second-order target of the $SU(3)$ level-30 BMN source gives $110-198+163-72+15=18$. The calculation combines both charge-conjugation parities and imposes the linear Taylor symmetry. Further conditions for a physical partner include charge conjugation and the other generators of the superconformal commutant. Obtaining a differential requires its action on the source. [Five rank records](su3_second.json) · [Prior superconformal analysis](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex)

At $USp(6)$ level 28 the nine-letter class lies in this diagram. Its ten-letter successor in the same diagram is empty. The multiplicity calculation is independent of applying the Ward generators to its 4,264-term representative. [Symplectic calculation](hwomega/resolved_Sp3_4_3.json) · [Orthogonal calculation](hwomega/resolved_SO7_4_3.json) · [Direct Ward test](../sp3_omega/outgoing/ward_source.json)

For the $SU(3)$ level-30 source the eighteen potential $d_2$ targets are all odd under colour charge conjugation in both sampled fields. Splitting the classical map by the evaluation at a traceless matrix and at its negative transpose gives the following multiplicities.

| Colour parity | Cochains | Incoming rank | Outgoing rank | Classical cohomology at both primes |
|---|---:|---:|---:|---:|
| Even | 3585 | 2093 | 1492 | 0 |
| Odd | 3683 | 2174 | 1491 | 18 |

The source is odd, so all eighteen remain possible at this stage. A8 jobs 5682907 and 5682908 used different seeds and primes 2147483647 and 1000033. [Parity-resolved sectors and receipts](colour_parity/) · [Native map](unitary_parity.cpp) · [Schur extraction](../unitary_hook_parity.py)

For the even $SU(3)$ level-34 BMN source of weight $(8,5,4)$ the order-two diagram has rows $(5,5,2,1,1,1)$. Restriction to two even and two odd Taylor directions gives 53,037 even cochains and 52,919 odd cochains. Independent permutation characters give the same integers. The original even target has 639,219 cochains. At third order the shifted Taylor highest weight would be $(6,6;3,0,-1)$. Since the last degree is negative, the third-order target is absent under filtered superconformal equivariance. The first and second correction equations for this source remain open. [Two-flavour branching and colour characters](../two_flavour_hooks.py) · [Twenty signed character terms](su3_j34_second_character.json)

At $USp(4)$ levels 40 and 42 the third-order diagrams for BMN weights $(10,5,5)$ and $(11,5,5)$ have rows $(6,6,1^5)$ and $(6,6,1^6)$. At the maximal target lengths 17 and 18, A8 job 5682624 gives $43-43-0=0$ and $124-124-0=0$ respectively in the field of 2147483647 elements. Every modular rank is taken on complete native weight spaces before extracting the Schur component. The corresponding rational target cohomologies vanish. Under filtered Taylor equivariance this excludes the transitions from source lengths 19 and 20. Shorter sources have different targets. [Twelve weight-space ranks, resolved blocks and receipt](b2_third_terminal/)

Under the second-order Taylor selection, the rational targets of dimensions 18, 28 and 61 at $SU(3)$ levels 30 and 32 and $USp(4)$ level 36 are all detected by Cartan restriction. Their incoming second differential has rank zero by the discriminant argument. [All-orders argument](../../perturbative_cartan_targets.md) · [Product spans and replay](cartan_products/README.md)

Compile the archived classical engine and run the same sectors from the repository root. The engine's level label refers to the substituted Taylor weight. The physical level of the original operator is given above.

```sh
d=algorithms/dual_singlets_20261118
c++ -O3 -std=c++17 -pthread "$d/taylor_hooks/native/engine.cpp" -o /tmp/n4-hook-engine
python3 "$d/taylor_hook_cohomology.py" --engine /tmp/n4-hook-engine --group SU3 --spin 5 --tail 3 --lengths 9 9 --threads 12 --out /tmp/su3-hooks
python3 "$d/taylor_hook_cohomology.py" --engine /tmp/n4-hook-engine --group Sp3 --spin 4 --tail 3 --lengths 8 10 --threads 12 --out /tmp/sp3-hooks
```
