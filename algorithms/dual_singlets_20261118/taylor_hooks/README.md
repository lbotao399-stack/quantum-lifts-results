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

Compile the archived classical engine and run the same sectors from the repository root. The engine's level label refers to the substituted Taylor weight. The physical level of the original operator is given above.

```sh
d=algorithms/dual_singlets_20261118
c++ -O3 -std=c++17 -pthread "$d/taylor_hooks/native/engine.cpp" -o /tmp/n4-hook-engine
python3 "$d/taylor_hook_cohomology.py" --engine /tmp/n4-hook-engine --group SU3 --spin 5 --tail 3 --lengths 9 9 --threads 12 --out /tmp/su3-hooks
python3 "$d/taylor_hook_cohomology.py" --engine /tmp/n4-hook-engine --group Sp3 --spin 4 --tail 3 --lengths 8 10 --threads 12 --out /tmp/sp3-hooks
```
