# A third supercharge differential into the sixth glueball power of SO(11)

In the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, let $S=\operatorname{Tr}_{11}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$ in the defining orthogonal representation. At level $j=2(E+J_L)=36$ its sixth power survives to the third page. The source below satisfies

$$
d_3[O_0]=8640[S^6]\ne0 .
$$

The coefficient uses the archived rational one-loop kernel. The comparison of this perturbative complex with physical local operators is the one in Appendix B of the all-loop manuscript. Dot's target-first approach motivated this calculation. [One-loop kernel](../../../su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-protected/q1d_exact.py) · [All-loop supercharge and comparison](https://arxiv.org/abs/2611.24160)

## Source and corrections

Use the polarized Pfaffian and the same auxiliary matrices, marker order and trace orientation as in the determinant sources. Take six single gaugino paths, one for each of the two dotted components and three auxiliary matrices. Apply the cube of the scalar-marker determinant of derivatives, which selects nine single scalar paths. The resulting coefficient has Pfaffian degree $15$. The Pfaffian of a $22$-by-$22$ matrix has degree $11$, so this is a finite-colour relation. Its stable trace expansion is fixed by $\tfrac12\operatorname{Tr}\log(1+A)$ in the auxiliary matrix $A$ of the source prescription. Every summand has an odd-length trace. Reversal of the first such trace, with each adjacent graded commutator replaced by its second-jet primitive, gives $O_0$. The factor $1/2$ comes from the sign of the reversed odd trace. Its tree differential is the degree-$15$ relation. [Pfaffian and orientation](../README.md) · [Reversal primitive](../../orthogonal_odd_trace_primitive.py)

| Cochain | Dotted degrees | Flavour degrees | Length |
|---|---|---|---:|
| $O_0$ | $(3,3)$ | $(3,3,3)$ | 14 |
| $O_1$ | $(4,4)$ | $(2,2,2)$ | 13 |
| $O_2$ | $(5,5)$ | $(1,1,1)$ | 12 |
| $S^6$ | $(6,6)$ | $(0,0,0)$ | 12 |

Each correction has precisely one second jet and otherwise first jets. The order-$r$ fan requires $3r$ flavour slots on $r+1$ inputs. On these cochains they have at most $r+2$ slots, so every bare $Q_r$ with $r\ge2$ vanishes.

Write $h$ for the explicit tree primitive with $Q_0h(w)=w-\mathrm{section}(w)$. On an even trace it sorts the letters by graded interchanges and removes repeated odd letters by the second ghost jet. On a product containing an odd trace its section is zero. Products of sorted even traces give the stable commuting quotient. Evaluating them on arbitrary numbers of distinct Cartan planes gives an invertible change of its 257 weight coordinates. [Tree homotopy](../../orthogonal_glueball_second_pair.py) · [Change of basis](../../orthogonal_glueball_cartan_basis.py)

In the zero-based order of the 56 determinant rows, choose rows $1,9,18$. Their first three distinct-plane coordinates and second pairings are

$$
\left(\begin{array}{rrr|r}
9216&-3072&0&288\\
9216&-3072&384&192\\
0&768&-384&-48
\end{array}\right).
$$

Let $R$ be $135$ times row $1$, minus $90$ times each of rows $9$ and $18$, as **ordered trace polynomials**. Each summand has Pfaffian degree $12$ and vanishes in eleven colours. All 257 commuting coordinates of $R-Q_1O_0$ vanish before finite-colour specialization. The full ordered corrections are therefore specified by

$$
O_1=h(R-Q_1O_0),\qquad O_2=-h(Q_1O_1).
$$

Every term of $Q_1O_1$ has odd length $13$, so the second formula has zero section term by term. After specialization,

$$
Q_0O_0=0,\qquad Q_0O_1=-Q_1O_0,\qquad Q_0O_2=-Q_1O_1 .
$$

The three determinant rows have scalar paths specified by their individual recipes. Rebuilding each row by sorting ordered traces and changing to distinct-plane coordinates gives the same 257 coefficients as evaluating its Pfaffian on the planes directly. The $3$-by-$3$ minor is $-2717908992$. A8 job `5693264` reconstructed the first action from 181596 quantum parents, with 130488 contributing, and matched the forward vector in all 257 coordinates. [First-action replay](first.json.gz) · [Determinant rows](second.json.gz) · [Ordered coefficient extraction](../../orthogonal_glueball_second_sources.py)

## Why the target survives both earlier differentials

At pure gaugino weight $(6,6)$ the stable tree quotient is one-dimensional. A sorted nonempty trace can use each odd letter only once and has even length, leaving just the bilinear. The complete finite-colour relation multiplicities at this weight and the preceding weight $(5,5;1,1,1)$ are both zero. Thus $[S^6]\ne0$ classically. All its gaugino letters have flavour degree zero, so every positive-order bare supercharge kills this representative.

Let $\lambda$ extract $[S^6]$. The first and second dual rows in the replay obey, on stable ordered traces,

$$
\lambda Q_1=\mathrm{first}\,Q_0,\qquad
-\mathrm{first}\,Q_1=\mathrm{second}\,Q_0 .
$$

The first row uses 13 ordered words, the second uses 251. Tree splitting and inverse one-loop reconnection enumerate all columns capable of entering either side. Two enumerations verify the identities, using respectively $(33,1266)$ and $(39,1362)$ columns. Extra zero columns are included in the latter test. The absence of relations at the first dual weight lets the first row descend directly to eleven colours.

From the second row subtract the functional on the stable commuting quotient whose first three distinct-plane coefficients are $(-1/32,-3/16,-1/4)$ and whose remaining coefficients vanish. This preserves the second adjoint equation since a commuting coordinate vanishes on tree boundaries. The corrected row vanishes on **all 29** finite-colour relations with balanced dotted weights and singlet flavour. The orthogonal invariant-tensor character gives 29 before imposing the dotted singlet. Independently, 29 explicit degree-$12$ or higher Pfaffian relations have the following nonzero determinant on ordered-word probes,

$$
5983560995463773616373296508630663618556816636006891520 .
$$

Each relation is flavour invariant by the square of the scalar derivative determinant. The auxiliary identity multiplies a symmetric colour path and each traceless auxiliary matrix multiplies an antisymmetric path, so every recipe is a coefficient of an actual skew-matrix Pfaffian. Spin averaging is unnecessary for this larger relation space. Gauge invariants of adjoint letters are identical for $SO(11)$ and $O(11)$ since the central element $-1$ lies in the other component. Flavour projection commutes with the supercharge. Thus the corrected row descends on every possible earlier source contributing to the singlet target. The two adjoint identities annihilate every $d_1$ and $d_2$ image on that target. [Character, Pfaffians and exact test](../../orthogonal_glueball_relation_span.py) · [Basis and all 29 residuals](relations11.json.gz)

## Third pairing

The second adjoint applied to $Q_0O_1=R-Q_1O_0$ gives

$$
\begin{aligned}
\lambda Q_1O_2
&=-\mathrm{first}\,Q_1O_1
=\mathrm{second}\,(R-Q_1O_0)\\
&=-17280+135(288)-90(192)-90(-48)=8640 .
\end{aligned}
$$

The forward pass tested 109525 ordered leading words. The transposed pass uses the corrected dual row directly. It cuts and reconnects the two colour edges, reconstructs the preimage of each merged scalar pair under reversal, and evaluates the ordered Pfaffian coefficients. Its 18035 quantum parents give 12362 nonzero contributions and the same sum. Both passes use the same physical one-loop trace kernel. The nonzero pairing on the third page also proves that $[O_0]$ survives to that page. [Forward A8 calculation](forward.json.gz) · [Transposed A8 calculation](adjoint.json.gz) · [Inverse reconnection and reversal](../../orthogonal_glueball_adjoint_replay.py)

At $SO(13)$ the same lower-page calculation gives three flavour-singlet relations at $(5,5;2,2,2)$, with minor $9289728$. The corrected second row vanishes on all three. The seventh glueball power at level 42 survives to the third page. [SO(13) certificate](relations13.json.gz)

## Replay

With Python 3.12 or later, NumPy, SciPy and python-flint, run from this directory. The existing repository supplies the trace runtime.

~~~sh
python3 verify.py --rebuild-relations --workers 4
python3 replay_forward.py --workers 8
python3 ../../orthogonal_glueball_adjoint_replay.py second.json.gz --out /tmp/so11-third.json --workers 8
python3 ../../orthogonal_glueball_adjoint_replay.py second.json.gz --out /tmp/so11-first.json --workers 8 --first
~~~
