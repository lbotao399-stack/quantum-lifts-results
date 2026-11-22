# Scalar multipliers of the orthogonal third differential

Use the defining trace $S=\operatorname{Tr}(\bar\lambda_{\dot +}\bar\lambda_{\dot -})$ and the campaign's ordered-pair normalization. Write $A=\phi^1$. The three trace coordinates with $p$ gauginos of each sign and two scalars, in order, are

$$
S^p\operatorname{Tr}(A^2),\quad
S^{p-1}\operatorname{Tr}(\bar\lambda_{\dot -}A)\operatorname{Tr}(\bar\lambda_{\dot +}A),\quad
S^{p-1}\operatorname{Tr}(\bar\lambda_{\dot +}\bar\lambda_{\dot -}A^2).
$$

All groups here are odd orthogonal unless the Pfaffian is explicitly specified. The full relative holomorphic jet complex includes the second-jet predecessors. Each terminal coordinate uses only first jets. The physical interpretation uses the comparison in Appendix B of the perturbative manuscript. [Supercharge and product](https://arxiv.org/abs/2611.24160)

| Group and terminal weight, dotted then flavour | Classical dimension | Incoming $d_1$ rank | Result on later pages |
| --- | ---: | ---: | --- |
| $SO(11)$, $(7,6;1,0,0)$ | 1 | 1 | Target exhausted |
| $SO(11)$, $(6,6;2,0,0)$ | 3 | 3 | Target exhausted |
| $SO(11)$, $(6,6;4,0,0)$ | 6 | 6 | Target exhausted |
| $SO(13)$, $(7,7;2,0,0)$ | 3 | 1 | $S^7\operatorname{Tr}(A^2)$ survives to $E_3$ and is hit by $d_3$ |

These four rows add eleven first-differential multiplets under the campaign's counting convention. Every target has a single flavour highest direction. For fixed dotted parity, each dotted-spin multiplet contains exactly one balanced component (an excess of one plus index for half-integer spin). Products and descendants in distinct target multiplets are included. The pure-power $SO(11)$ first map and the $SO(13)$ third source come from earlier work. This calculation also gives a $d_3$ into $S^7\operatorname{Tr}(\bar\lambda_{\dot+}A)$ in $SO(13)$. Both mixed third targets use Dot's known source.

## Exact first images

Use the polarized orthogonal Pfaffian and the tree primitive $h$ which reverses an odd trace by adjacent graded interchanges. Every first relation in this calculation has odd total length. Its stable tree primitive has tree derivative equal to the entire relation. At finite colour the relation vanishes, so the primitive specializes to a cocycle. Pulling a terminal coordinate back through the ordered $Q_1h$ evaluates its first image. The verifier checks every potentially contributing column of the transposed identity. [Pfaffian, colour conventions and orientation](../README.md) · [Odd-trace primitive](../../orthogonal_odd_trace_primitive.py)

For eleven colours the three Pfaffians specified literally in the verifier have degrees $12,13,12$. Their first pairings in the three displayed coordinates are

$$
\begin{pmatrix}40&-200&-880\\-40&80&640\\0&0&30\end{pmatrix},
\qquad \det=-144000 .
$$

Sorting first jets in an even orthogonal trace gives the graded-commuting trace, and sorting an odd-length trace gives zero. At this terminal weight there is no finite relation. In the orthogonal determinant with twelve left and twelve right endpoints there are only $14-12=2$ ordinary pairings. Six identical-type odd single paths, three per dotted sign, are the maximum. Two additional edges can cover at most four more odd letters. The twelve gauginos exceed this bound. The one-scalar target has thirteen gauginos at the same length and also has no finite relation. [Invariant tensor theorem](https://arxiv.org/abs/1102.3221)

At four scalars the stable sorted space has dimension nine. Its first three coordinates are the preceding three multiplied by $\operatorname{Tr}(A^2)$, so the page product makes them incoming first images. To find the remaining directions, replace one first-jet matrix $X$ in a degree-fifteen relation by $AXA$ or $A^2X+XA^2$ and sum over its occurrences. These substitutions preserve skew transposition. Differentiating the original finite-colour identity under this polynomial change of matrix variables gives another finite-colour identity. Apply the tree primitive and evaluate the first differential **after** the substitution. No covariance of the quantum kernel under nonlinear changes of variables is required.

Six such sources have determinant $843448320000000$ on the six remaining stable coordinates. Before imposing finite relations this fills every target coordinate modulo the preceding three. The complete classical finite-identity image has rank three. Thus the physical first image has dimension $9-3=6$. The eight reduced path types are enumerated from the four ordinary edges of a determinant with twelve endpoints on each side. Open paths sort by the ordinary first-jet tree commutators. Identical odd letters vanish within a sorted path. A closed cycle must be even. The only closed spectator that occurs is the two-gaugino bilinear. For an open path the two auxiliary endpoints are symmetric at odd length and antisymmetric at even length. This gives the three traceless axes or the single identity axis, respectively. Identical odd single paths exteriorize, and the scalar singletons are symmetric. Taking all zero-weight auxiliary polarizations gives 49 rows. For each row, ordered trace expansion followed by sorting agrees with direct coefficient extraction from the commuting trace logarithm in all nine coordinates. The three independent rows and the six-column first matrix are in the replay. [Exact matrices](first.json) · [Determinant path reduction](https://github.com/lbotao399-stack/quantum-lifts-results/tree/f05b289fc69fb5f65bd8286c09060156a554f23d/reviews/orthogonal-mixed-lifts)

For thirteen colours the $(7,7;2,0,0)$ terminal again has no finite relations. At the first dual weight the seventeen first jets include six gauginos of each sign. Removing the fourteen determinant edges leaves three ordinary edges. After sorting, all three must join a plus and a minus gaugino, and the other six gauginos must occupy their six distinct single paths. The five scalar singletons project to one scalar-marker determinant times a quadratic of flavour one. Its zero-weight auxiliary axes are $H,H$ and $E,F$. Their two first images are $(-60,120,1080)$ and $(-30,60,540)$. Sorting inside a relation is valid for this first pairing: the change is $Q_0$ of a finite-colour identity, and its $Q_1$ image is a finite-colour identity at the terminal weight. Stable sources give zero by the verified first transposed equation. Every finite-colour cycle lifts to a stable cochain with tree derivative a relation. This proves the rank-one upper bound as well as its attainment.

## Returning the product to a known third-page class

The strict odd substitution $P_\pm$ changes the flavour-one scalar first jet to $\bar\lambda_{\dot\pm}$ and annihilates the gaugino first jets. On a divided higher jet it removes flavour one, multiplies by one plus the existing Taylor degree in the indicated direction, and includes the exterior reordering sign. On every possible earlier chain here the total Taylor degree exceeds the source length by one. There is exactly one second jet. Thus a bare $Q_r$, $r\ge2$, vanishes by its $3r$-flavour requirement. For $Q_1$ the only active pair is a matter second jet and its complementary scalar. Substitution on either input destroys the pair and substitution on either gaugino output vanishes. On spectators the two odd actions anticommute. The verifier checks the classical action on all seventeen possible letters and the quantum condition on all 145 ordered pairs for each sign. This lets the substitutions act on the earlier-page closure and boundary equations in the full sector. [Native substitution check](raising.py) · [Prior strict-action verification](https://github.com/lbotao399-stack/quantum-lifts-results/blob/21ca40461b1bcacdf591f5f2fd85d85a5db6f8f0/reviews/orthogonal-mixed-lifts/bmn13/THIRD.md)

Applying $P_-$ and then $P_+$ to the three coordinates gives $(2,1,0)S^{p+1}$. Applying $P_-$ to $S^p\operatorname{Tr}(\bar\lambda_{\dot+}A)$ gives $-S^{p+1}$. At $p=6$ the first $SO(11)$ quadratic image already spans the chosen quadratic product, so one substitution also gives the entire one-scalar first image. At $p=7$ in $SO(13)$, $S^8$ survives to $E_3$ by the full 301-relation certificate. Both displayed mixed products therefore survive $d_1,d_2$. The two bilinear multipliers are closed gravitons. Multiplication on the third page by either of them sends Dot's $d_3$ source into the corresponding mixed product with the same coefficient $-921600$ when the even/odd factors are ordered source first. This uses the page Leibniz rule. [Pure third source](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md) · [Pure eighth-power lower-page certificate](../so13_j48/README.md)

The remaining potential $d_2$ image in the thirteen-colour quadratic block is confined to one line. It lies in the kernel of $(2,1,0)$ modulo the first image $(1,-2,-18)$. The independent covector $(18,0,1)$ tests that line. The stable first row has 71 words. Its second row has 2,809 words and satisfies the transposed identity on 15,928 columns.

Three families of finite identities give zero on this remaining line. From the three pure seventh-power relations insert two scalars in a single letter as above. From the pure eighth-power relations replace one gaugino of each sign by flavour one using the two ordered classical odd derivations. Finally use each of 76 earlier mixed polarizations and replace one scalar in turn by the plus gaugino, that same scalar, and the minus gaugino in matrix order, minus the reversed word. The composite is skew; the two odd entries contribute one extra sign under transposition. All three constructions preserve finite-colour identities. On a relation whose stable commuting section is zero, the odd-length intermediate output has a tree primitive. Applying the second row to the relation therefore gives the physical second pairing.

For the pure eighth-power transfer, the new score equals $-144$ times the archived pure second score on all 56 reduced polarizations. The two classical dotted-to-flavour derivations anticommute with the tree charge and preserve the finite-identity ideal, so this test extends through the previously proved reduction of all 301 pure relations modulo tree variations of identities. The three independent original commuting sections map to three independent mixed sections with minor $103680$. The first two families together have commuting rank eight and augmented rank eight. [Pure reduction and basis](../so13_j48/README.md) · [Second-row verifier](second.py)

On the 76 mixed polarizations the three new scores factor through the five earlier commuting coordinates at indices $0,1,2,10,11$, with pivot minor $19237055439568896000$. For the union of all three families the rank is 23 and appending the physical pairing leaves it 23. An explicit rational functional using 23 coordinates of the sorted mixed traces reproduces the pairing on the combined span. This excludes linear combinations across the three families. Other finite identities at the same weight remain open. [Base mixed rank proof](https://github.com/lbotao399-stack/quantum-lifts-results/tree/21ca40461b1bcacdf591f5f2fd85d85a5db6f8f0/reviews/orthogonal-mixed-lifts/bmn13) · [Joint verifier](joint.py) · [Polarizations and rational functional](mixed.json)

## Pfaffian multipliers

In the classical first-jet quotient put $X=\bar\lambda_{\dot+}\bar\lambda_{\dot-}$. The second-jet tree variations give graded commutation of first jets, including the squares of both gauginos. Hence $X^2=0$, $X^t=X$, $X$ commutes with every first-jet matrix, and $S=\operatorname{Tr}(X)$. These are matrix identities in the quotient by covariant tree variations. Contracting any such variation with invariant tensors and other first jets yields a gauge-invariant tree boundary.

In $SO(2r)$ every reflection-odd first-jet contraction reduces to one volume tensor, $r$ paths between its slots, and ordinary trace factors. The infinitesimal action of $X$ on the volume tensor is $\operatorname{Tr}(X)$ times that tensor. Inserting $X$ at either end of a given path has the same effect by symmetry and commutation. Applying the volume identity $r+1$ times puts two insertions on at least one path. Their product is zero. Thus $S^{r+1}$ annihilates **every** reflection-odd first-jet invariant already classically. In particular a reflection-odd first-jet multiplier of the known $SO(12)$ target $S^7$ cannot expose a second differential by the annihilator test.

For $SO(2r+1)$ take the standard vector of signed principal Pfaffian minors of a single skew scalar $A$. Its squared norm is the degree-$2r$ coefficient of the characteristic polynomial. Each gaugino annihilates this vector: equivariance differentiates the skew conjugation action, and its classical commutator with $A$ is zero. When the trace identity is applied to either volume tensor in the squared norm, the free slot therefore gives zero, and the $r$ paired slots again accept at most one $X$ each. This proves that the product of $S^{r+1}$ with this top characteristic coefficient is classically zero. It gives a separate exclusion in eleven and thirteen colours.

## Replay

Use the checkout's trace runtime with NumPy, SciPy and python-flint. From this directory run

~~~sh
python3 verify.py > /tmp/orthogonal-multiplier-first.json
cmp /tmp/orthogonal-multiplier-first.json first.json
python3 raising.py
python3 joint.py
~~~

Pass `--rebuild --workers 2` to `joint.py` to regenerate all 56 pure and 76 mixed second-pairing rows. The default regenerates the physical dual, the spanning-image pivots, the three small pure relations and the joint functional.
