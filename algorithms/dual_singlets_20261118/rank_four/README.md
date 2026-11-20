# The fifth gaugino power in $SO(9)$ and $USp(8)$

At level $j=2(E+J_L)=30$, put $S=\operatorname{Tr}(\lambda_+\lambda_-)$ with the defining trace and the two gaugino components of the source convention. The terminal singlet $S^5$ has ten letters and doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(0,0,10,10,10)$. Its first-order source has charges $(1,1,10,10,10)$.

| Terminal singlet | Formal dimension | Finite-colour dimension | Classical boundary rank | Rank with $S^5$ |
|---|---:|---:|---:|---:|
| $SO(9)$ | 20 | 20 | 19 | 20 |
| $USp(8)$ | 20 | 19 | 19 | 19 |

The formal boundary rank is 19 over the rationals. Each nine-letter classical primitive contains at most four traces. Replacing one letter by its classical variation preserves trace number, so the coefficient of the product of five bilinears vanishes on its image. The finite-colour cochain dimensions follow from independent singlet characters. The integral evaluations attain the indicated dimensions and ranks. Thus extraction of the canonical five-bilinear coefficient, divided by its sign in $S^5$, defines a rational functional on the entire $SO(9)$ target. It evaluates to one on $S^5$. [Target minors](pure.json)

At derivative bidegree $(4,4)$ and flavour degrees $(1,1,1)$, the ten-letter source has 1,904 formal singlets. Its two finite-colour fermion degrees contribute $1707+197=1904$. The classical supercharge $Q_0$ maps them to a singlet space of dimension 357 at length eleven. One first-order contraction changes the number of traces by at most one. Of 52,608 words appearing in the source basis, 1,576 have at least four traces. Exact pair contraction on these words gives only three nonzero singlet coordinates for the normalized fifth-power coefficient of $Q_1$. Their zero-based row numbers and coefficients are $(717,-84),(727,12),(1348,-6)$. [Recoupling rule](../../su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-protected/q1d_exact.py)

| Prime | Classical minor order | Determinant | Minor order after adjoining the $Q_1$ coefficient | Determinant |
|---:|---:|---:|---:|---:|
| 1000003 | 357 | 739063 | 358 | 973909 |
| 1000033 | 357 | 375392 | 358 | 637559 |

Each modular minor lifts to a nonzero rational minor by lifting the independent entries of the antisymmetric matrices. The classical rank upper bound is 357. The augmented map therefore has rational rank 358, and some rational vector in the classical kernel has a nonzero fifth-power coefficient under $Q_1$. The identity $Q_0Q_1+Q_1Q_0=0$ makes this coefficient vanish on classical boundaries in the source. The one-dimensional target cohomology then gives

$$
\operatorname{im}d_1=\operatorname{span}_{\mathbb Q}\{[S^5]\}\quad\text{in this terminal }SO(9)\text{ singlet},\qquad E_2=0.
$$

The source's independence of previously counted products and descendants is open. [All pivot indices, determinants and quantum coordinates](first.json) · [A8 job 5670470](a8_5670470.json)

From the repository root, with the Python requirements of the independent-prime runtime and a C++17 compiler, run

```sh
python3 algorithms/dual_singlets_20261118/rank_four_s5.py --work /tmp/so9-s5 --workers 24
cmp /tmp/so9-s5/pure.json algorithms/dual_singlets_20261118/rank_four/pure.json
cmp /tmp/so9-s5/first.json algorithms/dual_singlets_20261118/rank_four/first.json
```

## Rational sources at levels 30 and 36

A8 job `5670559` reconstructed the level-30 source $F$ with $[Q_1F]=[S^5]$. Its 11,891 coefficients have common denominator 12. Its first action contains just two trace monomials. With the three-term primitive $H$ the equality $Q_1F-Q_0H=S^5$ holds coefficient by coefficient before imposing finite-colour trace identities. [Source](source.json.gz) · [Primitive](primitive.json.gz) · [Coordinates on the 323 formal singlets used](reconstruction.json) · [Slurm receipt](a8_5670559.json)

To check $Q_0F=0$ in nine colours, use one fixed frame of integral antisymmetric matrices with entries $0,1,-1$. It has 405 evaluation coordinates and detects the full classical image of rank 357. After clearing the reconstruction denominator, the sum of the absolute formal coefficients in the classical image is 1097880. Each length-eleven trace product has at most $8!\,9^{11}$ summands per exterior coefficient. The four primes in the certificate give zero in every coordinate, and

$$
1000112004278059472142857
>2\cdot1097880\cdot8!\,9^{11}
=2778260706025372108800.
$$

The integer evaluations therefore vanish. Injectivity on the classical image gives the claimed equality of finite-colour polynomials.

The two gaugino letters in $S$ have no flavour derivative. Every letter in $F$ has at most two. The first-order pair rule needs all three flavours, so every contraction joining $S$ to $F$ vanishes. All 72 ordered pair checks pass. The gauginos are classically closed and $S$ is even. Hence

$$
Q_0(SF)=0,\qquad Q_1(SF)-Q_0(SH)=S^6 .
$$

At length twelve, 55 signed-partial-permutation evaluations give a faithful integer frame on the 47-dimensional pure $SO(9)$ singlet. The 146 predecessor rows have rational image rank 46. Adjoining $S^6$ gives rank 47. Its 55 coordinates also equal $6!$ times the Pfaffians of the corresponding skew trace-pairing matrices. Thus $d_1[SF]=[S^6]\ne0$ at level 36. The two results count distinct singlet pairings. The second source is an explicit product. [Integer matrices and both first actions](sixth.json.gz)

```sh
python3 algorithms/dual_singlets_20261118/rank_four_source.py --plan /tmp/so9-s5/plan.json.gz --first /tmp/so9-s5/first.json --work /tmp/so9-source --workers 24
python3 algorithms/dual_singlets_20261118/rank_four_sixth.py --source /tmp/so9-source/source.json.gz --primitive /tmp/so9-source/primitive.json.gz --out /tmp/so9-sixth.json.gz --workers 4
```

## Two additional orthogonal classical singlets at derivative bidegree $(3,3)$

At level 30 the twelve-letter terminal space contains six gauginos and six scalars with three flavour degrees $(2,2,2)$. Modulo the classical graded commutators, each trace is represented by its letter multiplicities. A nonzero trace has even length and contains each gaugino component at most once. The weight space contains 257 trace monomials. Quotienting by spin and flavour lowering gives a ten-dimensional formal singlet. Three independent singlets survive on a rank-four Cartan. In the Cartan restriction the term assigning distinct Cartan slots to the traces of a product is its orbit sum with nonzero coefficient. For at most four traces these orbit sums have distinct leading terms.

For $USp(8)$ cut open five traces of any product with more than four traces. Leave one pair of free defining indices at the ends of each open trace. The Pfaffian on these ten free indices vanishes. Pairing the two ends of each original trace reproduces the input with coefficient $+1$ or $-1$. Every other pairing joins at least two of the original traces. Iteration reduces the product to at most four traces. These shorter products form a spin- and flavour-invariant subspace, and their Cartan orbit sums give the matching lower bound at every weight. In derivative bidegrees $(t,t)$ for $t=0,1,2,3,4,5$, the level-30 terminal classical singlet dimensions are therefore $(0,6,0,3,0,0)$. Direct integer contraction verifies every leading coefficient for the 2,878, 150 and one excess-trace words at $t=1,3,5$. For $t=3$ seven projected Pfaffians give determinant $-1024$. Independent A8 sampling gives determinant 4. [Triangular reduction](pfaffian_reduction.json) · [Sampled Pfaffians](mixed_Sp.json) · [A8 job 5671626](a8_5671626.json)

For $SO(9)$ the action on adjoint matrices factors through $O(9)/\{I,-I\}$. The orthogonal second fundamental theorem generates every defining-index relation by antisymmetrizing ten slots against ten slots and pairing the rest. At twelve letters the remaining four slots give two fixed metric edges. Together with the two slots of each adjoint letter, these form ten paths and possibly a closed trace. The nonvanishing length patterns consist of one length-three path or two length-two paths among single letters, or ten single letters with a separate length-two trace. A closed single letter has zero trace.

Within a path the classical relations put its labels in graded order. Repeated identical gauginos vanish. Transposition of a path of length $k$ acts by $(-1)^k$ after graded reordering. Thus an even path with both ends on the same side of the antisymmetrization vanishes. Two identical odd paths of the same endpoint type also vanish, since exchanging them permutes two pairs of slots evenly and exchanges two odd coefficients. Permutations within the left or right slot set change only the sign. Interchanging the sets transposes the determinant. Enumerating the resulting path multisets and endpoint types gives 3,198 templates. Of these, 1,791 have nonzero projection. All $10!$ contractions per template give total rank five over $\mathbb Q$, with minor 262144. This exhausts the finite-colour relations after the classical quotient. [All pivot relations and projection digest](orthogonal_paths.json) · [A8 job 5671813](a8_5671813.json) · [Orthogonal second fundamental theorem](https://arxiv.org/abs/1102.3221)

$$
\dim H_{Q_0}^{\rm terminal}(SO(9))=10-5=5,\qquad
\dim H_{Q_0}^{\rm terminal}(USp(8))=10-7=3 .
$$

Two rational singlet representatives have 68 and 69 trace monomials. They come from symplectic identities, so both vanish on the common Cartan, and the two saved orthogonal functionals evaluate to the identity matrix on them. A second test replaces each scalar by a sum of the two gauginos with independent odd coefficients and extracts the product of all six auxiliary coefficients. The classical supercommutators still vanish. The two classes give $0$ and $-2[S^6]$, respectively, in the pure orthogonal quotient. Incoming quantum actions are required to compare the groups beyond the classical page. [Representatives and both projectors](classes.json.gz) · [Pure-gaugino substitution](scalar_substitution.json)

```sh
cd algorithms/dual_singlets_20261118
python3 rank_four_pfaffian_reduction.py --out /tmp/pfaffians.json
python3 rank_four_relations.py --work /tmp/rank-four-relations --workers 48 --orthogonal 2000 --symplectic 80000 --seed 1627924
python3 rank_four_orthogonal_paths.py --work /tmp/rank-four-paths --workers 48
python3 rank_four_scalar_substitution.py --symplectic rank_four/mixed_Sp.json --orthogonal rank_four/mixed_SO.json --output /tmp/substitution.json
python3 rank_four_mixed_classes.py --sp rank_four/mixed_Sp.json --so /tmp/rank-four-paths/orthogonal_paths.json --functional /tmp/substitution.json --output /tmp/mixed-classes.json.gz
```

## The first differential on the two mixed singlets

Use the rational ten-letter $SO(9)$ source $F$ saved above. On a Taylor letter with positive spin-derivative exponent, replace one such derivative by flavour derivative $i$. If $i$ already occurs the result is zero. Otherwise the sign is $(-1)$ to the number of existing flavour labels larger than $i$. The Taylor factorial cancels the exponent of the chosen spin derivative. Extend this odd operation to a product by the graded Leibniz rule. Each operation anticommutes with $Q_0$, by the linear substitution of the odd coordinates in the classical superfield equation. The verifier also expands all 102 anticommutators on the 17 letters in the source. A linear substitution of the letter coefficients preserves the nine-colour polynomial identities.

In the following table the replacements are performed in the stated order. First use the plus spin derivative for the listed flavours in increasing order, then the minus spin derivative for its listed flavours. The two columns evaluate $Q_1$ on the rational functionals dual to the 68- and 69-term classes. Since the functionals are singlets, these entries also hold after projecting the sources to their spin and flavour singlets.

| Bilinear multiplying $F$ | Plus flavours | Minus flavours | First functional | Second functional |
|---|---|---|---:|---:|
| $S$ | $0,1,2$ | $0,1,2$ | $-151200$ | $-12600$ |
| $\operatorname{Tr}(\lambda_+\phi_0)$ | $0,1,2$ | $1,2$ | $-25200$ | $-2100$ |
| $\operatorname{Tr}(\phi_0^2)$ | $1,2$ | $1,2$ | $-12000$ | $-600$ |

The first and third rows have determinant $-60480000$. A classical boundary has zero projected first action by $Q_1Q_0=-Q_0Q_1$. The two closed sources thus give independent first differentials. The second row is one sixth of the first. The sources are constructed from the previously certified $F$.

The contracted calculation collects the colour reconnections before applying the replacements. The independent replay instead canonicalizes each whole trace word after every replacement, applies the original $SO(9)$ pair operator and reduces the terminal traces afterwards. Each functional vanishes on products of fewer than five traces. A pair action changes the trace count by at most one. A pair of odd-length traces also requires joining to remove their odd lengths, leaving at most four traces for a twelve-letter input. These filters leave 1,273 possible parents for each independent row. The rational results agree entry by entry. [First row](lowered_gauginos.json) · [Mixed bilinear](lowered_mixed.json) · [Third row](lowered_scalars.json) · [Independent calculation](lowered_replay.json) · [A8 job 5671956](a8_5671956.json) · [A8 job 5671971](a8_5671971.json)

```sh
cd algorithms/dual_singlets_20261118
python3 rank_four_lowered_pair.py --channel gauginos --source rank_four/source.json.gz --classes rank_four/classes.json.gz --out /tmp/first-row.json
python3 rank_four_lowered_pair.py --channel scalars --source rank_four/source.json.gz --classes rank_four/classes.json.gz --out /tmp/third-row.json
python3 rank_four_lowered_replay.py --source rank_four/source.json.gz --classes rank_four/classes.json.gz --output /tmp/first-replay.json --workers 8
```

A8 job 5672555 appended three independent Cartan functionals to the calculation. They span the annihilator of the seven symplectic Pfaffian rows in the ten-dimensional formal singlet. Each of the three first actions has Cartan coordinates $(0,0,0)$. The two independent images therefore span the entire two-dimensional kernel of Cartan restriction in the orthogonal classical target. This removes both candidates from the second page. [Five coordinates for the first row](full_gauginos.json) · [Five coordinates for the third row](full_scalars.json) · [Slurm receipt](a8_5672555.json)

At $(t,t)=(1,1)$ there are fourteen terminal letters. Grouping the 6,702 weight words by the length and gaugino count of each trace gives the 18-dimensional formal singlet. Twelve orthogonal determinants have a minor of $17529450901606728990720$. The Cartan lower bound is six. The orthogonal classical dimension is therefore $18-12=6$, equal to the symplectic dimension. In the order $t=0,1,2,3,4,5$, the terminal classical dimensions are $(0,6,0,5,0,1)$ for $SO(9)$ and $(0,6,0,3,0,0)$ for $USp(8)$. All three excess orthogonal dimensions are hit by certified first differentials. [Integer relations](scalar_terminal.json) · [Replay](../rank_four_scalar_terminal.py)
