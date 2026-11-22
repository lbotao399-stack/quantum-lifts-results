# Classical targets for the $G_2$ level-36 singlet and the $USp(4)$ weight $(7,5,4)$

Resolving the entire high-axis pair, including its Taylor determinant and both colour lowering maps, improves the rational second-target bound from nine to **one**. Only the $(8,0)$ colour can contribute. [Determinant projection and certificates](determinant/README.md)

For the $G_2$ flavour singlet $(6,6,6)$, consider a possible 15-letter classical source. Under filtered Taylor equivariance its higher supercharge targets have the following rational bounds. The index places a generator at one or more of lengths 11, 13 and 15. Its precise length and the physical second differential remain open. [Source index and exclusions](../second_odd_generator.md) · [Supercharge selection rule](../../../perturbative_cartan_targets.md)

| Supercharge order | Taylor diagram | Target letters | Cochains | Upper bound on classical cohomology |
|---:|---|---:|---:|---:|
| 2 | $(5,5,3,3)$ | 14 | 38,435 | 9 |
| 3 | $(6,6,3)$ | 13 | 7,768 | 0 |
| 4 | $(7,7)$ | 12 | 229 | 0 |

Orders five and higher have no polynomial target for this source by the same selection rule. This leaves the first two orders for a 15-letter source. [Independent character counts](../5691072/) · [Rational ranks and characters](result.json)

## The three classical subcomplexes

Polarize the Taylor directions on ordinary even axes and take the dual Chevalley–Eilenberg homology of the Lie algebra of positive polynomial currents. The boundary merges two current letters by their Lie bracket and adds their Taylor indices. At each finite degree the rational multiplicities equal those in the dual cohomology. The multilinear Taylor module recovers the physical super Taylor multiplicities by signed polarization in characteristic zero. [Taylor functor and physical projection](../../taylor_hooks/README.md)

Separate two individual high-degree axes and one transverse space. Give a current supported in one of these three blocks cost zero and a current meeting different blocks cost equal to the square of its total Taylor degree. Merging letters preserves or increases cost. The cost-preserving boundary computes the three pure homologies independently. Each mixed Taylor monomial carries an exterior algebra on the adjoint colour representation. Taking colour invariants and projecting onto a permutation-group character are exact over the rationals. The resulting finite filtration bounds classical homology in each projected component. Permuting pure factors of homological lengths $y,t$ contributes the sign $(-1)^{yt}$.

For a diagram $(h,h,\ldots)$ with tail $\nu$, removing $(h,h)$ leaves the straight diagram $\nu$. Restriction to the two high axes and the transverse space therefore contains $\det^h\otimes\mathbb S_\nu$ once. On degrees $(h,h)$ its high-axis exchange has eigenvalue $(-1)^h$. Polarizing the transverse slots and pairing with the character of $\nu$ selects their Schur multiplicity. Other constituents can contribute to the upper bound.

In $G_2$ the nonzero single-axis homologies through energy seven are listed below. Each entry gives letter length and total colour dimension. The replay constructs the 14 colour basis matrices in the seven-dimensional representation, checks every bracket, builds the integer chain matrices by weight and checks consecutive products before computing rational ranks. Energy four has zero homology at every length. [Colour matrices](g2_colours.py) · [Chain calculation](transverse.py)

| Single-axis energy | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| Length and colour dimension | $0,1$ | $1,14$ | $2,77$ | $3,182$ | empty | $4,924$ | $5,1547$ | $5,1254$ and $6,748$ |

## The second target in $G_2$

At weight $(5,5,3,3)$ and length 14 there are exactly two occupations on the separate-axis page. Either the four pure lengths are $(4,4,3,3)$, or they are $(3,3,3,3)$ with two copies of the same quadratic mixed jet on the two high axes. To check exhaustiveness, choose pure lengths from the preceding table. Each remaining mixed letter consumes at least two slots and at least one slot outside each individual axis. Both surviving occupations use all transverse slots in pure letters. After the required axis exchanges their colour-singlet upper bounds are 79 and 1,155. Refining the transverse block to separate axes proves that terms meeting both a high and a transverse axis give zero at this length on the three-block page.

Now keep the transverse axes together. At transverse degree $(3,3)$ and length six, each factor lies in its 182-dimensional space of three-letter cycles. The exchange includes the sign from their odd lengths. The required negative exchange selects their ordinary symmetric square. The remaining boundary merges one first jet from each axis. There is no chain of length seven. Its exact kernel has dimension 5,720. The neighbouring $(4,2)$ length-six homology is zero by the single-axis test. The $(3,3)$ component therefore has precisely the required transverse Taylor highest weight. [Cross-bracket matrix](g2_terminal.py)

| $G_2$ colour highest weight, short-root label first | Multiplicity in the 5,720-dimensional kernel | High pair with zero mixed jets | High pair with two mixed jets |
|---|---:|---:|---:|
| $(0,0)$ | 1 | 0 | 0 |
| $(4,2)$ | 1 | 5 | 41 |
| $(6,0)$ | 1 | 2 | 39 |
| $(8,0)$ | 1 | 1 | 12 |

The last two columns come from the alternating square of the 924-dimensional module and from the symmetric square of the 182-dimensional module tensored with the alternating square of the adjoint, respectively. All $G_2$ modules are self-dual. This gives $8+92=100$.

In the high pair with two mixed jets, the cost is eight. Merging pure first jets from different axes raises cost by four. The three highest-weight kernel dimensions, in the order of the nontrivial colour rows, are $16,8,7$. The image lies in top-length pure-axis cycles, so its chain rank gives its rank on the auxiliary page. On its kernel, merging a first jet with a quadratic mixed jet raises cost by five. The two orientations are exchanged by swapping the axes. The cost-thirteen target again has pure axes at their maximal lengths. An earlier boundary would require cost nine with one cubic jet and residual pure degrees $(3,4)$ or $(4,3)$. The energy-four homology is zero. The joint kernel dimensions are $0,0,1$. Bracketing the two mixed jets into a quartic jet leaves these kernel dimensions unchanged. These ranks concern the classical boundary. [Exact brackets and weight ranks](mixed_bracket.py)

The successive upper bounds on the physical second target are

$$
38{,}435\ \longrightarrow\ 1{,}234\ \longrightarrow\ 100\
\longrightarrow\ 8+16+8+7=39\ \longrightarrow\ 8+1=9 .
$$

The upper bounds by colour are $5,2,2$. Further classical relations, earlier physical $d_1$ images, and the pairing with a source can reduce them.

## The third and fourth targets in $G_2$

The $(6,6,3)$ diagram contains the weight $(5,5,5)$ once by the tableau rule. At that weight and length 13, three pure energy-five axes supply only 12 letters. With one mixed quadratic jet, both incident pure axes have energy four. With two quadratic jets every pure axis must have energy at most three to meet the length, requiring at least six slots in mixed jets, which use only four. A single cubic mixed jet uses three slots and also cannot reduce all three energies to at most three. This exhausts the two units of excess Taylor degree. Thus the classical target is zero.

At $(7,7)$ and length 12 only the two pure energy-seven, length-six factors contribute. The traces of the identity and of the exchange on their colour singlets are both one. The required negative eigenspace has dimension $(1-1)/2=0$. [Occupations, tableaux and Weyl averages](exceptional.py)

## The second target in $USp(4)$

At BMN flavour $(7,5,4)$ the second target has diagram $(5,5,2,1,1)$. Its last four slots transform as the permutation representation of Young shape $(2,1,1)$. The table gives its multiplicity after the high-axis sign projection, first keeping only the two single-axis differentials, then retaining the full transverse differential. At target length 12 the signed permutation-class traces give $(2838+1740-102+228)/48=98$. A four-slot control independently computes the full rational cohomology and compares it with the filtration in every length. [Permutation and colour projection](paired_axes.py)

| Source letters | Target letters | Cochains | Separate-axis upper bound | Permutation projection | With transverse homology |
|---:|---:|---:|---:|---:|---:|
| 11 | 10 | 49,802 | 349,460 | 50,376 | 36,576 |
| 13 | 12 | 8,771 | 1,128 | 225 | 98 |
| 15 | 14 | 99 | 0 | 0 | 0 |

The unprojected bound can exceed the cochain dimension of the selected representation because the weight contains other representations. With Python 3.12 or later and python-flint, replay from this directory.

~~~sh
python3 verify.py | diff -u result.json -
~~~
