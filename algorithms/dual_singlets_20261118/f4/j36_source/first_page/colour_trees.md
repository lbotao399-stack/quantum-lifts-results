# Colour contractions in the second differential of the sixth glueball

Consider $F_4$ at $j=36$ in the full relative holomorphic jet complex. The two first dotted gauginos are $U,V$. In the source with one matter fermion there are four copies of each dotted gaugino, four scalars and one matter fermion. The terminal component of its classical differential has the same gauginos and six scalars. On a classically closed source the corrected second obstruction is the negative of the intermediate covector $\eta$ applied to the one-loop supercharge $Q_1$. The joint test keeps the classical image modulo the inactive splittings and the scalar $\eta Q_1$. Ghost slides past the matter fermion preserve both entries by the two adjoint ghost actions on the exchanged colour edge. [Joint test and signs](trace_obstruction.md)

## Eliminating adjoint contractions on a tree

Form a multigraph with one vertex for each defining trace and one edge for each contraction of two Lie algebra slots by the inverse of $B=\operatorname{Tr}_{26}/6$. Physical letters decorate the vertices. If this multigraph is a forest, every source reduces to the previously tested span of 19 defining-trace singlets for the joint test. Every terminal singlet in the forest reduces to ordinary defining traces as well. This is a spanning statement. Additional native identities in the trace span can still change the obstruction.

Cut a leaf at its unique contracted slot. All its physical letters now form a single ordered word. Sliding $U,V$ through the other letters leaves at most one of each on this leaf. Contracting the edge replaces the word by six times its trace-orthogonal projection to the adjoint. A leaf with no physical letters is zero. For one Lie input the projection is the identity. For two or four symmetrized colour inputs it is zero by metric transposition against the skew test generator. For three symmetrized inputs $A,B,C$ the projection is

$$
\frac{A\operatorname{Tr}_{26}(BC)+B\operatorname{Tr}_{26}(AC)+C\operatorname{Tr}_{26}(AB)}{36}.
$$

This follows by differentiating $12\operatorname{Tr}_{26}(X^4)=(\operatorname{Tr}_{26}X^2)^2$ and polarizing the remaining three slots. These are identities of multilinear colour tensors. Coefficients in the exterior algebra acquire the corresponding Koszul signs. An ordered word differs from its symmetrization by adjacent exchanges. Each exchange replaces two Lie inputs by their Lie bracket and reduces the number of arguments by one. A leaf with at most four physical inputs therefore contracts to Lie words carried by its neighbour and ordinary trace spectators. Expanding each bracket back into its two matrix words leaves one fewer contracted edge. The verifier checks the resulting projection against all 52 Lie generators on every ordering through degree four for two integer assignments. [Exact recursion](adjoint_trace_trees.py) · [Projection tests](adjoint_trace_trees.json)

At the source, a leaf that resists this operation would need at least three of the five nonghost letters. Each nontrivial tree has at least two leaves. Thus there is always a reducible leaf. Induction on the number of contracted edges gives the source assertion, including disconnected forests.

At the terminal weight, two resistant leaves must each contain three scalars and both dotted species. Together they use all six scalars. At either leaf exchange the three scalar arguments and average over their permutations. Every error term has a scalar bracket and at most four Lie inputs on that leaf, so contracts by the preceding formula. The average has three symmetric scalar colour slots. The flavour singlet has slot diagram $(2,2,2)$. Its character on the identity, a transposition and a three-cycle is $5,-1,-1$. Its multiplicity invariant under permuting this triple is

$$
\frac{5+3(-1)+2(-1)}{6}=0 .
$$

This eliminates every remaining terminal tree by induction. The distinction at degree five can be seen on the Cartan input $(2,3,5,7)$. Pairing the input and its fifth matrix power with Cartan generators one and three gives the determinant $12(-8718)-(-18)(10932)=92160$. The primitive quintic covariant really is present before the flavour projection. [Native Cartan pairings](adjoint_trace_trees.json)

## A contraction through the 26-dimensional Jordan representation

Remove two spectator factors $\operatorname{Tr}_{26}(UV)$ and keep all ten remaining terminal letters. The full colour character in this dotted and flavour singlet is 162. For a genuine exceptional contraction use the covariant obtained by projecting a matrix word onto traceless Jordan multiplication operators in the Albert algebra. This projection vanishes on words of length at most three. On four Lie inputs it is symmetric as a colour tensor. On five symmetrized inputs it vanishes by transposition because the multiplication operators are self-adjoint. [Quadratic component identity and quartic symmetry](../mixed_colours.md)

Take a nine-input polynomial consisting of two words contracted through this single Jordan projection, with any ordinary defining-trace spectators. Each word needs at least four inputs, and a one-letter spectator trace vanishes. The two lengths must therefore be four and five. Commuting adjacent arguments reduces the five-input word to quartics with one bracket. Performing a classical split on one second jet inserts one further bracket. It either nests inside the existing bracket, shares its quartic, or belongs to the other quartic. At least one innermost bracket contains a dotted gaugino when the split is inactive. Conversely each such bracket has a second-jet primitive. Its $Q_1$ vanishes even after restoring the two spectator factors, since every individual jet in the primitive has at most one flavour label.

Contract the four dotted labels in pairs with the invariant alternating tensor and the six scalar labels in triples with two flavour determinants. These contractions span both singlet spaces by the first fundamental theorem for the defining representations of $SL(2)$ and $SL(3)$. Enumerating the layouts modulo the graded quartic symmetries gives seven independent formal inactive rows. The nested, same-quartic and opposite-quartic cases supply one, two and four, respectively. Five additional rows with unrestricted brackets or a quadratic spectator serve as controls. The 12 rows use 354 ordered graphs. [Singlet enumeration](jordan_bracket_singlets.py)

Evaluate all 370 stable defining traces and the new colour polynomials on the same exterior-algebra assignments in the Albert representation. Its invariant unit line contributes zero to every trace with an adjoint input. The denominator of the inverse Jordan trace Gram matrix is 72. On the joint evaluation the ranks of the trace rows, of the trace rows with seven inactive rows, and of all rows are respectively $153,157,160$. Every rational relation amongst the first $370+7$ evaluated rows has zero image after restoring the two glueball factors and imposing the ordinary inactive and quartic relations. Hence every genuine relation of this type vanishes in the five-coordinate test, including any that mix the two kinds of colour tensors. This implication only uses containment of the true relations in the kernel of evaluation. [Rational elimination](jordan_bracket_test.py) · [Pivot rows and determinants](jordan_bracket_test.json)

## The missing nine-input colour tensor

The nine-letter inactive primitive has a single second jet. For a jet with two dotted labels the two other dotted slots have symmetric spin and alternating colour. The scalar colour diagram is $(2,2,2)$. For a differentiated scalar, the three other dotted colours have diagram $(2,1)$ and the five other scalar colours have diagram $(2,2,1)$. The second jet supplies one further adjoint in both cases. The Weyl character and direct native evaluations give

| Second jet | Native singlets | Rank attained by defining traces | Rank attained with one Jordan contraction |
|---|---:|---:|---:|
| Two dotted labels | 85 | 84 | 84 |
| One dotted and one scalar label | 299 | 294 | 298 |

The classical primitive in the physical jet convention is fixed by expanding its graded bracket with the native $Q_0$ and accounting for the odd differential passing the preceding letters. Seven Jordan primitives occur with a differentiated scalar. A8 job 5704736 used a fresh integer frame and reproduced both ranks and the terminal ranks $153,157,160$. [Source characters and evaluator](jordan_inactive_sources.py) · [Local pivots and A8 receipts](jordan_source_certificates.json)

On three ordinary adjoint variables, project the nine colour slots to the Young diagram $(3,3,3)$ by alternating down each of three columns. The invariant character is 11. All 18 formal defining-trace singlets give an integer evaluation matrix of rational rank 10. A single Jordan contraction contributes zero to this diagram. Its two quartics and one bracket leave four individual slots symmetric on one of the quartics. A diagram of first-row length three has zero projection on four symmetric slots.

To supply the eleventh tensor, leave a free adjoint input on each quartic. Join the three physical arguments of the first quartic to those of the second by an adjoint edge carrying three further inputs. In order along the edge these act by three successive Lie brackets. Use the existing inverse Jordan trace Gram matrix between the quartics and the existing inverse $B$ on the adjoint edge. Alternate each column containing one input at each endpoint and one on the edge. This is a polynomial of degree three in each of the three colours, transforming by the third power of the determinant. Its 216 signed terms reduce to 156 ordered terms. On exactly the same 16 integer assignments, adding it gives a minor of order eleven with residues $180$ modulo $1009$ and $4$ modulo $1013$. The evaluation is therefore injective over the rationals on the whole 11-dimensional native multiplicity space. [Construction and exact trace rank](nine_colour_tensor.py) · [Both minors](jordan_source_certificates.json)

The restriction of the slot representation $(3,3,3)$ to the other eight slots contains $(3,3,2)$ once. Restricting further gives $(1,1)\otimes(2,2,2)$ once for the two-dot jet and $(2,1)\otimes(2,2,1)$ once for the differentiated scalar. The rational permutation representation is semisimple, and the span of defining traces and single Jordan contractions is stable under every slot permutation. The missing copy consequently contributes at least one missing primitive in each physical sector. The rank bounds in the table leave room for exactly one in each. This identifies their colour type and supplies the tensor from which both primitives can be projected.

## Replay

Run these commands from this directory with NumPy, SciPy and python-flint.

~~~sh
python3 adjoint_trace_trees.py
python3 jordan_bracket_test.py --out /tmp/f4-jordan.json
python3 jordan_inactive_sources.py --kind 2 --out /tmp/f4-second-ghost.json
python3 jordan_inactive_sources.py --kind 1 --out /tmp/f4-differentiated-scalar.json
python3 nine_colour_tensor.py --prime 1009 --out /tmp/f4-nine-1009.json
python3 nine_colour_tensor.py --prime 1013 --out /tmp/f4-nine-1013.json
~~~
