# The sixth glueball power of F4 survives the first supercharge differential

In the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, the glueball power $S^6$ at $j=36$ is nonzero on $E_2$, the cohomology of $d_1$ on the classical page. Dot's earlier $E_2$ claim survives independent recalculation. The proposed $d_2$ source has a classical primitive. The incoming $d_2$ from other sources remains open. [Classical primitive of the proposed source](../README.md) · [Dot's withdrawal](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1#issuecomment-7074126654)

## The classical target

Use Dot's $B=\operatorname{Tr}_{26}/6$, with the first dotted gauginos ordered before the second. The cohomological target is the pure twelve-gaugino sector. After imposing the separate same-dot brackets, the two six-ghost factors have two common maximal-Casimir representations. Enumerating the six-root abelian ideals of a Borel and applying the Weyl dimension formula gives dimensions $7113106,1850212$. This uses Kostant's decomposition of the exterior algebra modulo the bracket ideal. The exceptional chiral-ring theorem also identifies the classical invariant quotient as $\mathbb C[S]/(S^9)$. [Chiral-ring theorem](https://doi.org/10.1002/prop.202400027)

For the mixed classical bracket, keep a highest vector of the first five-ghost module and a lowest vector of the second. The remaining colour lies in the adjoint representation. Testing the full allowed weight gives eight rows. One row is $(-1672/273,1472/231)$. The covector with branch weights $(2392,2299)$ annihilates all eight. Denoting it by $\lambda$, its value on $S^6$ is

$$
\lambda(S^6)=-6!\big(2392\cdot7113106+2299\cdot1850212\big)=-15313094596800 .
$$

An independent integer covector, evaluated on 55 exterior probes in the embedded $\mathfrak{so}_9$, kills all 146 classical predecessor rows and gives $720$ on the same target. On this subalgebra $B=\operatorname{Tr}_9/2$, checked on all 666 unordered pairs of tangent generators. Restriction commutes with the classical bracket. The orthogonal invariant tensor theorem spans the restricted predecessors by defining traces, since these are adjoint inputs in odd dimension. These two covectors use different scales. [Integer weights](../../../subgroup_normal_jets/f4_sixth_dual.json) · [Restriction and integer replay](../../../subgroup_normal_jets/f4_sixth.py)

## The complete first-page test

On an incoming classical cycle $X$, we need $\lambda(Q_1X)=0$. Dot's 40 rational coefficients determine an invariant functional $\eta$ satisfying $\lambda Q_1=\eta Q_0$. We used the original Lie brackets, the 40 extremal seed tensors and their coefficients. The local evaluator computes the exterior actions, determinant pairings and Casimir recursions from those inputs. It reconstructs the target's branch weights directly from the mixed classical bracket.

| Incoming profile | Exact test |
|---|---|
| Matter fermion and complementary scalar | On all $52^2=2704$ ordered colour pairs, the scalar lift composed with the classical splitting equals the symmetrized one-loop colour tensor. All 40 homogeneous tensors annihilate the split. |
| Differentiated scalar | All $326+437=763$ residuals vanish. The inhomogeneous term is nonzero in $311+409=720$ entries and cancels the homogeneous term in each. |
| Second ghost jet | Same-dot splittings vanish in the exterior quotient. For a mixed splitting all 6,920 raw zero-weight components of the homogeneous term vanish. The particular term inserts extra ghosts into the classical bracket ideal, which $\lambda$ annihilates. |

For the scalar identity, let $k$ be classical splitting on the symmetric scalar pair and $C$ its quadratic Casimir. On all 1,378 symmetric basis vectors,
$$
2k^t k=C-54,\qquad C(C-26)(C-40)=0.
$$
This inverts $k^t k$ on its complete domain. Each homogeneous seed lies in the kernel of the ordered bracket contraction. We checked the raising and lowering powers on every seed and checked the exchange phase on the two diagonal blocks.

The extremal tests span each invariant residual. An irreducible highest- or lowest-weight module is generated from the indicated vector. In an invariant trilinear form, move lowering operators off the highest-vector leg, then move raising operators off the lowest-vector leg. The latter kill the fixed highest vector. Values on the resulting entire third weight space therefore determine the original functional. With a differentiated scalar the old four-ghost bracket ideal remains an ideal after insertion of the next ghost. With a mixed second ghost this holds on both four-ghost legs. This reduces the tests to the abelian-ideal components in exterior degrees four and five. Full raw weight enumeration in the last two rows eliminates dependence on the choice of reduced tests.

The bare positive-order supercharges vanish on the zero-flavour target. Together with $\lambda Q_1=\eta Q_0$ and $\lambda(S^6)\ne0$, this proves $[S^6]\ne0$ on $E_2$. A8 job 5701911 regenerated Dot's two smaller spanning frames and all 145 and 185 rational values with their original hashes. A8 job 5702028 ran the separately written implementation over all the raw weights above. [Both Slurm results and hashes](replays.json) · [Colour algebra](colour.py) · [Full dual evaluation](ward.py)

## Testing the second differential

For one second jet and otherwise first jets, the bare $Q_2$ requires six flavour slots on three letters. The letters have at most four. If $Q_0X=0$ and $Q_0X_1=-Q_1X$, the target coefficient is determined by
$$
\lambda(Q_2X+Q_1X_1)=-\eta(Q_1X).
$$
This evaluates candidate sources before computing $X_1$. A nonzero test still needs a source on $E_2$.

Fix the eight ghost spectators, two first-flavour scalar colours, and a matter fermion carrying flavours two and three. Removing these spectators leaves one fermion colour and one scalar of each of flavours two and three. The classical split maps the 140,608-dimensional raw colour space into the tensor product of two symmetric squares. Exact elimination over the rationals on 17 Weyl orbits gives
$$
\dim\ker Q_0=1430={52\choose2}+2\cdot52.
$$
The first 1,326 vectors are differentials of two matter fermions. The remaining 104 multiply one of the two invariant bilinears $B(\psi_{23},\phi^2),B(\psi_{23},\phi^3)$ by a free scalar. Every listed vector is checked in the computed kernel, and their joint span has full dimension. Arbitrary additional first-jet spectators commute with the classical differential. Gauge averaging a bilinear times such spectators preserves the invariant bilinear factor. Its bare supercharge vanishes by flavour support. The spectator invariant has only first jets, so its bare positive-order supercharges also vanish by flavour counting. Under the marked product theorem these classes admit a perturbatively closed product. This excludes this entire separately closed single-flavour family as a source hitting $[S^6]$ at $d_2$. [Rational kernels](one_fermion.json) · [Independent basis calculation](one_fermion.py)

As a direct residual test, taking both first-flavour colours equal to the first Cartan basis element gives 1,212 nonzero scores on 1,456 raw columns. All 38 independent classical cycles have score zero. Four further colour choices also have zero score on their entire classical kernel. [Inputs and scores](second_samples.json) · [Direct evaluation](second_test.py)

Allow all three matter-fermion flavours and factor common scalar spectators out of the source. At residual scalar degrees $(1,1,1)$, the classical kernel has dimension $8112-8058=54$. At $(2,1,1)$ it has dimension $352872-347257=5615$. The corresponding spaces are exhausted by explicit classical boundaries and scalar multiples of the traceless bilinears. For the two diagonal traceless bilinears the colour coefficients of $Q_1$ are $36(1,-1,1)$, proportional to the classical coefficients. The two classically closed combinations are also quantum closed. [Three bilinear columns](../../../subgroup_normal_jets/f4_first_columns.py)

Each proposed cycle was differentiated over the rationals. Modular lower bounds on the rank of the classical map and on the rank of those cycles sum to the exact domain dimension. Flavour permutations give the same result at permuted degrees. Cancellations involving differentiated scalar or ghost sources, and mixed-flavour cycles of higher residual scalar degree, require further tests. [Joint calculation](mixed_syzygies.py) · [Three scalar slots](mixed_111.json) · [Four scalar slots](mixed_211.json)

## Five active Taylor slots of arbitrary parity

At $j=36$ an incoming second-order source for $S^6$ has dotted degrees $(4,4)$, flavour degrees $(2,2,2)$ and thirteen letters. Its classical differential splits its one second jet. Call a slot active if it belongs to the factor in which the classical cancellation takes place. The remaining first jets can have arbitrary colours. Every source whose closure factors over at most five active slots is a sum of classical boundaries and products with a closed bilinear. This includes cancellations between matter fermions, differentiated scalars and second ghost jets.

| Active slots after the classical split | Source columns | Classical image rank | Boundary rank | Additional cycles from the bilinear |
|---|---:|---:|---:|---:|
| Scalar degrees $(2,2,1)$ | 14,763,840 | 14,262,508 | 487,968 | 13,364 |
| Four same-dot ghosts and one scalar | 4,734,704 | 4,325,125 | 408,305 | 1,274 |
| Three same-dot ghosts and one scalar | 209,560 | 201,449 | 8,059 | 52 |
| Two same-dot ghosts and one scalar | 5,408 | 5,355 | 52 | 1 |

Each row meets the rank bound separately at every dominant colour weight. The integer kernel checks and the complementary modular ranks prove saturation over the rationals. The first two rows used 32 Weyl orbits each. In the dotted slice the additional cycles multiply $2B(\phi,c_2)+B(D\phi,c_1)$ by undifferentiated first ghosts, where $c_1,c_2$ are the divided first and second jets in that dotted direction. Both terms of its classical variation cancel. A8 job 5702371 completed both five-slot rows in 101 seconds on 16 allocated CPUs. The sparse solver pivots on the shortest available row or column. [A8 hashes and accounting](polarization_a8.json) · [Scalar case](mixed_221.json) · [Dotted case](dotted_41.json) · [Dotted differential](dotted_syzygies.py) · [Native solver](sparse_rank.cpp)

To see why the two slices suffice, polarize $n$ distinct odd Taylor coordinates and quotient the one-second-jet cycles by classical boundaries and all degree-three bilinear cycles times first jets. Tree splitting is natural in the Taylor variables, so the permutation group of the $n$ labels acts on this quotient. Equal scalar labels select the trivial character of their permutation subgroup. Equal dotted labels select the sign character there. This last sign comes from changing the parity of the repeated coordinate. Equivalently, polarize an even coordinate as a sum of odd coordinates times distinct auxiliary odd parameters and take their top coefficient. Over the rationals the averaging projections are exact and commute with taking the quotient.

At $n=5$ the scalar $(2,2,1)$ slice detects the partitions $(5),(4,1),(3,2),(3,1,1),(2,2,1)$. Inducing the sign of the four ghost slots with one scalar detects $(2,1,1,1)$ and $(1,1,1,1,1)$. At $n=4$ the $(2,1,1)$ slice and the three-ghost slice likewise cover every partition. At $n=3$ the all-distinct scalar slice is the full multilinear complex. At $n=2$ the colour cobracket is injective because the Lie bracket spans the simple algebra. Thus the quotient vanishes across the stated range. The signed bilinear is the degree-three mixed component. By exactness of polarization it lifts from the degree-three multilinear bilinear space. [Character calculation and rank checks](polarization.py)

After the spectator colours are restored, gauge averaging moves onto the spectator coefficient of each invariant bilinear. The bilinear has fewer than three flavour slots unless all three are scalar flavours, in which case the already checked traceless combination has zero first image. An invariant formed entirely from first jets has too few flavour slots on the active inputs of every bare positive-order supercharge. The marked product thus supplies a perturbatively closed representative of each surviving product. A source outside this exclusion requires a cancellation across at least six active slots. The full source has fourteen slots after splitting.

~~~sh
c++ -O3 -std=c++17 sparse_rank.cpp -o /tmp/f4-colour-rank
export N4_SPARSE_RANK=/tmp/f4-colour-rank
python3 mixed_syzygies.py 221 | diff -u mixed_221.json -
python3 dotted_syzygies.py 5 | diff -u dotted_41.json -
python3 polarization.py
~~~

## Replay

Python 3.10 or later suffices for the native Ward calculation. Installing python-flint accelerates the rational arithmetic. The input file contains only native colour data, extremal seeds and coefficients. The optional extraction command recreates it from the original zip.

~~~sh
python3 extract.py /path/to/F4-j36-d2-verification.zip /tmp/f4-inputs.json.gz
python3 ward.py target
python3 ward.py seeds
python3 ward.py psi
python3 ward.py dphi
python3 ward.py mixed
python3 one_fermion.py | diff -u one_fermion.json -
python3 mixed_syzygies.py 111 | diff -u mixed_111.json -
python3 mixed_syzygies.py 211 | diff -u mixed_211.json -
python3 second_test.py 0 0 0 1 7 13 19 29 19 12
~~~

## A single Taylor representation can hit the sixth glueball

At length thirteen the dotted-spin and flavour singlet contains 70,238 invariant cochains. Only the 3,720 cochains of Taylor partition $(4,4,3,3)$ can contribute to an incoming $d_2$ on $S^6$. This restriction holds for the native fan supercharge in the relative jet complex. The classical cohomology and corrected pairing in this last component remain to be determined. [Integer characters and native component checks](source_symmetry.json) · [Replay](source_symmetry.py)

The forward odd Taylor shears commute in the graded sense with the fan supercharge. For each reverse shear $\bar P_0$ take the first-order correction $\bar P_1$ in the [relative Ward calculation](../../../../su3_j30_classical_ranks/bmn_return/README.md#the-nine-letter-target-survives-every-perturbative-page). Both terms preserve the excess of total Taylor degree over letter number. The two-colour identity used there extends to an arbitrary simple metric Lie algebra. For a single input the contraction of two structure constants is its adjoint Casimir times the identity. For two inputs the connected invariant four-leg contractions with two structure constants satisfy Jacobi and have dimension at most two. In the three-dimensional cross-product algebra they are independent: with external colours $(0,0,1,1)$ and $(0,1,0,1)$ the first two pairings $([a,b],[c,d])$ and $([a,c],[b,d])$ take values $(0,1)$ and $(1,0)$. These evaluations also detect every subspace obtained by graded input or output symmetry. The colour permutations used for the coefficient solve span all coordinate and pair tests in this algebra. Thus the exact solve establishes the identity modulo simultaneous gauge rotations in any colour algebra here. The replay independently substitutes the six shear choices in native $F_4$ root and Cartan coordinates, passing 30,715 component equations each through total pair degree four.

All letters of a source with fourteen slots and thirteen letters are first jets except one second jet. Its first correction has thirteen slots and twelve letters, hence also has excess one. On three letters of excess at most one there are at most four Taylor slots. The bare $Q_2$ requires six flavour slots. The anticommutator $\{Q_1,\bar P_1\}$ has differential order at most three. Its terms on disjoint inputs cancel by odd parity. On any surviving set of inputs the two operations together consume five flavour slots. Consequently both $\{Q_2,\bar P_0\}$ and $\{Q_1,\bar P_1\}$ vanish at this excess. Applying $\bar P_0+u\bar P_1$, where $u$ is the perturbative parameter, to the first two descent equations therefore gives $d_2\bar P_0=-\bar P_0d_2$ on these sources. The first-order identity also applies to the excess-two incoming first-page cochains, whose pairs use at most four slots.

Order the two even dotted coordinates before the three odd flavour coordinates. For a Taylor partition $\lambda$, its covariant highest weight is $(\lambda_1,\lambda_2\mid\max(\lambda'_1-2,0),\max(\lambda'_2-2,0),\max(\lambda'_3-2,0))$, where $\lambda'$ is the transpose partition. The classical complex is a subquotient of tensor powers of the Taylor superspace. Its covariant tensor powers are completely reducible in characteristic zero by Schur–Weyl duality. The first Ward identity for the deformed shear gives the induced Taylor action on the second page. The two shear identities for $d_2$, and its even weight shift $(2,2\mid-2,-2,-2)$, make its image transform with the square of the Berezinian character. The target $S^6$ has type $(6,6)$ and highest weight $(6,6\mid0,0,0)$. Subtracting the shift gives $(4,4\mid2,2,2)$ and hence the unique diagram $(4,4,3,3)$.

Independently, the even subgroup permits ten diagrams, each with multiplicity one, from the product of the dotted singlet $(4,4)$ and the transposed flavour singlet $(3,3)$. The integer colour character uses Adams operations on the $F_4$ adjoint and Weyl reflection. In the selected diagram the length-$12,13,14$ cochain dimensions are $10007,3720,704$. In the target diagram the length-$10,11,12$ dimensions are $129,56,25$. These are cochain dimensions before taking the classical differential.

~~~sh
python3 source_symmetry.py | diff -u source_symmetry.json -
~~~

## The six-scalar singlet and the dotted colour filtration

Keep the eight first-ghost spectators and require the six-scalar factor to close separately. At flavour $(2,2,2)$ take the $SL(3)$ singlet. The one-fermion source consists of a fermion colour and a scalar tableau of shape $(2,1,1)$. The dimension before imposing gauge invariance or adding the spectators is

$$
52\left(52\cdot53\cdot51\cdot50/8\right)=45{,}680{,}700.
$$

At every nonzero dominant colour weight the classical kernel is spanned by checked boundaries and products with the closed bilinear. A8 job 5702426 returned all 54 of those weight files. With Weyl orbit multiplicities they cover 45,503,928 columns. The classical rank lower bound is 40,848,456 and the verified cycle-span lower bound is 4,655,472. Their sum is the full dimension there. The job reached its wall-time limit during the remaining 176,772-column zero-weight block. [Rank and closure calculation](six_scalar.py) · [Completed weight files, hashes and final accounting](six_nonzero_weights.json)

The quotient by the known cycles is an $F_4$ representation over the rationals. If only weight zero occurs, all root operators act as zero and the representation is trivial. Since invariant projection is exact, it remains to test invariant source cochains. Their number is three. The invariant cochain dimensions at lengths four, five and six are $(0,3,4)$ by the Weyl character in Taylor shape $(3,3)$. Write the fermion with complementary scalar index and take the invariant flavour contraction with three alternating scalars. For three independent choices of colour tensor, the values of the classical split on three explicit scalar triples are

$$
\det\begin{pmatrix}
72&-12744&5064\\
30258&-34260&6776\\
0&10968&-3608
\end{pmatrix}=292863838464 .
$$

Thus the invariant source has zero classical kernel. This proves saturation of the full rational colour module without the missing numerical zero-weight rank. Restoring the spectator colours and averaging over $F_4$ leaves boundaries and marked products as in the five-slot argument. Hence this entire separately closed scalar-singlet family has zero incoming $d_2$ on the sixth glueball. [Three rational colour tensors, native test points and matrix](six_invariants.json) · [Independent invariant replay](six_invariants.py)

For a second reduction, take one scalar and $n-1$ equal dotted slots. Filter the tree differential by the number of dotted derivatives on the scalar. At excess one this is zero or one. The associated graded differentiates only the pure ghost factors. At scalar order zero the ghost cochain is the adjoint times exterior degree $n-3$ of the adjoint, of dotted degree $n-1$. At scalar order one it is exterior degree $n-2$ and dotted degree $n-2$. In both cases tensor the pure-ghost cohomology with the scalar adjoint. Garland's positive-current Laplacian retains exactly the finite-colour components whose quadratic Casimir equals 18 times the dotted degree in the long-root-square-two convention. Newton's exterior-power identities and Weyl reflection then bound the colour types on every later filtration page. [Garland's formula, Theorem 3.7](https://arxiv.org/pdf/math/0309232#page=22)

Every surviving irreducible contains its indicated test weight. The certificate gives explicit simple-root lowering paths. At each step the current coroot label is positive, so the next weight occurs. Testing full weight spaces thus avoids constructing highest-weight bases. The integer closure checks and modular rank lower bounds exhaust every tested column.

| Active slots | Columns over all colours | Tested columns | Classical rank | Rank of boundaries and products |
|---:|---:|---:|---:|---:|
| 6 | 73,836,100 | 13,461 | 11,999 | 1,462 |
| 7 | 867,186,320 | 39,842 | 33,775 | 6,067 |
| 7 | same sector | 39,813 | 33,873 | 5,940 |

The six-slot test detects 17 possible colour highest weights and the two seven-slot tests detect 23. Both separately closed dotted families have zero quotient by boundaries and products. For this glueball the selected scalar-polarized diagram is $(4,4,4,2)$, of height four. It already excludes any alternation over five dotted labels; the two dotted calculations also apply to other target diagrams. [Proof and replay](dotted_filtration.py) · [All lowering paths and local ranks](dotted_filtration.json)

## Products inside the remaining diagram

Polarize the selected classical diagram with four ordinary scalar variables. Its highest weight is $(4,4,4,2)$ in this convention. Every classically closed linear combination of products in which each factor is itself a gauge-invariant highest-weight cochain is a classical boundary. Products projected from other weight components require a separate test.

For two highest weights adding to $(4,4,4,2)$, each factor has shape $(a,a,a,b)$, with $0\le b\le a$. The possible proper complementary pairs, recorded by $(a,b)$, are $(1,0)+(3,2)$, $(1,1)+(3,1)$, $(2,0)+(2,2)$, and $(2,1)+(2,1)$. At excess one only one factor can contain a second jet. For $(2,1)$ the terminal invariant cochain space is zero, so the last pair supplies no source. For $(1,1)$ the terminal space and one-excess cohomology are both zero. This follows from the pure-current calculation at dotted degree four. Both of these statements are independent of any quantum kernel.

For $(2,0)$ the consecutive cochain dimensions are $(0,3,4)$. The integer three-by-three minor above proves injectivity at the middle term. For $(2,2)$ the dimensions are $(6,7,6)$. Two explicit predecessors give a minor $-294982331904$. Five independent images into the six terminal cochains leave one Cartan-visible direction. Thus $7-2-5=0$ at the middle term. To construct the two predecessors, use a third jet and five scalars, alternate separately in both four-slot columns, and contract either two invariant cubic forms or two nested brackets. The replay differentiates the third jet by the three Jacobi terms before evaluating the probe. [The two rational predecessors](plane_predecessors.py) · [Minor](plane_predecessors.json)

For the terminal evaluation, take the exterior product of the first three or four colour vectors. On each pair of exterior slots apply the colour bracket followed by its metric transpose, with one common integer clearing of denominators, and sum over pairs. Pair successive powers of this operator with the original exterior vector. The zeroth power gives its Gram determinant. Each positive power has an explicit classical primitive since the first application contains a scalar commutator. The powers produce respectively four and six independent highest-weight terminal polynomials. All 24 pair products are independent; an integer evaluation minor reduces to $955$ modulo $1009$. Exactly 23 of the products are boundaries. The two Gram determinants give $32$ on the specified four Cartan vectors. The rational tensor-product formula together with injectivity on the terminal product then kills every separately closed source in this product channel. [Exterior calculation](plane_products.py) · [Both evaluation minors and all factor characters](plane_products.json)

To test cancellations between channels, put $\omega=B(X_1,[X_2,X_3])$, the only terminal cochain for $(1,0)$. It has a primitive obtained by replacing the commutator with its matter fermion and alternating the three flavours. The one-excess and terminal spaces of this factor both have dimension one, so this two-term factor is acyclic. Multiplication by its nonzero scalar polynomial is injective. The terminal space in the $(3,2)$ factor has dimension 36.

The 24-product space meets the 36-dimensional space of these multiples of $\omega$ in exactly six dimensions. On the rational divisor $\omega=0$, an 18-by-18 minor of the 24 products reduces to $644$ modulo $1009$. The probes solve one linear equation for a coordinate of $X_3$ with denominator invertible at $1009$. The six products of $\omega^2$ with the six quartic-plane terminal cochains are independent and lie in the intersection, proving equality. Each of these six has the same primitive in both product channels: replace one factor of $\omega$ by its primitive. For any cancellation subtract this common primitive from one source and add it to the other. Each source now closes separately and has a primitive by the factor calculations. This excludes the entire stated source family. The joint terminal span has dimension $24+36-6=54$, including $23+36-6=53$ independent classical boundaries.

## The fifth glueball survives to perturbative infinity

Let $S$ be the native glueball and $p$ its exponent. Every incoming perturbative differential to its terminal class starts at excess one. In these two excesses the full fan supercharge is exactly $Q_0+uQ_1$. Indeed $Q_k$ for $k\ge2$ requires $3k$ flavour slots on $k+1$ letters, which together have at most $k+2$. Both supercharges vanish at the terminal excess. Every correction to a source retains excess one. The reverse shear $\bar P_0+u\bar P_1$ therefore anticommutes with the entire differential here. Its coefficient at order two is $\{Q_1,\bar P_1\}=0$ by the five-slot test above. All its higher coefficients vanish. The order-one Ward check needs pairs of at most three Taylor slots, independently of $p$.

For the incoming image one may use the two-term complex of excesses one and zero. A terminal class has an incoming image at order $r$ precisely when an excess-one cochain with corrections solves the descent through order $r-1$ and has that residual modulo earlier terminal images. These equations use only the indicated two terms. Classical boundaries further back change the representatives of the sources and preserve these terminal image spaces. The shears induce their leading Taylor actions on the graded pages of this two-term complex. The even weight shift at order $r$ is $(r,r\mid-r,-r,-r)$. Thus a nonzero incoming $d_r$ to the $(p,p)$ target requires source highest weight $(p-r,p-r\mid r,r,r)$. A polynomial Taylor diagram with a positive third odd label must have at least three boxes in its second even row. Hence

$$
1\le r\le p-3,\qquad
\lambda=(p-r,p-r,\underbrace{3,\ldots,3}_{r\ \mathrm{rows}}).
$$

This bound uses the relative fan complex and the gauge Ward identity. It imposes no instanton selection rule.

At $F_4$, $p=5$, the product derivation on the first page and the survival of $S^6$ exclude an incoming $d_1$ on $S^5$. If $d_1[X]=[S^5]$, then $d_1([X][S])=[S^6]$, since $d_1[S]=0$. The only possible incoming $d_2$ has shape $(3,3,3,3)$, whose transpose is the scalar-flavour weight $(4,4,4)$ of the BMN calculation. The adjacent classical ranks give

$$
\dim H(Q_0)_{\text{eleven-letter source}}=141-80-61=0 .
$$

The Taylor highest-weight component here lies entirely in the three physical scalar directions after transposition. Exactness of the covariant Schur decomposition identifies its classical multiplicity complex with that BMN complex. The bound excludes every $r\ge3$. The zero-flavour target itself is strictly closed. Therefore $[S^5]\ne0$ on perturbative $E_\infty$ at $j=30$. Its Cartan restriction vanishes since the fifth exterior power of the rank-four Cartan is zero. This supplies a perturbative target for a physical nonperturbative test. [Classical source proof, minors and replay](../../../bmn_primary/f4_scalar/README.md) · [Allowed orders and characters](source_symmetry.json)

For $S^6$ the only orders still possible are two and three. The respective source cochain dimensions are 3,720 in $(4,4,3,3)$ and 1,862 in $(3,3,3,3,3)$. The corrected pairings on these diagrams remain open.

## A nonzero charge-one endpoint for the fifth power

Embed the charge-one instanton along the highest long root of $F_4$ and choose the orientation supporting the barred Weyl zero modes. In the native colour file the raising generator has index 45 and its $\mathfrak{sl}_2$ lowering partner is minus generator 12. The Cartan element of this triple has adjoint weight multiplicities $1,14,22,14,1$ at weights $-2,-1,0,1,2$. This gives one triplet, fourteen doublets and twenty-one singlets. Each doublet contributes one Weyl zero mode of the chosen R-flavour. In a doublet the fundamental BPST solution has an invertible two-by-two spin–colour matrix away from the centre. Denote this matrix, in the zero-mode convention used for the moduli, by $M(x)$. The opposite instanton orientation exchanges chiralities. [Fundamental Dirac zero mode, Eqs. (5.35)–(5.36)](https://arxiv.org/pdf/0802.1862#page=43)

Take the fourteen weight-one basis vectors $v_i$ and the corresponding lower vectors $[f,v_i]$. The colour form on the doublet moduli $\nu^i$ is $B(v_i,[f,v_j])$ with the original $B=\mathrm{Tr}_{26}/6$. It is alternating and has seven disjoint upper-triangular entries

$$
(19,50;-1),\ (20,49;-1),\ (22,48;2),\ (25,47;-1),\
(27,46;-1),\ (29,40;-2),\ (31,36;2).
$$

On the linear doublet zero modes the glueball is $\det M(x)$ times the contraction of this two-form with $\nu^i\nu^j$. In its fifth power the coefficient of the ordered moduli with indices $(19,20,22,25,27,46,47,48,49,50)$ is

$$
2^5\,5!\,(-1)(-1)(2)(-1)(-1)\,(\det M(x))^5
=7680\,(\det M(x))^5 .
$$

Direct exterior multiplication and the five-pair minors agree on all 21 nonzero coefficients. This verifies a bare target insertion at a fixed charge-one background. The physical supercharge matrix element still requires a perturbative source, the supercurrent, the remaining fermionic and bosonic moduli integrals, and operator matching. [Native colour and exterior replay](fifth_instanton.py) · [All 21 coefficients](fifth_instanton.json)

## One row decides the second differential

Every classical cycle in the selected 3,720-dimensional source admits a first correction. Dot obtained this surjectivity in the original $F_4$ package. The following independent calculation starts with the native colour brackets and Dot's extremal tensor seeds. The intermediate terminal cochains at degree $(5,5\mid1,1,1)$ have dimension $150+238+209+36=633$ in the spin and flavour singlet. A covector annihilating their classical boundaries must kill the brackets within each five-ghost leg and within the three alternating scalar colours. The separate exterior quotients in degree five have highest weights $(0,1,0,3)$ and $(0,0,3,0)$. In degree three the highest weight is $(0,0,2,0)$. Interchanging the two ghost legs in a spin singlet gives positive sign. The Weyl character formula gives the following multiplicities of the scalar module in the product of the ghost modules.

| Ghost modules in the displayed order | Ordered multiplicity | Exchange trace | Allowed multiplicity |
|---|---:|---:|---:|
| First with first | 23 | 13 | 18 |
| First with second | 13 | | 13 |
| Second with second | 10 | 8 | 9 |

This bounds the possible annihilator by $18+13+9=40$. On Dot's forty extremal tensors the classical split of a differentiated scalar has a minor of order forty. The same raw colour probes give determinants $194027$ modulo $1000003$ and $30888$ modulo $1000033$. The program computes every entry by Lie-algebra lowering from the seeds. The tensors are independent, exhaust the character bound, and have zero common annihilator on these probes. Hence the full classical map onto the intermediate 633 cochains is surjective over the rationals. [Independent characters, modular Ward evaluator and probes](intermediate_rank.py) · [Both minors](intermediate_rank.json)

Write $C$ for the selected source and identify its spin/flavour singlet copy in the physical cochains as above. Its classical split $Q_0|_C$ has 704 rows. For $X\in\ker(Q_0|_C)$ choose a first correction $X_1$ with $Q_0 X_1=-Q_1 X$, possible by surjectivity. The established target identity $\lambda Q_1=\eta Q_0$, the vanishing of the bare $Q_2$ on this source, and the one-dimensional nonzero target on $E_2$ give

$$
\lambda(d_2[X])=-\eta Q_1 X,\qquad
S^6\in\operatorname{im}d_2
\quad\Longleftrightarrow\quad
\operatorname{rank}_{\mathbb Q}
\begin{pmatrix}Q_0|_C\\-\eta Q_1|_C\end{pmatrix}
=\operatorname{rank}_{\mathbb Q}(Q_0|_C)+1.
$$

Earlier incoming boundaries have zero obstruction on $E_2$. Thus this test needs neither a basis of the length-twelve predecessors nor separate first-correction solves. A rational row identity certifies vanishing. A nonzero augmented minor certifies a lift once its size exceeds a proved rational upper bound on the classical rank. The augmented rank remains open.

Every product of two positive-degree $E_2$ classes has zero pairing with this target, including products projected from arbitrary weights. Indeed the derivation rule acts on one homogeneous factor at a time. For its output to have zero excess and zero flavour, the spectator must already have zero excess and zero flavour. The invariant two-gaugino ring makes any such surviving positive-degree spectator a multiple of $S^k$ with $k\ge1$. The acted-on factor would have to map to $S^{6-k}$. For $1\le k\le6$, multiplying this equation by $S^{k-1}$ would make $S^5$ a second boundary. Its perturbative survival excludes this. Linear combinations and weight projections preserve the conclusion. This leaves cycles whose constituents do not separately define products on $E_2$.

For the classical split, polarize the four columns of heights $(4,4,3,3)$ on four ordinary scalar colours. Alternate all columns in a product of adjoint traces with one marked scalar commutator. Replacing this commutator by the corresponding second jet gives the classical primitive, including when the factors are away from highest weight. A8 job 5703179 evaluated 6,656 such recipes. An order-135 minor is $76$ modulo $1009$ and $142$ modulo $1013$ with the same rows and columns. The rank lower bound applies over the rationals. On the four-dimensional Cartan, the scalar Taylor module of shape $(4,4,4,2)$ is the fourth determinant power times the dual symmetric square. The Weyl reflection representation has one invariant quadratic form. The earlier product of the two Gram determinants attains this Cartan direction. It follows that

$$
135\le\operatorname{rank}_{\mathbb Q}(Q_0|_C)\le 704-1=703.
$$

The 135 rows and the earlier 53 can overlap. [Flag evaluation and integer control](flag_boundaries.py) · [Selected recipes and points](flag_135.json) · [A8 accounting and second prime](flag_a8.json)

## Integrating one flavour of lifted instanton modes

At fixed charge-one bosonic moduli, take the linear gaugino profiles in the orientation above. The embedded triplet is orthogonal to the doublets. Its exact supersymmetric and superconformal modes of the glueball flavour enter through a two-component spinor at the insertion point. The triplet contribution to $S$, denoted $S_{\rm trip}$, is quadratic in these two odd components and squares to zero. Write $\omega$ for the doublet quadratic polynomial in the native fourteen modes, so that the doublet contribution is $\det M(x)\,\omega$. Then the profile obtained from these linear fields is

$$
S^5=(\det M(x))^5\omega^5+
5 S_{\rm trip}(\det M(x))^4\omega^4 .
$$

Here the nilpotence statement is pointwise. The exact modes at separated insertions involve different linear combinations. [Adjoint zero modes and fundamental doublets, §5](https://arxiv.org/abs/0802.1862) · [Exact-mode combination, Eq. (2.5)](https://arxiv.org/abs/1605.06346)

Use the ascending fourteen-colour order given above for Berezin integration, meaning extraction of the coefficient of the ordered product. For each listed pair $(i,j;b)$, the corresponding term in $\omega$ is $2b\nu^i\nu^j$. Let $\Lambda$ contract the pairs by applying the left derivatives $\partial_{\nu^j}\partial_{\nu^i}/(2b)$ and adding the seven results. The Pfaffian in this order is $-1024$ and the coefficient of $\omega^7$ is $7!(-1024)=-5160960$. If $P_m$ is the part of the remaining integrand of degree $m$ in these fourteen modes, integration of the displayed target insertion gives

$$
-61440\,(\det M(x))^5\Lambda^2 P_4
-20480\,S_{\rm trip}(\det M(x))^4\Lambda^3 P_6 .
$$

Indeed the coefficient for $\omega^{7-k}P_{2k}$ is $(7-k)!(-1024)\Lambda^k P_{2k}/k!$. Direct multiplication checks all $1001+3003=4004$ exterior basis monomials at $k=2,3$. This integrates fourteen modes for arbitrary remaining insertions at this order in the target fields. The physical source, supercurrent, other modes and matched operator contributions remain to be evaluated. [All-monomial verifier](fifth_berezin.py) · [Coefficients](fifth_berezin.json)

~~~sh
python3 verify_six_family.py
python3 intermediate_rank.py | diff -u intermediate_rank.json -
python3 fifth_berezin.py | diff -u fifth_berezin.json -
OPENBLAS_NUM_THREADS=1 python3 flag_boundaries.py --replay flag_135.json --check --out /tmp/f4-flag-replay.json
~~~

## The corrected second row on trace sources

At level 36 the 19 reduced trace singlets have a classical matrix of rank five after the polarized quartic identity. The corrected quantum row vanishes on its 14-dimensional kernel. Eight of the 19 individual scores are nonzero. A8 gives a rational lower bound of 542 for the classical rank in the full selected diagram. The remaining finite-colour test in the trace presentation has primitive row $(0,16,36,0,3)$. [Derivation, explicit cycles and replay](trace_obstruction.md)

## Exact-mode selection for the fifth glueball

Take the separated one-instanton Ward correlator with the radial conjugate of $S^5$, in the orientation where its ten gauginos have linear zero modes. The R-charges of the incoming source, supercurrent, and radial target are respectively $(9/2,9/2,9/2)$, $(1/2,1/2,1/2)$, and $(-5,-5,-5)$. Each of the first three zero-mode flavours has charge $(1/2,1/2,1/2)-e_I$ for its corresponding matter label $I$. Flavour four has charge $(-1/2,-1/2,-1/2)$. The minimum occupation vectors in this order are

$$
(9,9,9,0),\qquad(1,1,1,0),\qquad(0,0,0,10).
$$

Every higher occupation at the same charge adds a common nonnegative integer to all four entries. The source has $E=29/2$ and $J_L=1/2$. It contains exactly one chiral matter fermion or one covariant derivative. The charge dictionary gives fifteen flavour assignments in the former case and five in the latter. Each has the displayed source minimum from its constituent fields. [Charge dictionary](../../../../../supplement/supplement.tex) · [All 20 assignments and the complete occupation test](fifth_exact_selection.json)

The one-instanton measure requires four exact and fourteen lifted coordinates per flavour. The four-fermion action uses one lifted coordinate of each flavour and carries zero exact modes. At a fixed insertion the linear radial term supplies two exact fourth-flavour modes or zero. The source and current must therefore supply at least two or four additional fourth-flavour modes, respectively. If their total common shift is $k$, the required number of quartic vertices is $18-10-k=8-k$. The two radial terms permit at most six or four vertices, respectively. At $k=0$ the fourth-flavour Berezin integral vanishes in both terms. This argument covers the linear radial field with separated classical profiles at zero scalar vev. Higher profiles, fluctuation contractions and current contact terms enter the physical Ward matrix element. [Exact-mode combinations, Eq. (2.5)](https://arxiv.org/abs/1605.06346) · [Four-fermion action, Eq. (19)](https://arxiv.org/pdf/hep-th/0009150#page=4) · [Selection verifier](fifth_exact_selection.py)

The fourteen lifted colours form the primitive three-form representation of $Sp(3)$. Its symmetric fourth power contains one singlet. The quartic obtained by polarizing the native nested bracket restricts to $6x^2y^2$ on the embedded $SU(3)$ doublets. All 64 brackets and all 64 metric entries of that embedding agree with the three-dimensional defining matrices. For a test insertion consisting of the fifth power of the native doublet two-form in **each** flavour, the four-vertex exterior integral over all 56 lifted modes is $8135671362984345600000$. Wick contraction and direct exterior multiplication give the same integer. The physical normalization and the projected source–current insertion remain open. [Native tensor, Wick contractions and embedded group](fifth_quartic.py) · [Independent exterior implementation](fifth_quartic_direct.cpp) · [All coefficients](fifth_quartic.json)

~~~sh
python3 fifth_exact_selection.py | diff -u fifth_exact_selection.json -
python3 fifth_quartic.py --direct | diff -u fifth_quartic.json -
~~~
