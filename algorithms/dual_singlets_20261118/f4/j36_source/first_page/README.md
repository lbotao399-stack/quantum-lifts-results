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
