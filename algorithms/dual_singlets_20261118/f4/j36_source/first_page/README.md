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
