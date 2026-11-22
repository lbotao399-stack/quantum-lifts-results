# The first exceptional colour tensor and the level-36 source

For eight separately labelled adjoint inputs of $F_4$, the space of invariant tensors has dimension **2891**. Products of defining $26$-dimensional traces span **2877** dimensions. Their quotient is one copy of the $14$-dimensional irreducible representation of the slot permutation group $S_8$ with Young diagram $(4,4)$. One Albert-algebra contraction and its slot permutations supply the entire quotient. [Colour calculation](degree_eight_colours.py) · [All selected words and rational-rank witnesses](colour-certificates.json)

In the projected level-$36$ source for the sixth glueball power, take two quartic Albert covariants contracted with one another and multiply by an arbitrary separate invariant of the remaining letters. Include every placement of the one second jet. This gives exactly **four independent source cochains**. Their four classical images are independent modulo the 61 earlier images. This family contains **zero nonzero classical cycles**. Other colour contractions and the incoming second differential on the glueball remain open. [Source calculation](quartic_arm_complex.py) · [Two-prime minors](arm-certificates.json)

## Eight labelled adjoint inputs

Write $J(X_1,X_2,X_3,X_4)$ for the projection of the matrix product of four adjoint inputs, acting as derivations on the Albert algebra, onto its $26$ traceless Jordan multiplication operators. The pairing on those operators is the matrix trace. The derivations kill the identity line, so their ordinary matrix traces are the defining traces. For each multiplication operator its trace against a product of two derivations vanishes. Exchanging adjacent derivations in a cubic inserts a Lie bracket and reduces to this quadratic identity. Metric transposition kills the resulting symmetric cubic. Exchanging two inputs in a quartic now inserts a cubic, so $J$ is symmetric. The replay tests the quadratic statement on all $26\cdot52^2$ components. [Albert products, derivations and trace form](../../f4_colour.py)

In two colour variables $X,Y$, polarize the scalar product of $J(X,X,X,X)$ and $J(Y,Y,Y,Y)$ by alternating between the rows in each of the four columns. Using repeated arguments as superscripts inside $J$, this gives

$$
2\langle J(X^4),J(Y^4)\rangle
-8\langle J(X^3,Y),J(X,Y^3)\rangle
+6\langle J(X^2,Y^2),J(X^2,Y^2)\rangle .
$$

For each of the 22 partitions of eight, the Weyl character of the adjoint and the symmetric-group character give the exact invariant multiplicity. Column alternation on products of defining traces attains this multiplicity on every diagram except $(4,4)$. On $(4,4)$ it gives five of six directions. Adjoining the displayed tensor gives a minor of order six with determinants $942890\bmod 1000003$ and $75741\bmod 1000033$. The bound of five is exact. At two-variable weight $(4,4)$ there are 33 formal traces and at least 17 independent polarized quartic identities. The other two-row diagrams already supply $3+1+6+1=11$ independent weight vectors, leaving at most $33-17-11=5$. Summing each multiplicity times the dimension of its slot representation gives 2891 and 2877. Independently, 2891 is the sum of squared multiplicities in the fourth tensor power of the adjoint. [Full block certificates](colour-certificates.json)

The $(2,2,2,2)$ weight space of this single copy has dimension three. The coefficient of degree four in $(1+x+x^2)^4$ minus its coefficient of degree three gives $19-16=3$. Independent direct evaluation gives ranks $224$ and $227$ in this weight. The native defining traces span every invariant at degrees four, five, six and seven, of dimensions $5,16,80,436$. [Direct weight and lower-degree checks](mixed_colour_relations.py) · [Weight minors](colour-certificates.json)

## Where the four level-36 directions come from

Use the four auxiliary Taylor directions in the level-36 projection. Their ordinary first jets are even colour variables. The second jet is antisymmetric in its two Taylor labels and has classical differential $Q_0$ equal to the bracket of its first jets, with the native orientation. The Young diagram is $(4,4,4,2)$. Substituting the physical divided jets and extracting the six scalar-to-dotted parameters and two marked parameters gives the physical dotted and flavour degrees $(4,4,2,2,2)$. Pullback in the Taylor variables commutes with $Q_0$. On terminal first-jet tensors the extraction is injective on this diagram by the nonzero highest-weight coefficient 50400. [Projection and native tree signs](quartic_sources.md#physical-closure-from-ordinary-colour-polynomials)

* **Second jet outside the two quartics.** The symmetric pair of quartic Taylor tensors decomposes into shapes $(8),(6,2),(4,4)$. Only $(4,4)$ fits in $(4,4,4,2)$. Removing its top two rows leaves the ordinary diagram $(4,2)$, with multiplicity one. To compute the residual source it therefore suffices to use two Taylor directions with first jets $X,Y$ and one second jet. At weights $(4,2)$ and $(5,1)$ the remaining first jets have multiplicities $(3,1)$ and $(4,0)$. The respective invariant dimensions including the single second jet are 2 and 0. Their difference gives an upper bound of two.
* **Second jet inside one quartic.** The five remaining first jets have invariant multiplicities two in shape $(3,1,1)$ and one in $(2,1,1,1)$. Add a horizontal four-strip for the unmarked quartic, a horizontal three-strip for the ordinary inputs of the marked quartic, and a vertical two-strip for its second jet. Pieri's rule permits exactly one path,

$$
(3,1,1)\longrightarrow(4,3,1,1)\longrightarrow(4,4,3,1)
                  \longrightarrow(4,4,4,2).
$$

This gives an upper bound of two more sources. For a finite lower bound choose the product of a two-letter and a three-letter trace, and then a single five-letter trace, in the five residual slots. For each selected source the replay replaces the second jet by its matrix commutator before alternating the four columns of heights $4,4,3,3$. The two exterior placements and two interior placements give

| Prime | Earlier rank | After the two exterior placements | After the two interior placements | Determinant of order 65 |
|---:|---:|---:|---:|---:|
| 1000003 | 61 | 63 | 65 | 951023 |
| 1000033 | 61 | 63 | 65 | 438599 |

All denominators are invertible at the displayed primes. Either determinant supplies the rational lower bound. The character upper bound is $61+2+2=65$. Consequently adding this whole projected family leaves the image of the combined source kernel in classical cohomology equal to that of the original family, which has dimension zero. [Exact proof for the earlier 19 closed combinations](quartic_sources.md) · [Minors and tree orientations](arm-certificates.json)

A second evaluation reads the cells directly from the four diagrams in the Pieri path. At the first colour assignment it sums all $20736$ signed column permutations before grouping equivalent words. Its 22 five-letter trace values agree with the grouped computation at both primes. [Literal check](quartic_arm_complex.py)

There is a direct check on the three original exterior placements. In their order in the replay, their columns satisfy $12\,\mathrm{column}_3+24\,\mathrm{column}_2-\mathrm{column}_1=0$. Polarize $12\operatorname{Tr}_{26}(X^4)=(\operatorname{Tr}_{26}X^2)^2$ on the two brackets and the two single letters. Its two crossed pairings contain the three-form $\operatorname{Tr}_{26}([X,Y]Z)$. After column alternation all 64 formal coefficients of the crossed tree terms cancel. Replacing the first bracket by its antisymmetric second jet gives cancellation in all 444 formal source coefficients. This relation already holds on the indicated representatives. [Three columns](jordan_source.py) · [Integer cancellation](quartic_arm_complex.py)

## Size of the remaining calculation

With every colour contraction admitted in this auxiliary Taylor diagram, the classical cochain dimensions at lengths $12,13,14$ are $10007,3720,704$. The length-13 count was also reconstructed by removing a vertical two-strip and inserting the adjoint character. Every classical image at length 14 vanishes on commuting colours. The product of the size-four and size-three Gram determinants of the auxiliary colours has this diagram. Its value on the first four native Cartan elements in the archive's metric $B=\operatorname{Tr}_{26}/6$ is $4\cdot8=32$. Hence the full outgoing rank in this diagram lies between 65 and 703. Determining the incoming rank at length 13 and the second differential on surviving sources requires further calculation. [Character and Cartan computation](quartic_arm_complex.py)

## Replay

Use Python, NumPy, SciPy and python-flint. Supply Dot's original zip with SHA-256 `0a06b306a12abbe77cadbcfab616c08c579d82086e9798aece27b8c83a0f7b52`. From this directory run

~~~sh
for p in 1000003 1000033
do
  python3 degree_eight_colours.py --prime "$p"
  python3 mixed_colour_relations.py /path/to/F4-j36-d2-verification.zip --paired-eight --prime "$p"
  python3 quartic_arm_complex.py /path/to/F4-j36-d2-verification.zip --prime "$p"
done
python3 quartic_source_test.py /path/to/F4-j36-d2-verification.zip
~~~
