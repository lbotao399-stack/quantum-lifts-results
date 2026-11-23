# The full classical two-colour block of R-007

In derivative-free classical $SU(3)$ cohomology, the restriction of the 6,139-term level-30 operator $F$ to the $S(U(2)\times U(1))$ block is

$$
[F|_{\rm block}]=-\frac{5184}{5}\,i\, h_1^3[O_0].
$$

Here $O_0$ is Dot’s frozen 117-term $SU(2)$ operator in unscaled Pauli coordinates. The first scalar has central part $\operatorname{diag}(h_1,h_1,-2h_1)$. The other two scalar central parts and all three abelian matter fermions are arbitrary. The abelian field-strength letter is absent from the input. On the formal classical neighbourhood of any semisimple joint eigenvalue collision of type $2+1$, the ideal of scalar dressings annihilating $F$ is the **reduced collision ideal**. In particular a single power of any projected scalar discriminant annihilates the germ there. The bound nine for a global polynomial primitive and the nonperturbative current question have different hypotheses. [Global support](collision_support.md) · [Dot’s original cubic coefficient](https://github.com/lbotao399-stack/quantum-lifts-results/blob/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-exact-input/block_witness.py)

## Computation with abelian spectators

First take the scalar background in the first flavour only. In every scalar or fermion word expand both the two-colour matrix trace and the one-colour trace. Put the abelian fermions on the right in flavour order with their exterior signs. The 16 nonzero coefficients have 8,790 monomials in the traceless $SU(2)$ fluctuations. Each has the common imaginary unit from the Pauli trace. Fifteen coefficients are literal classical boundaries. Their primitives use a total of 1,012 rationally weighted orbits under rotations of the coordinate axes of the Pauli matrices. The remaining coefficient, with three first-scalar shifts and zero central fermions, is exactly $-(5184/5)O_0$ after removing the common $i$. Its separating functional gives $31104 i$ with the original factor restored. [Full expansion, elimination and direct differentiation](levi_block.py) · [All rational primitives](levi_primitives.json)

For example the coefficient linear in the scalar shift has 759 monomials and a primitive with 136 orbits. In the same charge block the polynomial $(\tfrac12\operatorname{Tr}_2 (X^1)^2)O_0$ has 351 monomials and a primitive with 57 orbits, where $X^1$ is the fluctuating traceless first scalar. The computation starts with 13,734 raw predecessors at length eight. Projection to the finite rotation group gives 568 orbit columns. The classical image on them has rank 439. For each reported primitive the verifier expands the rotations and compares **all** rational monomial coefficients of its classical derivative to the prescribed raw target.

Colour invariance of each target is tested with all three infinitesimal rotations. To obtain a fully gauge-invariant primitive from any archived primitive of length $L$, apply the polynomial
$$
\prod_{\ell=1}^{L}\left(1-\frac{C}{\ell(\ell+1)}\right),
$$
where $C$ is the positive $SU(2)$ quadratic Casimir with eigenvalue $\ell(\ell+1)$ on integral spin $\ell$. A tensor product of $L$ adjoint vectors has maximum spin $L$. The polynomial fixes the invariant target and commutes with $Q_0$.

The six flavour root operators annihilate $O_0$ as a literal polynomial. Flavour covariance therefore gives
$$
\left[\tfrac12\operatorname{Tr}_2(X^iX^j)\,O_0\right]=0
\qquad (1\leq i,j\leq 3)
$$
by rotating and polarizing the scalar square. Here $i,j$ label the three scalars. The stabilizer of the highest component $F$ contains the determinant-one shears adding the first scalar to either of the others. Whenever $h_1\ne0$ these send a general scalar triple of central shifts to $(h_1,0,0)$. They act invertibly on the three central fermions, which were arbitrary in the replay. This proves the displayed block identity on a dense open set of background triples. In each fluctuation degree the classical block differential is independent of the central variables. Passing to its complex vector-space cohomology coefficientwise therefore extends the polynomial identity to every triple. [Flavour action on the parent](discriminant_pencil.py)

## Formal neighbourhood of a separated two-colour block

Fix a scalar combination with nonzero gap between the two-colour eigenspace and the remaining colour. A formal off-diagonal conjugation block-diagonalizes this scalar degree by degree. The linear equation for each conjugation coefficient divides its two off-block matrix components by the two nonzero opposite gaps. The conjugation depends only on a closed scalar and commutes with $Q_0$. The residual gauge group is the complexified block group.

For this choice of first scalar, in the ideal of off-block fluctuations the linear part of $Q_0$ pairs the other two scalar entries with their matter-fermion partners and pairs the remaining matter fermion with the field-strength letter, with invertible coefficients given by the gap. The anticommutator with the inverse linear pairing counts the off-block slots. Divide by this positive count on the kernel of restriction to the block. The remaining terms of the classical differential raise total fluctuation degree. Subtract the primitive at the lowest residual degree and repeat. This gives a contracting homotopy on that kernel in the completed complex. The construction commutes with the residual colour action, or one may project each finite degree to invariants. Restriction thus detects classical exactness in this formal neighbourhood. This argument concerns the derivative-free classical complex. [Native commutators and regular scalar homotopy](discriminant_ward.md)

Locally the third joint eigenvalue is separated. The two remaining joint eigenvalue vectors have a central part and opposite fluctuation parts, and their interchange changes the sign of the latter. The ideal of their reduced coincidence is generated in the invariant scalar ring by the six pairwise products of these three fluctuations. Their gauge-invariant lifts to the two-colour block are $\tfrac12\operatorname{Tr}_2(X^iX^j)$. Higher symmetrized scalar traces there are polynomials in these quadratic traces. Their vanishing on $[O_0]$ gives one inclusion of ideals. For the reverse inclusion, remove these quadratic terms from a dressing. Its remaining central series multiplies the separating functional to give $31104 i\,h_1^3$ times that series. The ring of formal central scalar series over the complex numbers is an integral domain, including when $h_1$ has zero background value. Every annihilating series consequently has zero restriction to the collision component.

For a scalar projection with central part $h$ and traceless part $X$, put $s=\tfrac12\operatorname{Tr}_2 X^2$. The three projected eigenvalues give
$$
\Delta=4s(9h^2-s)^2 .
$$
This gives the single-power annihilator at a $2+1$ background. The formal slice proof requires a separated third colour. The triple collision at the conformal origin and the physical instanton differential remain open.

Run the following in this directory. Both commands reconstruct the native operator. The second also solves every primitive afresh.

~~~sh
python3 -B levi_block.py | diff -u levi_block.json -
python3 -B levi_block.py --rebuild | diff -u levi_block.json -
~~~
