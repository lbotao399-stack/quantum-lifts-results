# Three ordered Albert exchanges in the second source for the seventh F₄ glueball

In the $(5,5,3,3)$ Taylor diagram, take the same two alternated Albert endpoints as in the [two-exchange calculation](../seventh_exchange_kernel.md). Allow up to three ordered pair exchanges on four open colour slots. The number of words is $1+6+36+216=259$. Either endpoint can carry the free second jet. The classical split on the resulting 518 formal recipes has rational rank 191. Every closed combination in this family is a classical boundary or an identically zero cochain.

## Constraints on the native endpoints

Use slot numbers $0,1,2,3$ as in the native evaluator. A pair on a single Albert cubic, with the third index fixed, represents a symmetric Jordan multiplication matrix $L$ in the invariant form. The 52 native Lie matrices are skew. All $26\cdot52^2=70304$ traces $\operatorname{Tr}_{26}(LUV)$ vanish. Also $\operatorname{Tr}_{26}(LUVW)=0$ for three Lie matrices $U,V,W$. Indeed transposing reverses the three factors with a minus sign, and

$$
UVW-WVU=[U,V]W+V[U,W]+[V,W]U .
$$

Each term on the right is a product of two Lie matrices. This identity applies to arbitrary native Lie inputs, including a freely assigned marked colour and the commutator replacing it. [Native matrices and tests](verify.py)

On the right endpoint, pair $(0,2)$ carries two actions in the first alternation and one in the last. Its trace vanishes by the three-factor identity. On pair $(1,3)$ the first alternation gives a skew matrix orthogonal to the Lie algebra. The next alternation gives a symmetric matrix of trace zero. [Pair calculation](../seventh_exchange_kernel.md)

On the left, pair $(0,1)$ carries the mark and a scalar in the first alternation. In the final alternation the trace with two further actions pairs their commutator with this intermediate matrix, giving a trace of $L$ times three Lie matrices. On pair $(2,3)$ the first alternation again gives a skew matrix orthogonal to the Lie algebra. The final alternation acts once on this pair. Its trace vanishes. Its pairing with a Lie matrix reduces to the pairing of the skew intermediate matrix with a product of two Lie matrices. Only their commutator contributes, giving zero. Thus both left pairs are trace-free, and pair $(2,3)$ is also orthogonal to the adjoint representation.

The remaining constraint contracts slots $(1,2,3)$ of the right endpoint with one symmetric Albert cubic. Fix the complementary pair of flavour labels of its mark on slot $0$. On the other three slots each of the two alternating groups acts by a separate alternating triple of Lie matrices. The mark commutes with all these slot actions. Permuting the three slots multiplies each alternating operator by the same permutation sign. Their product therefore commutes with the permutation. Since the output contraction is symmetric, only the fully symmetric part of the input corner on those three slots contributes.

In the native integer tensors the sum over the six permutations of the upper corner equals four times the sum over the three inverse-metric pairings, with the stored metric factor $3h^{-1}$. The code checks this equality on all $26^4$ entries. Consider any one of these pairings and fix both actions on the strand paired with slot $0$. Antisymmetrizing the two inner Lie actions on the other metric pair gives their commutator $K$. Antisymmetrizing the two outer actions gives, up to the common sign, $UKV-VKU$. The preceding identity writes this as three products of two Lie matrices. The remaining cubic tests it against $L$ and gives zero. This proves the cubic constraint for both choices of marked colour.

## All closed recipes in the family

The verifier inserts the six native exchange matrices on the four tensor slots. It imposes the preceding constraints on both endpoints. On the zero-weight subspace of dimension 5512 the last cubic contraction has rank two in the right parametrization. Eliminating two pivot coordinates gives 5510 spanning test columns. For each of the 68 supplied independent rational relations it tests the full operator between the constrained endpoints on every such column. Integer bounds certify each multiplication. For larger intermediate sums, coprime moduli whose product exceeds twice the entry bound certify exact equality over the integers.

All operators and all endpoint subspaces are $F_4$-equivariant. Every irreducible $F_4$ representation contains weight zero. Indeed the weight and root lattices coincide, and zero is in the convex hull of every Weyl orbit. The weight-polytope theorem gives zero as a weight. Consequently a nonzero equivariant map between these subspaces would act nontrivially at weight zero, since its image is a direct sum of irreducibles. This proves the 68 relations on the full colour tensors, for either free mark. Their coefficient rows already have rank 68 modulo 1009. [Relations](relations.json) · [Complete native verifier](verify.py)

Splitting either mark gives the same outgoing column. The same ordered minor of order 191 gives $663\bmod1009$ and $49\bmod1013$. The 68 relations give the matching upper bound $259-68=191$. Each of the 259 differences of the two markings is the classical derivative of its double-mark cochain. Together with the 68 identities on one marking these are independent as formal coefficient vectors. Hence

$$
518-191=327=259+68 .
$$

The contribution of this family to the incoming $d_2$ is zero. The entire selected source has dimension 62866, and the remaining independent sources require the [mixed-jet Ward test](../seventh_ward.md). [Native ordered evaluations](seventh_cubic_paths.py) · [A8 runs and complete result](accounting.json)

With NumPy and SciPy, the following verifies all rational identities and both saved minors. The two final commands regenerate both evaluation matrices from the native colour tensors.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 -B verify.py --workers 4
python3 -B seventh_cubic_paths.py --prime 1009 --workers 4 --output /tmp/third-exchange-1009.json
python3 -B seventh_cubic_paths.py --prime 1013 --workers 4 --fixed /tmp/third-exchange-1009.json --output /tmp/third-exchange-1013.json
~~~
