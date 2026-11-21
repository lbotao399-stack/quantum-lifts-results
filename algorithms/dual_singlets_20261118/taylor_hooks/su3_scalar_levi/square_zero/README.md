# The first correction on every fixed three-colour scalar background

For the archived $SU(3)$ sources at levels 32 and 34, fix the three undifferentiated traceless scalar matrices and leave all other matrix entries free. The first image is a $Q_0$ boundary in this component complex on every such fibre. The unrestricted, gauge-invariant polynomial correction remains open at both levels. [Level-32 source and first image](../../su3_level32/README.md) · [Level-34 source and first image](../../../../su3_j34_20261118/README.md)

## A cyclic scalar

Separate the first odd coordinate in the source superfield as $\Psi=A+\theta^1 B$. Here $A$ is odd with zero constant coefficient and $B$ is even with constant coefficient $S=\phi^1(0)$. The remaining two even and two odd Taylor variables are free. The convention $Q_0\Psi=\Psi^2$ gives

$$
Q_0 A=A^2,\qquad Q_0 B=[A,B],\qquad Q_0 S=0.
$$

Choose a colour vector $v$ such that $\det(v,Sv,S^2v)\ne0$. With the three columns $T=(v,Bv,B^2v)$, Cayley–Hamilton gives

$$
BT=TM,\qquad
M=\begin{pmatrix}0&0&\det B\\1&0&\operatorname{tr}(B^2)/2\\0&1&0\end{pmatrix},\qquad Q_0 M=0 .
$$

On this chart $T$ is invertible as a Taylor series. In the coordinate frame $v=(1,0,0)^{\mathsf T}$ the eight scalar coordinates $Bv,B^2v,\operatorname{tr}(B^2)/2,\det B$ have Jacobian determinant $(\det T)^2$ in the eight row-ordered traceless matrix coordinates. The inverse is $B=TMT^{-1}$. At the principal nilpotent companion matrix $\det T=1$. The symbolic calculation also reproduces all 180 ordered splits of the 31 letters needed for the two first images. [Polynomial identities and signs](cyclic.json) · [Replay](../../../su3_cyclic_scalar.py)

The odd matrix $T^{-1}AT-T^{-1}Q_0T$ commutes with $M$. A matrix commuting with a cyclic matrix is a polynomial of degree at most two in that matrix, as its action on the cyclic vector determines all three columns. The odd polynomial squares to zero. Its $Q_0$ variation is its square by the change-of-frame identity. Its two nonscalar coefficients are therefore closed. Its trace is $-\operatorname{tr}(T^{-1}Q_0T)$ since $A$ is traceless. These data reconstruct $A$ uniquely. In every positive Taylor degree the six orbit coordinates in $Bv,B^2v$ pair with their six classical variations. The two coefficients of the characteristic polynomial and the two nonscalar odd coefficients have zero differential. After deleting the six contractible pairs in each degree, the resulting algebra is free supercommutative over the scalar chart. In particular multiplication by the discriminant of $S$ is injective there, including on a chart containing a principal nilpotent scalar.

For a Cartan-zero classical trace cocycle some power of that discriminant kills its class. Freeness makes its restriction to every cyclic chart exact. This applies to both first images and to a closed level-34 second residual if a global first correction exists. The earlier Levi calculation treats each nonnilpotent scalar with a simple eigenvalue, including the diagonalizable double-root case. Flavour transformations commute with the classical restriction. A remaining fixed triple must have every scalar combination nilpotent of rank at most one. [Discriminant annihilator](../../../../perturbative_cartan_targets.md) · [Levi calculation](../README.md#arbitrary-matrices-above-a-scalar-with-a-simple-eigenvalue)

## Rank-one scalar pencils

If two rank-one matrices have independent image lines and independent covectors, their sum has rank two. A linear space of matrices of rank at most one therefore has a common image line or a common kernel plane. Indeed, once two matrices have common image and independent covectors, any third with different image would need its covector proportional to both. The dual argument gives the common kernel case. Tracelessness then makes every pairwise scalar product zero. The family with common image has dimension $2+3\cdot2=8$ in the 24-dimensional space of scalar triples. The common-kernel family has the same dimension and is obtained by negative transposition. Both archived first images have even colour parity.

Fix the common image as the first coordinate line. Each scalar has arbitrary entries in positions $12,13$ and zero elsewhere. All other matrices are generic and traceless, with independent exterior coordinates for their odd entries. Expansion over the rationals gives the following counts. The raw polynomial hashes at level 34 agree between the local replay and A8 job 5686273. [Level 32](j32.json) · [Level 34](j34.json) · [Cluster receipt](a8.json)

| Level | Input trace words | Words with no adjacent scalars, including cyclic adjacency | Nonzero matrix monomials |
|---:|---:|---:|---:|
| 32 | 331201 | 2585 | 1620 |
| 34 | 517468 | $24756+16294=41050$ | $31852+9538=41390$ |

At level 32 each surviving word contains five scalars, four undifferentiated matter fermions and the two first gauginos $\lambda_+,\lambda_-$. These letters are $Q_0$ closed on this locus. Some trace contains more scalars than matter fermions. If it contains two scalars separated only by gauginos, commute one scalar through the gap. Each commutator is the signed $Q_0$ variation of the corresponding first scalar derivative. The final adjacent scalar product vanishes. The other case is a trace of one scalar $S$ times a word $W$ in gauginos. For the projector $P$ onto the fixed image line along the other two coordinates, $PS=S$ and $SP=0$, so

$$
\operatorname{tr}(SW)=\operatorname{tr}(P[S,W]).
$$

The rational first correction has 2447 trace terms, 68 of which use $P$. Its classical variation cancels all 2585 trace coefficients after the scalar and projector relations. This calculation is polynomial in all six scalar entries. [Correction](j32_correction.json.gz) · [Exact replay](j32_replay.json) · [Construction](../../../su3_square_zero_primitive.py)

At level 34 fix the scalar triple. For each spin, the variations of the first scalar derivatives give the entries of $[\lambda_\pm,\phi^i]$ up to the archived signs. Select independent such components, using their derivative parents and these variations as linear coordinates. Setting the selected pairs to zero restricts the two gauginos to the joint centralizer. This quotient induces an isomorphism on $Q_0$ cohomology. To see this, filter by the number of unfixed letters. All constant variations vanish for a commuting scalar triple. The linear differential splits off the selected contractible pairs, so the induced map on the cohomology of the associated graded free supercommutative algebra is an isomorphism. At each fixed spacetime bidegree and classical cochain degree the filtration is finite. In fact the sum of the odd Taylor degrees minus the letter count gives the latter grading up to sign. Negative contributions require positive spacetime degree, and every unfixed letter of zero grading does as well. This proves the assertion for the full differential. Applying the centralizer restriction to the gauginos alone suffices for the polynomial evaluations below.

For two independent scalars the joint traceless centralizer is spanned by the matrix units $E_{12},E_{13}$. Restricting both gauginos to it gives zero in every matrix coordinate. For a one-dimensional nonzero scalar span take $\phi^i=x_i E_{12}$. A traceless matrix commuting with $E_{12}$ has entries $21,23,31$ zero and equal entries $11,22$. Restricting the gauginos by these four equations leaves **five** monomials. The field strength $f$ has

$$
Q_0 f=[E_{12},-x_1\Psi_{23}+x_2\Psi_{13}-x_3\Psi_{12}],
$$

where the subscripts on $\Psi$ specify its two odd Taylor directions. The two entries $(Q_0 f)_{13}$ and $(Q_0 f)_{32}$ give a five-row, two-column system of rank two. Its correction coefficients are both $-5012280$. Its residual has zero terms, identically in $x_1,x_2,x_3$. Substitution on each of the three charts $x_i\ne0$ independently gives zero. At the zero triple the first image already vanishes. Negative transposition transfers the calculation to common kernel. [Two-term polynomial and three rational checks](j34_field_strength.json) · [Exact system](../../../su3_square_zero_field_strength.py)

The fixed-scalar calculations permit colour frames and corrections depending on the fibre. For a global lift test one still needs a single invariant polynomial first correction or an unrestricted obstruction functional. Pointwise exactness alone cannot supply the polynomial, since a torsion class may vanish on every residue field. The level-34 second equation remains conditional on that first correction.

## Replay

From the repository root, with Python and SymPy installed, run the following commands. The direct expansions use the local rational input streams and the vendored Taylor alphabet.

~~~sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_scalar_levi/square_zero"
python3 "$d/su3_cyclic_scalar.py" --out /tmp/n4-cyclic.json
cmp /tmp/n4-cyclic.json "$p/cyclic.json"
for j in 32 34; do
  python3 "$d/su3_square_zero_image.py" --level "$j" --workers 4 --centralizers --out "/tmp/n4-square-$j.json"
  cmp "/tmp/n4-square-$j.json" "$p/j$j.json" || exit 1
done
python3 "$d/su3_square_zero_primitive.py" --correction "$p/j32_correction.json.gz" --check /tmp/n4-square-32-replay.json
cmp /tmp/n4-square-32-replay.json "$p/j32_replay.json"
python3 "$d/su3_square_zero_field_strength.py" --input "$p/j34.json" --out /tmp/n4-square-34-field.json
cmp /tmp/n4-square-34-field.json "$p/j34_field_strength.json"
~~~
