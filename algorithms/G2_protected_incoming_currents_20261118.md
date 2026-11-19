# Two incoming nonperturbative channels in `G2` SYM

In the holomorphic twist, the bottom components of the stress-tensor currents $\mathbf L_{Ii}$ and $\mathbf G^I_j$ are closed under the physical supercharge. Here $I,J,K$ are the three flavour components, and $i,j,k,h$ are the two holomorphic spin components. The two-operator product at separated points is closed by the supersymmetry Ward identity. Its degree-zero cohomological OPE extends over coincidence in complex dimension two. The constant term therefore gives a closed local product with the ordinary product as its leading classical term. Holomorphic derivatives commute with the supercharge. These physical facts would exclude an outgoing lift if the archived polynomials are the classical limits of the indicated products. The component identification is tested below. [Semi-chiral ring, §2 and currents (4.81), (4.83)–(4.85)](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

## The Gram determinant and the square of the gaugino bilinear

The doubled physical charges are ordered as $(2J_1,2J_2,2q_1,2q_2,2q_3)$. At incoming charge $(1,1,4,4,4)$ the only classical singlet cohomology line has dimension $5-1-3=1$. Its representative contains 18 trace terms. With the alternating tensor in the saved spin-component ordering, the proposed physical product is

$$
\frac16\epsilon^{ij}\mathbf L_{Ii}\mathbf G^I_j .
$$

Take the two diagonal Cartan elements in the seven-dimensional representation used in the verifier. Their trace Gram matrix is $\operatorname{diag}(4,12)$. Set the first scalar and the plus gaugino jet along the first element, the second scalar and the minus derivative of the third scalar along the second element, and the matter fermion to zero. With $\eta$ the exterior coefficient of the plus gaugino, the only contributing term is

$$
\frac16(4\eta)(6)=4\eta .
$$

The saved representative evaluates to the same $4\eta$. All 217 classical letter variations at this weight evaluate to zero. This comparison uses the scalar term of (4.84) on the chosen jet, so the matter-fermion component sign drops out. This scalar-only evaluation does not test the relative sign of the matter term. On a second Cartan jet, with the plus gaugino, minus gaugino and first matter fermion in that exterior order, the archived coefficient is $-16/3$, whereas the direct substitution in (4.81) and (4.84) with both flavour orientations positive gives $+16/3$. The exact-current identification, and hence this proposed exclusion of an instanton differential, is unresolved. [Source ranks and representative](verify_g2_np_partner.py) · [Direct physical-current evaluation](verify_g2_current_product.py)

## The cube of the gaugino bilinear

At incoming charge $(1,1,6,6,6)$, equality of the two spin components forces equal derivative degrees. For common derivative degree $t$, the three flavour degrees are $3-t$ and the letter length is $8-t$. Nonnegativity gives exactly four possibilities. The classical differential $Q_0$ preserves the five degrees and increases length by one. The ranks below are lower bounds from integral evaluations at each of the primes $1{,}000{,}003$ and $1{,}000{,}033$. Nilpotence gives the upper bounds in the last column.

| $t$ | Length | Cochains | Incoming rank | Outgoing rank | Cohomology upper bound |
|---:|---:|---:|---:|---:|---:|
| 0 | 8 | 15 | 6 | 9 | 0 |
| 1 | 7 | 48 | 28 | 20 | 0 |
| 2 | 6 | 30 | 17 | 11 | 2 |
| 3 | 5 | 3 | 0 | 3 | 0 |

At length six take the following two proposed classical limits of exact current products, using the alternating flavour and spin tensors of (4.83)–(4.85) and the saved normalization of $S$.

$$
\frac{S}{6}\epsilon^{ij}\mathbf L_{Ii}\mathbf G^I_j ,
\qquad
\epsilon^{IJK}\epsilon^{ik}\epsilon^{jh}
\mathbf L_{Ii}\mathbf L_{Jj}\partial_k\mathbf L_{Kh}.
$$

Direct substitution of the printed components with the above orientations into two diagonal jets gives the following submatrix, with a column per product. The first row selects the exterior coefficient with four gauginos and one matter fermion. The second selects a three-gaugino coefficient with one scalar derivative. All 994 tested classical letter variations vanish on these jets.

$$
\det\begin{pmatrix}-64&0\\16&576\end{pmatrix}=-36{,}864 .
$$

The closed, archived superfield representatives instead give

$$
\det\begin{pmatrix}64&0\\16&576\end{pmatrix}=36{,}864 .
$$

This last test certifies two rational classical classes and saturates $30-17-11=2$. The relative sign on the matter row is unresolved as a physical dictionary; neither determinant computes an instanton differential or proves that all incoming sources have exact completions. [Independent superfield and component calculation](verify_g2_superfield_sign.py) · [Both signed matrices](../data/N4_G2_superfield_sign_20261118.json) · [Two-prime classical ranks](../data/N4_G2_S3_incoming_currents_20261118.json).
