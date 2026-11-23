# Classical boundaries in the source of the sixth F₄ glueball

In the relative holomorphic jet complex, the length-thirteen source of the incoming second differential to $S^6$ at $j=36$ has **at most 15 classical cohomology directions** in Taylor diagram $(4,4,3,3)$. The source cochain dimension is 3,720. The classical differential out of it has rational rank 703. The certificates below bound the incoming classical rank from below. [Outgoing rank and character](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank) · [Reduction to the corrected quantum row](../first_page/README.md#one-row-decides-the-second-differential)

| Actual length-twelve predecessors included | Certified incoming rank | Ordered minor modulo 1009 | Upper bound on classical source |
|---|---:|---:|---:|
| One second jet on each side of a colour contraction | 2398 | 676 | 619 |
| Add a third jet or both second jets inside one word, and further split pairs | 2748 | 873 | 269 |
| Add projected subwords and all five endomorphism channels | 2879 | 431 | 138 |
| Connect two words by a defining Jordan leg and an adjoint leg | 3002 | 833 | 15 |

A separate row-pivoting elimination on the same ordered minors gives determinants 174 and 100 modulo 1013 for orders 2748 and 2879 respectively. Each nonzero modular minor already proves the indicated lower bound over the rationals. All reductions use the same rational maps, rational colour representation and integer points. Multiplicities, chosen row and column numbers, and determinants are in the three certificates. The order-3002 minor has determinant 615 modulo 1013. The leading classical kernel has dimension $3720-703=3017$.

## Construction

Polarize with four ordinary scalar colours and alternate the four Taylor columns of heights $4,4,3,3$. The transposed highest weight is $(4,4,4,2)$. In this polarization the scalar colours $X_i$ are even, the second jets $Y_{ij}$ are odd and alternating, and the third jets $Z_{ijk}$ are even and alternating. In the exterior-jet convention the classical differential is

$$
Q_0X_i=0,\qquad Q_0Y_{ij}=[X_i,X_j],\qquad
Q_0Z_{ijk}=[X_i,Y_{jk}]+[X_j,Y_{ki}]+[X_k,Y_{ij}].
$$

Starting from two second jets or one third jet in a colour contraction gives a length-twelve predecessor. Apply the graded Leibniz rule, alternating each column onto its first three or four scalar labels. This gives a genuine highest-weight boundary linear in one remaining second jet. In the 26-dimensional representation we join matrices by the trace, the trace pairing on the adjoint, the trace pairing on the Albert multiplication operators, the product of the scalar traces and the invariant metric transpose pairing. With $H$ the invariant form of the representation, the transpose used here is $H^{-1}A^{\mathsf T}H$. The symmetric and alternating matrix spaces decompose as

$$
\operatorname{Sym}^2(26)=1\oplus26\oplus324,\qquad
\Lambda^2(26)=52\oplus273 .
$$

Inside a word, subtracting the 1 and 26 projections from the symmetric part gives the 324 projection. Subtracting the adjoint projection from the alternating part gives the 273 projection. Invariance of the Albert cubic and the metric is checked against every one of the 52 generators. Repeated or proper disjoint projected subwords produce further invariant contractions. The formulas use the fixed native representation. [Integer colour data](representation26.json.gz) · [Invariant cubic](global_boundary_cubic.json)

## Computation and verification

The expanded evaluator splits every differentiated jet into words before alternating. The second evaluator forms the six scalar commutators and four Jacobi matrices once per colour point. It computes just the prescribed columns of the earlier minors and evaluates the new predecessors with a signed subset recursion. On each prime all thirteen expanded banks agree with the direct substitution on independent selected entries. An additional verifier enumerates complete Taylor permutations literally and differentiates the predecessor by the product rule. Its nonzero controls include every terminal colour channel. Specializing every $Y_{ij}$ to $[X_i,X_j]$ annihilates every boundary bank by $Q_0^2=0$. The finite-field products use binary64 only within the exact integer range and reduce between matrix products. The selected columns in the last stage come from a fixed nonzero-column screen; every rank quoted is obtained on separate colour assignments.

From this directory with Python 3 and NumPy, set the BLAS thread counts to one. Set the worker count in the full replays to the number of available processes. The first three commands are short.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -B literal_channels.py | diff -u literal-channels.json -
python3 -B colour_channel_search.py controls --prime 1009
python3 -B colour_channel_search.py controls --prime 1013
python3 -B global_source_scan.py replay --prime 1013 --certificate jacobi-3088.json --workers 48
python3 -B colour_channel_search.py replay --prime 1013 --certificate channels-3060.json --workers 48
~~~

The last two commands rebuild the minor from the native colours. An uncompressed list of colour points is unnecessary because each recorded row deterministically specifies its integer assignment. The corrected quantum functional on the residual classical quotient remains to be determined. [Ranks and replays](results.json)


## Two virtual colour representations

Leave a 26-dimensional Albert multiplication operator $L$ and a 52-dimensional generator $T$ open in a gauge-invariant colour trace $\operatorname{Tr}_{26}(LAT)$, where $A$ is a matrix word. Append a physical colour $X$ by sending the adjoint argument to $[X,T]$ or the defining argument to its native $X$ action. Pair two such tensors through the inverse trace forms of both representations. The replay verifies the bracket representation on every native generator and the intertwining relation between the defining action and the Albert multiplication operators. It then applies the same third-jet and two-second-jet differentials as above. A separate calculation evaluates nonzero open entries by literal Taylor permutations, fresh native brackets and a single defining-representation trace. [Mixed colour maps](mixed_colour_edges.py) · [Literal entry controls](literal-mixed.json)

Use the previously certified 2879 evaluation rows and columns as an invertible square block. Write $A$ for that matrix, $C$ for the old columns on new points, and $B,D$ for the new columns at old and new points respectively. Solve $W A=C$ at the current prime and verify every entry. The new rank is measured on $D-WB$. A minor with the old rows and columns first has determinant $\det(A)\det((D-WB)_{I,J})$ for the selected extra rows $I$ and columns $J$. This preserves a single rational source evaluation throughout. On 168 new points, 123 extra boundary columns are certified independently of the first 2879. The determinant products are $431\cdot325=833\pmod{1009}$ and $100\cdot462=615\pmod{1013}$. In particular, the residual bound $3017-3002=15$ is over the rationals. The 168 evaluation rows in the complement permit a strict upper bound only from the previously known classical kernel dimension. [Selected extra columns and rows](mixed-minor.json)

~~~sh
python3 -B literal_mixed_edges.py | diff -u literal-mixed.json -
python3 -B mixed_source_search.py replay --workers 48 --prime 1013
~~~
