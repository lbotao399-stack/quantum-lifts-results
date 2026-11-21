# An odd BMN generator in $G_2$

The two remaining $(6,6,4)$ ranks vanish. The generator accounting for the level-36 primary has weight $(6,6,6)$, and a second odd generator is required by level 44. [New rank certificates and index inequalities](second_odd_generator.md)

At level $j=36$ and flavour partition $(6,6,6)$ the classical BMN cohomology contains an odd-length primary. Here a primary is a class modulo the images of the three BMN raising supercharges $Q_+^i$, and length $Y$ counts letters. At least one odd generator contributing to this class has the following level and flavour partition.

| $j$ | Flavour partition |
|---:|---|
| 30 | $(6,5,4)$ |
| 32 | $(6,6,4)$ |
| 36 | $(6,6,6)$ |

The generator count is relative to the BMN algebra of products and raising descendants. The two primitive scalar invariant degrees of $G_2$ are 2 and 6. Their raising descendants preserve length, and their products have even length. The calculation establishes a generator outside this algebra. A rational representative and its quantum survival remain open. [Rank and index verifier](verify.py) · [Product channels](product_channels.py)

## The primary index

On a BMN letter with $k=1,2,3$ antisymmetric flavour slots, $Q_+^i$ acts by exterior multiplication in direction $i$. The lowering charge acts by contraction with coefficient $k-1$. Their diagonal anticommutator counts slots in the other two directions. These operators act by graded derivations. The lowering charges anticommute with $Q_0$, and the anticommutator of a raising charge with $Q_0$ is an infinitesimal gauge transformation. Thus both act on gauge-invariant classical cohomology. At each fixed length there is a positive inner product for which the lowering charges are the adjoints of the raising charges. It is induced from the exterior flavour basis, with the relative norms fixed by the letter actions. The orthogonal complement of the boundaries inside the cycles is invariant under both actions. [BMN letters and supercharges](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/code/fortuity/bmn/words.py) · [Lowering action and letter checks](../su3_level36_return.py)

Let $v$ have flavour highest weight $(a,b,c)$, $a\ge b\ge c$, and be annihilated by all three lowering charges. In the order $Q_+^1$, $Q_+^2$, $Q_+^3$ the successive squared-norm multipliers are

$$
b+c,\qquad a+c+1,\qquad a+b+2.
$$

The off-diagonal anticommutators needed when commuting the adjoints to the right are upper-triangular flavour generators and annihilate $v$. For $b+c>0$ the threefold raising map is injective on the whole bottom flavour representation because it intertwines flavour up to a determinant. Any kernel in its exterior module would meet the top exterior degree by multiplication by the complementary raising charges. Consequently the full character of that multiplet is its bottom character times $(1+x_1)(1+x_2)(1+x_3)$, where $x_i$ records the number of slots of flavour $i$. Orthogonal complements give the decomposition into such multiplets and the short multiplets with one-row bottom. In a one-row multiplet the second row of any descendant has length at most one by the exterior Pieri rule.

To extract the primary index at $(6,6,6)$, divide by these three exterior factors. The inverse in degree $r$ has character $(-1)^r$ times the symmetric power of degree $r$. The horizontal-strip rule uses exactly the seven ordinary indices at $(6,6,k)$, $0\le k\le6$. One-row multiplets cannot enter these seven coefficients. With the index signed by $(-1)^Y$, the two evaluations give

$$
1-3+7-11+11-9+3=-1,\qquad 34961-34962=-1.
$$

There is also a check using only the ordinary index at level 36 and descendants from lower levels. The primary indices at $(5,5,5)$, $(6,5,5)$ and $(6,6,5)$ are $1,1,2$. The exterior Pieri rule supplies one copy of $(6,6,6)$ from each corresponding multiplet, at the same letter length. These supply at least four even classes. The ordinary index at $(6,6,6)$ is $+3$, so the odd dimension there is at least $4-3=1$. [A8 job 5690901](5690901/receipt.json) · [Integer character algorithm](../bmn_primary_index.py)

## Products and lower levels

The public $G_2$ matrices give lower bounds on the ranks entering and leaving each length. For **all 293** nonempty odd-length highest-weight blocks through level 28 their sum equals the independently computed integer cochain dimension. Since $Q_0^2=0$, the ranks exhaust those blocks over the rationals. [All 293 equalities](odd_through_28.json) · [Upstream rank data](https://github.com/davidsd/fortuity-ai/tree/4256f207986029540060fb44dd494746828f65ce/data/g2)

An odd product through level 36 must contain an odd factor of level at least 30. Every raising action on that factor can be moved onto its dressing modulo total descendants, by the graded Leibniz rule. The available dressing levels are $0,2,4,6$, with characters $1,0,s_{(2)},s_{(2,1)}$, respectively. Here $s_\lambda$ denotes the flavour character of partition $\lambda$. A product containing the determinant representation $(6,6,6)$ pairs each dressing with its unique complementary partition in the three-by-six rectangle. This gives precisely the three rows above. In particular level 34 has no compatible dressing. This excludes level 34 as the origin of this signal while allowing unrelated level-34 generators.

At levels 44 and 50 the primary indices at $(8,7,7)$ and $(9,8,8)$ are $-2$ and $-1$. One generator in each of the three candidate rows can supply at most the numbers below. Each entry is the multiplicity of the target in the product of the proposed bottom flavour representation and the entire allowed lower cohomology. All required dressing levels are at most 20, inside the odd-vanishing census. In particular, if the first two candidate rows supply no odd generator and there is a single generator in the level-36 row, at least one further odd generator is required by level 44. These higher weights add no unconditional generator to the count. [Character multiplication](product_channels.py)

| Target | Required odd primaries | From $(6,5,4)$ | From $(6,6,4)$ | From $(6,6,6)$ |
|---|---:|---:|---:|---:|
| $(6,6,6)$ | 1 | 1 | 1 | 1 |
| $(8,7,7)$ | 2 | 7 | 4 | 1 |
| $(9,8,8)$ | 1 | 16 | 8 | 2 |

## Higher-supercharge targets

For the level-36 class the odd source lengths allowed by its primary cochains are $9,11,13,15,17$. Under filtered Taylor equivariance, the second and third targets have Young diagrams $(5,5,3,3)$ and $(6,6,3)$, and lengths $Y-1$ and $Y-2$. The entries count target cochains before taking classical boundaries or lower supercharge images. Orders five and higher have a negative odd highest-weight entry and hence no polynomial target under this selection rule. [A8 job 5691072](5691072/receipt.json) · [Taylor algorithm](../g2_primary_targets.py) · [Selection hypothesis](../../perturbative_cartan_targets.md)

| Source length | $d_2$ target cochains | $d_3$ target cochains |
|---:|---:|---:|
| 9 | 67070 | 21151 |
| 11 | 303806 | 80487 |
| 13 | 262452 | 61625 |
| 15 | 38435 | 7768 |
| 17 | 594 | 104 |

## $SU(3)$ descendant control

At level 36, flavour weight $(8,5,5)$, the even charge-conjugation, eleven-letter classical space has dimension two. Use the ordered basis of descendants from levels 32 and 34. Applying the third lowering charge, raising again along the third axis and projecting onto the flavour highest weight gives the matrix

$$
\begin{pmatrix}35/3&221/108\\0&25/4\end{pmatrix},
\qquad \det=875/12.
$$

The return acts through level 34 and has rank two. This exhausts the entire even space with raising descendants. [Original and returned rational pairings](5690811/receipt.json) · [Trace evaluation](../su3_level36_return.py) · [Classical dimension](../taylor_hooks/su3_bmn_sources/source_lengths.json)

From the repository root, the index and rank audit runs with the Python standard library.

```sh
python3 algorithms/dual_singlets_20261118/bmn_primary/verify.py
python3 algorithms/dual_singlets_20261118/bmn_primary/product_channels.py
```

Passing `--upstream /path/to/fortuity-ai/data/g2` to the verifier also checks the six public Git blobs and regenerates the full odd-length table.

## A third odd generator in $USp(4)$

At weight $(8,8,8)$, level 48, the primary cochains of even and odd letter length give $878313-878315=-2$. Thus at least two odd primary classes occur. The two independent generators already proved at $(10,4,4)$ cannot produce either class. A polynomial Littlewood–Richardson product contains the diagram of each factor, and the raising operations add boxes. The first row of a descendant or dressing of either generator is therefore at least ten. At least **one further odd generator** is required. The two units of deficit by themselves allow the same new generator to contribute twice. [Integer cochains](5690901/symplectic.json) · [Two level-36 generators](../b2_level36/README.md)

If a single additional generator accounts for all deficits through level 108, its only possible weights are $(6,5,4),(6,5,5),(7,4,4),(7,5,3),(7,5,4),(8,5,2)$. To bound odd cohomology, multiply each proposed odd generator, including its raising descendants, by the entire even cohomology of the complementary degree. An odd product always has such a factor. The even dimension equals the ordinary index plus the odd dimension; recursion terminates because each generator has positive degree. For the primary quotient, the graded Leibniz rule moves raising operations off the odd factor and onto the dressing. Each step gives an upper bound, allowing relations and overlapping images. A candidate is excluded when that upper bound falls below the negative primary index. [Six surviving alternatives](rank_two_indices/C2_bound_results.json) · [Character calculation](generator_bounds.py) · [Ordinary indices](rank_two_indices/C2_ordinary_exact.json) · [Primary indices](rank_two_indices/C2_primary_exact.json)

A6 job 7658467 evaluated the torus constant terms for six further groups. At weight $(18,18,17)$, level 106, the $SO(7)$ and $USp(6)$ primary indices are $-6$ and $-59$. At $(3,3,3)$ the $F_4$ entry is $-1$. Only this last entry is negative in the tested $F_4$ range through level 72 with third part at least two. The integer reconstruction uses a prime product greater than twice the absolute coefficient bound. The integration grid exceeds the maximum Fourier exponent of the integrand in every colour direction. [All six tables and hashes](7658467/receipt.json) · [Torus summation](torus_index.cpp) · [Integer reconstruction](torus_index.py)

## Shorter $G_2$ source lists

The following rational cohomology dimensions are zero. Each row gives letter length, exact cochain dimension and the two independently evaluated adjacent ranks. They saturate the exact cochain dimension. A terminal target needs only enough evaluation functionals to witness a nonzero image minor, so the newer jobs omit construction of its trace basis. [Complete edge job](5691136/receipt.json) · [Lower-weight job](5691142/receipt.json) · [Nine-letter $(6,5,4)$](5691420/receipt.json) · [Nine-letter $(6,6,4)$](5691422/receipt.json)

| Flavour weight | $(Y,\dim C,\operatorname{rank}_{\rm in},\operatorname{rank}_{\rm out})$ |
|---|---|
| $(6,5,4)$ | $(7,159,9,150)$, $(9,5038,1109,3929)$, $(15,230,230,0)$ |
| $(6,6,4)$ | $(7,66,0,66)$, $(9,3555,622,2933)$, $(15,2445,2091,354)$ |
| $(6,6,6)$ | $(7,4,2,2)$, $(9,1271,109,1162)$, $(17,1946,1660,286)$ |

There is also an exclusion in the order-three target for the level-36 singlet. Use the ordinary-even Taylor-coordinate realization of the polynomial Schur functor used in the character calculation. Its classical generators are odd colour ghosts indexed by positive Taylor multi-indices. The differential splits a multi-index into two positive multi-indices with the Lie bracket on colours. At fixed total weight filter by the sum of the following costs over the letters. Pure first and pure second jets along each coordinate axis cost zero. Every other jet costs the square of its total Taylor degree. A split of a positive-cost jet strictly lowers cost. The associated graded differential replaces a pure second jet by the bracket of two first jets on its own axis. The axis complexes tensor together with the remaining generators, whose graded differential is zero. This is a finite filtration in each weight.

For one axis with no second jet, its cohomology is the exterior algebra of the 14 first-jet components modulo the quadratic bracket relations. Multiplying those relations by exterior degree two gives a matrix onto the whole degree-four space. Its rank is $\binom{14}{4}=1001$ and the determinant in pivot order is $350171$ modulo $1000003$. Multiplication gives vanishing at every higher exterior degree as well. [Rows and pivots](g2_exterior_degree4.json.gz) · [Integer Lie brackets and elimination](exterior_bound.py)

If an axis has weight four, at least one of its slots must belong to a letter of degree at least two. If it has weight five or greater, at least two slots must do so. Otherwise that axis has at least four pure first jets, no pure second jet, and zero cohomology in the tensor product. A letter of degree $k\ge2$ carries $k\le2(k-1)$ such slots. The weight $(6,5,4)$ occurs in the Schur diagram $(6,6,3)$ by dominance. It requires at least $2+2+1=5$ slots in higher jets, hence Taylor degree minus letter length at least three. At source length 15 the third target has Taylor degree 15 and length 13, with excess $15-13=2$. Its selected Schur component therefore vanishes in classical cohomology. Under filtered Taylor equivariance the remaining source lengths in the level-36 singlet are $11,13$ for $d_3$ and $11,13,15$ for $d_2$. The inference about the physical page uses the same Ward hypothesis as the target table above. [Schur character](../b2_schur_complex.py) · [Selected diagrams and bounds](g2_exterior_degree4.json.gz)

The $F_4$ index at level 18 is realized explicitly. The scalar generator, full-letter partner calculation and level-24 exceptional source are given in the [Albert-algebra calculation](f4_scalar/README.md).

## Two remaining G2 weights

At weight $(6,5,4)$ the table now exhausts every odd letter length. The rank lower bounds add to the independent character dimension at each length. Nilpotence of the classical supercharge $Q_0$ gives equality over the rationals.

| Length | Dimension | Incoming rank | Outgoing rank |
|---:|---:|---:|---:|
| 7 | 159 | 9 | 150 |
| 9 | 5038 | 1109 | 3929 |
| 11 | 15608 | 7487 | 8121 |
| 13 | 6777 | 5065 | 1712 |
| 15 | 230 | 230 | 0 |

Thus the odd generator responsible for the primary index $-1$ at $(6,6,6)$ can originate at $(6,6,4)$ with length 11 or 13, or at $(6,6,6)$ with length 11, 13 or 15. In the former case the level-four quadratic dressing into the level-36 primary adds two letters. For the latter weight, length 15 already has zero third-order target under filtered Taylor equivariance. [Length eleven](5691483/chain_6_5_4_11_11.json) · [Length thirteen](5691423/chain_6_5_4_13_13.json) · [Receipt and character replay](verify_triality_and_sources.py)

For $SO(8)$, A6 job 7658707 integrated the identity, transposition and three-cycle triality traces. Write these signed indices as $I_e,I_{(12)},I_{(123)}$ in this order. Multiplicities in the trivial, sign and two-dimensional representations of the permutation group $S_3$ are
$$
\frac{I_e+3I_{(12)}+2I_{(123)}}6,\qquad
\frac{I_e-3I_{(12)}+2I_{(123)}}6,\qquad
\frac{I_e-I_{(123)}}3.
$$
All three are nonnegative integers at all 3571 partitions of level at most 144 with first part at most 28 and third part at least two. At $(3,3,3)$ the traces $(2,0,-1)$ give $(0,0,1)$, the known triality doublet. Each full trace uses nine primes with integer reconstruction modulus greater than twice the absolute coefficient bound. An index of either sign permits cancellation between opposite letter lengths. [Integer traces](7658707/receipt.json) · [Folded torus integrand](torus_index.cpp) · [Triality doublet in Miao and Tang](https://arxiv.org/abs/2611.23173)

The four two-scalar exceptional BMN cores, the two-variable index calculation and their relation to the previously published Higgs-branch multiplicities are given in the [exceptional scalar calculation](exceptional_scalar_pairs.md).
