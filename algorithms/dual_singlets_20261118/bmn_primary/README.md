# An odd BMN generator in $G_2$

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
