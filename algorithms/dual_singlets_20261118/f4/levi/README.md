# Restrictions of the two exceptional BMN sources

In the level-24 source for the second supercharge differential in $F_4$, every classical cocycle restricts to a boundary on every proper Levi subalgebra and on $\mathfrak{so}_9$. Its restriction is also exact after fixing any scalar linear combination with a nonzero semisimple part. The full $F_4$ source cohomology still has dimension between zero and three. [Source bound and supercharge selection](../README.md#two-more-classical-boundaries-at-level-24)

The flavour highest weight is $(m,m,m)$, with $m=4,5$ at levels 24 and 30. The letter number is $3m-1$. Each cochain contains one matter fermion. At level 24 the central Weyl transformation $-1$ kills the restriction of every invariant to eleven Cartan letters. At level 30 the two root equations on the Cartan coefficients have determinant $-3780$, so every closed fourteen-letter source restricts to zero. [Root equations](../README.md#the-cartan-restriction-of-the-level-30-third-order-source)

## Branching the classical complex

The classical supercharge $Q_0$ acts by commutators. On a direct sum of Lie algebras its complex factors, and its cohomology is the tensor product over the rationals. A polynomial flavour representation of highest weight $(i,j,k)$ can contribute to the rectangle only with the complementary weight $(m-k,m-j,m-i)$. This follows by dualizing the representation and tensoring with the $m$th power of the determinant. The multiplicity is one.

For one abelian spectator, $s$ scalar letters give the symmetric power of the flavour fundamental. Adding a single matter fermion tensors it with the exterior square. The admissible active weights and active lengths are

| Spectator | Active highest weight | Active length | Range |
|---|---|---:|---|
| $s$ scalars | $(m,m,m-s)$ | $3m-s-1$ | $0\le s\le m$ |
| $s$ scalars and one fermion | $(m,m-1,m-s-1)$ | $3m-s-2$ | $0\le s<m$ |
| Second summand with one fermion | $(m-1,m-1,m-s)$ | $3m-s-2$ | $1\le s\le m$ |

Thus precisely $3m+1$ active blocks occur. The independent abelian character gives the same 13 and 16 blocks. For two simple factors the replay uses the three-row Weyl character formula and interlacing two-row weights. In each simple block, the cochain dimension minus the incoming rank lower bound and the rank lower bound of the stacked outgoing and Cartan maps bounds the rational kernel of Cartan restriction. The stored determinants are nonzero modulo 2147483647. The abelian complex restricts identically. [Branching and replay](../../f4_bmn_subgroups.py) · [Simple factors](../../f4_bmn_root_levi.py)

| Level | Subalgebra | Cochains | Classical dimension at most | Cartan kernel |
|---:|---|---:|---:|---:|
| 24 | $B_3+\mathfrak t^1$ | 779 | 4 | 0 |
| 24 | $C_3+\mathfrak t^1$ | 774 | 4 | 0 |
| 24 | $A_2+A_1+\mathfrak t^1$ | 1012 | 18 | 0 |
| 24 | $A_2+A_2$ | 818 | 6 | 0 |
| 24 | $C_3+A_1$ | 633 | 0 | 0 |
| 30 | $A_2+A_1+\mathfrak t^1$ | 5193 | 46 | 0 |
| 30 | $A_2+A_2$ | 5288 | 14 | 0 |

Here $\mathfrak t^1$ is the one-dimensional centre. Deleting one simple root of $F_4$ gives the first two rank-three types and the two embeddings of the third type. Restriction to a smaller Levi factors through a maximal Levi. The Cartan vanishing above and the zero kernels in the table make the indicated restrictions boundaries. [All component records](subgroup_bounds.json)

## The full-rank orthogonal subgroup

At $(4,4,4)$, the independent orthogonal character gives $559,341,121$ cochains at lengths $10,11,12$. The incoming and outgoing minors have orders 222 and 119. Their residues modulo 2147483647 are 1328936368 and 778753830. A second calculation modulo 1000033 gives both orders again. Hence

$$
\dim H(Q_0)_{\mathfrak{so}_9,(4,4,4),11}=341-222-119=0
$$

over the rationals, since consecutive maps compose to zero. The defining nine-dimensional matrices satisfy all $9^4$ completeness checks in the native engine. The computation supplies exact characteristic-zero colour multiplicities for all 28 letter patterns. [First minors](Spin9_2147483647_minors.json) · [Second minors](Spin9_1000033_minors.json) · [Character and reconstruction](../../f4_spin9_kernel.py)

## Fixing a scalar

Choose the first scalar $B$ to be a linear combination with nonzero semisimple part. Its Jordan decomposition is $B=B_s+B_n$. In the complexified Lie algebra take the zero eigenspace of $\operatorname{ad}_{B_s}$ as the Levi subalgebra and the sum of the nonzero eigenspaces as its complement. On the complement $\operatorname{ad}_B$ is invertible since $\operatorname{ad}_{B_n}$ is nilpotent and commutes with $\operatorname{ad}_{B_s}$.

At fixed $B$, the linear part of $Q_0$ pairs the two other scalars with the corresponding two matter fermions, and the remaining matter fermion with the field strength, via the signed map $\operatorname{ad}_B$. On each transverse pair use the signed inverse map in the reverse direction. The anticommutator counts the transverse entries. Dividing by their number contracts the ideal of positive transverse degree for the linear differential. The linear differential on the Levi entries anticommutes with this contraction since it acts on a disjoint set of generators.

Extend a Levi primitive by projecting each remaining letter to the Levi. Its closed residual belongs to the transverse ideal. The quadratic part of $Q_0$ preserves the ideal and raises the number of unfixed letters by one. Cancel the lowest such degree by the linear contraction and iterate. Every unfixed BMN letter carries at least one of the two other flavour degrees, whose sum is $2m$, so the procedure terminates. This proves exactness at every specified $B$ at level 24. It bounds the possible support at commuting triples to those for which every scalar pencil is nilpotent. Exactness on the entire unreduced $F_4$ complex and a quantum supercharge pairing require further calculation. [Signed classical rules and transverse contraction](../../taylor_hooks/su3_scalar_levi/README.md#the-full-two-colour-block-and-the-nilpotent-scalar-locus)

## The centre at level 30

At flavour weight $(5,5,5)$ and length 14, restriction to either rank-three simple Levi factor with its one-dimensional centre splits into 16 summands. The number of central letters is 14 minus the active length. The eight summands with odd central degree vanish on any restricted $F_4$ invariant. In the Euclidean root realization the reflection normal to $(0,0,0,1)$ fixes the $B_3$ root subsystem pointwise. The reflection normal to $(0,0,1,1)$ fixes $C_3$ pointwise. Each sends the centre to its negative. On the derived Lie algebra its representative fixes the Cartan pointwise and scales each root vector. Multiplicativity on root brackets makes this an inner torus action, which fixes the invariant cochains.

For the seven positive even summands the rational kernel of restriction to the Cartan has the following bounds. The stacked rank is a lower bound for the rank of the outgoing classical map together with Cartan evaluation.

| Derived group | Cochains | Incoming rank at least | Outgoing rank at least | Stacked rank at least | Cartan kernel at most |
|---|---:|---:|---:|---:|---:|
| $SO(7)$ | 2226 | 1833 | 380 | 393 | $2226-1833-393=0$ |
| $USp(6)$ | 2208 | 1814 | 381 | 394 | $2208-1814-394=0$ |

The degree-zero summand is the pure $(5,5,5)$ block of length 14 in the derived algebra. Every positive-central-degree coefficient of a closed $F_4$ source is a boundary because its Cartan restriction vanishes by the two root equations. The largest positive even symplectic summand has $1653-1322-331=0$ Cartan kernel. A8 jobs 5688192 and 5688957 reached their time limits after saving thirteen of the fourteen sector records used here. The final symplectic record was computed locally. [Individual ranks](level30_central_bounds.json) · [Minor rows and columns](level30_central_minors.json.gz) · [Provenance](level30_central_receipt.json) · [Cartan equations](../root_restriction.json)

For two or three central directions, Cauchy's formula gives the multiplicity of a central scalar partition as the dimension of that Schur representation of the central directions. One central matter fermion adds a vertical strip of two boxes, with one copy per central direction. Complementation in the $(5,5,5)$ rectangle gives the active weight. Successive branching over one central direction at a time reproduces every multiplicity. The simple factors have zero Cartan kernel in every participating summand. Tensoring their injections over the rationals gives the following bounds on the entire Levi source.

| Semisimple part and centre | Cochains | Classical dimension at most | Cartan kernel |
|---|---:|---:|---:|
| $B_2+\mathfrak t^2$ | 6528 | 72 | 0 |
| $A_2+\mathfrak t^2$ | 5258 | 91 | 0 |
| $A_1+A_1+\mathfrak t^2$ | 2850 | 130 | 0 |
| $A_1+\mathfrak t^3$ | 1708 | 251 | 0 |

Thus the full level-30 source restricts to a boundary on every proper Levi with derived rank at most two. The central calculation uses only classical restrictions. [Branching and replay](../../f4_level30_central.py)

## Replay

Use the Python environment for the [rank-four calculation](../../rank_four/README.md). From the repository root run

```sh
d=algorithms/dual_singlets_20261118
p="$d/f4/levi"
python3 "$d/f4_bmn_subgroups.py" --SU2 "$p/SU2.json.gz" --SU3 "$p/SU3.json.gz" --SO7 "$p/SO7.json.gz" --Sp3 "$p/Sp3.json.gz" --out /tmp/f4-subgroups.json
cmp /tmp/f4-subgroups.json "$p/subgroup_bounds.json"
python3 "$d/f4_level30_central.py" --out /tmp/f4-central-30.json
cmp /tmp/f4-central-30.json "$p/level30_central_bounds.json"
for prime in 2147483647 1000033; do
  python3 "$d/f4_spin9_kernel.py" --width 4 --prime "$prime" --work "/tmp/spin9-$prime"
  cmp "/tmp/spin9-$prime/minors.json" "$p/Spin9_${prime}_minors.json" || exit 1
done
```

## The remaining rank-three blocks at level 30

For the pure flavour weight $(5,5,5)$ let the stacked map consist of the outgoing classical supercharge $Q_0$ and restriction to the Cartan. Boundaries are in its kernel. In both derived algebras the outgoing minor fills its codomain. The stacked minor gains one rank, and the independent incoming minor exhausts the remaining kernel. Thus the full rational classical cohomology is one-dimensional and the Cartan detects it.

| Group | Cochains at lengths 13, 14, 15 | Incoming minor | Outgoing minor | Stacked minor | Cartan kernel dimension |
|---|---|---:|---:|---:|---:|
| $SO(7)$ | 6650, 2323, 418 | 1904 | 418 | 419 | $2323-1904-419=0$ |
| $USp(6)$ | 6499, 2276, 415 | 1860 | 415 | 416 | $2276-1860-416=0$ |

At the prime 2147483647 the three determinant residues are 1471201234, 1984488694, 413726119 in the orthogonal block and 786583585, 1816544542, 648149037 in the symplectic block. The characteristic-zero colour character independently gives all six cochain dimensions. [Orthogonal minors](SO7_level30_pure.json.gz) · [Symplectic minors](Sp3_level30_pure.json.gz) · [Dimensions and branches](level30_pure_bounds.json)

The two root equations in $F_4$ annihilate its full Cartan restriction. Combining these pure blocks with the positive central degrees above shows that the level-30 source restricts to a boundary on every proper Levi. The transverse contraction in the fixed-scalar calculation therefore applies also at level 30 whenever the fixed scalar has a nonzero semisimple part. The unrestricted exceptional source and its third differential remain open.

For the full-rank subgroup with Lie algebra $C_3+A_1$, the rectangular flavour representation branches with multiplicity one. Over the rationals the cohomology of each product is the tensor product of the two cohomologies. The $A_1$ factor has zero cohomology in 103 branches. Its other seven factors each have a Cartan-detected space of dimension at most one. The corresponding $C_3$ factors have Cartan-kernel dimension zero and classical dimensions bounded by $2,1,1,0,5,3,1$. Therefore the product has classical dimension at most 13 and Cartan-kernel dimension zero. Every restricted level-30 exceptional source is a boundary on this subgroup as well. [Factor calculation](../../f4_level30_pure.py)

From the repository root, set up the Python path as for the rank-four calculations. The native records can be regenerated with the pinned headers.

```sh
d=algorithms/dual_singlets_20261118
p="$d/f4/levi"
python3 "$d/f4_level30_pure.py" --out /tmp/f4-level30-pure.json
cmp /tmp/f4-level30-pure.json "$p/level30_pure_bounds.json"
c++ -O3 -std=c++17 -pthread -I "$d/taylor_hooks/native" "$d/taylor_hooks/su2_cartan_kernel.cpp" -o /tmp/f4-levi-kernel
for group in SO7 Sp3; do
  PAR_THREADS=3 /tmp/f4-levi-kernel 5 5 5 0 0 14 2147483647 2527911005270 "$group" > "/tmp/$group-level30.json" || exit 1
  python3 -c 'import gzip,json,sys; assert json.load(open(sys.argv[1])) == json.load(gzip.open(sys.argv[2]))' "/tmp/$group-level30.json" "$p/${group}_level30_pure.json.gz" || exit 1
done
```

## The $G_2+A_1$ restriction and parabolics

Let $A_1$ denote a rank-one simple factor, with Lie algebra $\mathfrak{sl}_2$. In a commuting product with another Lie algebra $\mathfrak h$, replacing this factor by its Cartan $\mathfrak t_1$ gives a subalgebra of the centralizer of $\mathfrak t_1$ in $F_4$. This centralizer is a proper Levi. Both exceptional sources already become boundaries there.

On the product the classical differential is the sum of the two differentials, with the exterior sign. Over the rationals its cohomology is the tensor product. In each complementary pair of flavour weights, restriction on this tensor product is the identity on the $\mathfrak h$ factor tensored with Cartan restriction on $A_1$. The rank-one minors give the following exact bounds, reproduced with primes 2147483647 and 1000033.

| Level | Product branches | Rank-one branches with zero cohomology | Other rank-one branches | Maximum Cartan-kernel dimension on any branch | Cochains and cohomology for $A_1$ alone |
|---:|---:|---:|---:|---:|---:|
| 24 | 68 | 62 | 6 | 0 | 1, 0 |
| 30 | 110 | 103 | 7 | 0 | 1, 0 |

Thus restriction from the product to $\mathfrak h+\mathfrak t_1$ is injective on cohomology in these source sectors. Its composition with restriction from $F_4$ vanishes. Every source therefore restricts to a boundary on $G_2+A_1$ at both levels. The last column also covers every simple $A_1$ embedding. [Rank-one summands](rank_one_factor.json) · [Native rank and character replay](../../f4_rank_one_factor.py) · [Draper, Albert-algebra realization of the subgroup](https://agt2.cie.uma.es/~loos/jordan/archive/F4/F4.pdf)

For a proper parabolic subgroup choose a cocharacter that acts trivially on its Levi and with positive weights on its nilradical. On a tuple of parabolic letters, conjugation by this cocharacter multiplies each nilradical entry by a positive power of its parameter. An $F_4$-invariant polynomial is constant under this conjugation. Its value at parameter zero is its restriction to the Levi pulled back by the parabolic projection. This projection preserves brackets, hence commutes with $Q_0$, and pulls back a Levi-invariant primitive to a parabolic-invariant primitive.

In characteristic zero the reductive maximal connected subgroups of $F_4$ have types $B_4$, $A_1+C_3$, $A_2+A_2$, $A_1+G_2$ and $A_1$. The remaining maximal connected subgroups are parabolic. Combining the restrictions above proves that the level-24 source restricts to a classical boundary on every proper connected algebraic subgroup over $\mathbb C$. This leaves the unrestricted $F_4$ source dimension between zero and three. [Thomas, Theorem 3.1, $F_4$ table](https://arxiv.org/pdf/1501.04858)

From the repository root, the second prime regenerates all 108 distinct positive-length rank-one sectors used in the two decompositions.

```sh
d=algorithms/dual_singlets_20261118
c++ -O3 -std=c++17 -pthread -I "$d/taylor_hooks/native" "$d/taylor_hooks/su2_cartan_kernel.cpp" -o /tmp/f4-rank-one
python3 "$d/f4_rank_one_factor.py" --engine /tmp/f4-rank-one --work /tmp/f4-rank-one-minors --out /tmp/f4-rank-one.json
cmp /tmp/f4-rank-one.json "$d/f4/levi/rank_one_factor.json"
```
