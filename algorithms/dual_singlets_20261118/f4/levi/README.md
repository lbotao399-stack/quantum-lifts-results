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
