# The level-32 USp(4) source at flavour (6,5,5)

At $j=32$, the flavour partition $(6,5,5)$ admits both a second-order and a third-order target under filtered Taylor equivariance. Among the six single-additional-core explanations of the symplectic primary index in the archived search, it is the only one whose third flavour has at least five slots. A6 job 7658750 extended the ordinary and primary indices to levels 162 and 192 with first part at most 50. All six explanations still satisfy the index upper bounds. [Exact indices](../7658750/receipt.json) · [Six alternatives](../C2_level192_bounds.json) · [Product test](../generator_bounds.py)

For the classical supercharge $Q_0$, each nonzero minor modulo $2147483647$ bounds the rational rank from below. Since $Q_0^2=0$, the sum of the adjacent lower bounds is bounded by the independent rational cochain dimension. A6 job 7658772 attained equality in the following five odd letter lengths. [Complete rank files and receipt](../7658772/receipt.json) · [Classical engine](../g2_chain.cpp) · [Reproducible seeds and independent dimensions](../usp4_655_chain.py)

| Letters | Cochains | Incoming rank | Outgoing rank | Classical cohomology |
|---:|---:|---:|---:|---:|
| 7 | 84 | 10 | 74 | 0 |
| 9 | 3805 | 722 | 3083 | 0 |
| 11 | 13616 | 6240 | 7376 | 0 |
| 13 | 8431 | 5670 | 2761 | 0 |
| 15 | 739 | 687 | 52 | 0 |

For each even Taylor axis, filter the positive Taylor jets by assigning cost zero to every pure jet and the square of its total degree to each mixed jet. The classical differential preserves the pure-axis complexes. Every split of a mixed jet strictly lowers cost. The first page is a tensor product of the pure-axis cohomologies and the exterior powers of the adjoint for each remaining mixed jet. At any fixed weight and length the finite filtration gives an upper bound on classical cohomology. Computing the dual one-axis Chevalley–Eilenberg boundary over the integers gives these total colour dimensions, indexed by axis degree and letter length. [Integer bracket matrices](../symplectic_axis_bounds.py) · [Every dominant weight and adjacent rank](symplectic_axis_bounds.json)

| Axis degree | Length | Homology dimension |
|---:|---:|---:|
| 0 | 0 | 1 |
| 1 | 1 | 10 |
| 2 | 2 | 35 |
| 3 | 3 | 30 |
| 4 | 3 | 105 |
| 5 | 4 | 238 |
| 6 | all | 0 |

A weight dominated by the target Young diagram occurs in its Schur module. Computing the colour-singlet upper bound at one such weight bounds the multiplicity of the whole diagram. The second-order diagram is $(5,5,3,1)$ and the third-order diagram is $(6,6,1)$. The target length is the source length plus one minus the order. [Weyl singlet calculation and dominance checks](../symplectic_axis_colours.py) · [Numerical bounds](symplectic_axis_colours.json)

| Order | Source length | Target length | Weight used | Target cochains | Upper bound |
|---:|---:|---:|---|---:|---:|
| 2 | 13 | 12 | $(5,5,3,1)$ | 6948 | 92 |
| 3 | 11 | 9 | $(6,6,1)$ | 4569 | 554 |
| 3 | 13 | 11 | $(6,5,2)$ | 616 | 0 |
| 2 | 15 | 14 | $(4,4,4,2)$ | 57 | 0 |

At source length eleven the remaining third-order target has length nine. Restrict the Taylor functor to two even directions and one odd direction. The diagram $(6,6,1)$ occurs once in the component with even spin weight $(6,6)$ and odd flavour degree one. This follows by removing the bottom box in the odd direction, leaving $(6,6)$ on the even directions. The classical differential acts only on the multiplicity space and commutes with this restriction. The positive multiplicity makes vanishing of the whole containing sector sufficient. [Taylor branching](../../taylor_hooks/README.md)

A6 job 7658777 computes the rational upper bound in that containing sector,

$$
6741-4386-2355=0.
$$

The complete adjacent cochain dimensions at lengths $8,9,10$ are $9259,6741,3082$. Both rank lower bounds were obtained over the same prime $2147483647$. The positive-characteristic ranks saturate the characteristic-zero dimension and $Q_0^2=0$. This closes the third-order channel at length eleven, and hence at every odd length of the $USp(4)$ source $(6,5,5)$ under filtered Taylor equivariance. The five direct source equalities above also exclude every perturbative obstruction originating in this odd-length column, independently of filtered Taylor equivariance. [A6 rank record](../7658777/one_scalar.json) · [Remote result hash](../7658777/receipt.json)

From the repository root, the containing-sector calculation can be reproduced with the archived C++17 engine. Its printed $j$ is the level of the restricted Taylor weight.

```sh
d=algorithms/dual_singlets_20261118
c++ -O3 -std=c++17 -pthread "$d/taylor_hooks/native/engine.cpp" -o /tmp/usp-taylor-engine
A_NORHO=1 A_CHECK=1 PAR_THREADS=8 /tmp/usp-taylor-engine chain Sp2 1 0 0 6 6 9 9 831277 2147483647
```

## Second-order target on proper connected subgroups

For source length eleven, the second target has Taylor diagram $(5,5,3,1)$ and length ten. In two even and three odd Taylor directions, its component of spin partition $(2,2)$ and flavour partition $(4,3,3)$ occurs once. Indeed, adding the transpose $(3,3,3,1)$ of the odd partition to $(2,2)$ row by row gives the diagram itself, the multiplicity-one highest summand of that tensor product.

| Compact subgroup | Cochains in the containing component | Classical dimension at most | Kernel of restriction to its Cartan at most |
|---|---:|---:|---:|
| $U(2)$, even number of central letters | 24078 | 726 | 0 |
| $SU(2)\times SU(2)$ | 137907 | 694 | 0 |
| Principal $SU(2)$ | 1708 | 0 | 0 |

The normalizer of $U(2)$ reverses its centre and acts on its simple algebra by an inner automorphism. Only the even central degrees occur in a restricted symplectic invariant. Decomposing them gives 210 terms. For the two commuting simple factors the Künneth formula and flavour and spin branching give 417 unordered blocks. On a tensor product, if each factor injects into its Cartan cohomology, their tensor product does also. The calculation includes cases in which the other factor's cohomology vanishes. Direct Weyl integration on the six-dimensional adjoint reproduces 137907. The three rows together use 529 different nonvacuum rank-one blocks, 464 previously recorded and 65 from A6 job 7658782. [Character, ranks and Cartan test](../rank_two_second_wall.py) · [New minors and receipt](../7658782/receipt.json) · [Bounds](second_wall.json)

These subgroups suffice for compact connected proper subgroups. A reducible action on the defining complex symplectic four-space preserves a quaternionic line or a complex Lagrangian plane and is contained in $Sp(1)\times Sp(1)$ or $U(2)$. An irreducible action has semisimple Lie algebra since its connected centre would act by symplectic scalars. Two two-dimensional simple factors have a symmetric invariant form on their tensor product. The only proper simple possibility with an irreducible symplectic four is the cubic symmetric power of the defining $SU(2)$. Thus any classical Cartan-kernel class in this target restricts to a boundary on every proper connected compact subgroup. The five source ranks already settle the specified outgoing physical transitions.

From the repository root, using NumPy, SciPy and python-flint, run

```sh
d=algorithms/dual_singlets_20261118/bmn_primary
python3 "$d/rank_two_second_wall.py" --out /tmp/usp4-second-wall.json
cmp "$d/symplectic_655/second_wall.json" /tmp/usp4-second-wall.json
```
