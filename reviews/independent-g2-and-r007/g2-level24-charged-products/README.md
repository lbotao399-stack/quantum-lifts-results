# Two Cartan-invisible level-24 graviton classes in the defining seven of G₂

At derivative-free weight $(6,3,3)$ and letter number $Y=10$, both classical highest-weight classes of undeformed four-dimensional $G_2$ Yang–Mills vanish on the Cartan. Both survive every page of the perturbative relative jet complex. One is the restriction of the previously published $SO(7)$ polynomial $P_0$. The calculation determines the four classical blocks from which a first physical instanton differential could hit the two classes. [Public polynomial and its distinct symplectic specialization](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/taylor_hooks/bmn_633)

## Source and perturbative completion

The classical differential $Q_0$ increases $Y$ by one. At lengths nine through eleven the exact cochain dimensions, followed by lower bounds on the two adjacent ranks, are

| $Y$ | Dimension | Incoming | Outgoing |
|---:|---:|---:|---:|
| 9 | 733 | 307 | 426 |
| 10 | 629 | 426 | 201 |
| 11 | 226 | 201 | 25 |

The first and last rows saturate their dimensions. From $Q_0^2=0$ the two middle ranks are exact over $\mathbb Q$. The resulting cohomology dimension is $629-426-201=2$. The Cartan multiplicity is zero by the Weyl-character calculation for the rank-two dihedral reflection representation in `verify.py`.

In the defining seven, adjoining the 22 highest-weight graviton polynomials to the boundary columns raises the rank from 426 to 428. The selected columns are 3 and 5, counted from zero. Appending $P_0$ alone raises it from 426 to 427. The stored square minors give

| Prime | Joined order-428 determinant | Joined order-427 determinant for $P_0$ |
|---:|---:|---:|
| 2147483647 | 706905894 | 1075960324 |
| 1000033 | 240965 | 257117 |

Each displayed determinant is nonzero modulo its prime and hence certifies a rational lower bound. The rank theorem above supplies the rational upper bounds. The marked product completes both basis elements to closed perturbative classes. Every higher incoming differential would require a negative dotted Taylor degree. Their leading classes therefore survive every page and give two independent classes after inverting the perturbative parameter. [Li–Liu–Chang, Theorem 9.1 and quantum product](https://arxiv.org/abs/2611.24160)

## Matching the incoming physical charges

The target has doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(2,2,12,6,6)$. For an incoming supercharge matrix element the doubled charges are $(3,3,12,6,6)$. At equal dotted degrees $t$ the flavour weight must be $(6-t,3-t,3-t)$ and $Y=9-t$. Nonnegative degrees give exactly $0\le t\le3$. In the table the joined columns are classical boundaries together with rational graviton products.

| $t$ | Dimension | Incoming rank at least | Outgoing rank at least | Joined rank at least |
|---:|---:|---:|---:|---:|
| 0 | 733 | 307 | 426 | 307 |
| 1 | 3215 | 1331 | 1825 | 1390 |
| 2 | 2205 | 923 | 1282 | 923 |
| 3 | 235 | 72 | 151 | 84 |

In each row the joined and outgoing bounds sum to the cochain dimension. All joined columns lie in the kernel. Every possible incoming classical class is therefore a linear combination of graviton products. The marked completions also span the surviving perturbative space, by subtracting a linear combination of completed products at each leading grade. If the same generators have closed products for the physical supercharge with these leading terms, the first instanton differential has zero incoming rank at the two targets and zero outgoing rank on their span. The perturbative theorem supplies no microscopic instanton coefficient.

## BMN census

The public Fortuity calculation already reported graviton saturation through level 28. The replay combines its modular ranks for classical boundaries and rational graviton products with an independent characteristic-zero $G_2$ Weyl-character count. In each of the 617 nonvacuum highest-weight blocks the joined lower bound plus the outgoing lower bound equals the exact cochain dimension. Thus the graviton products span the rational cohomology in these blocks. Every odd-$Y$ block has rational cohomology zero by its two adjacent boundary minors. For identifying the dimension in an even-$Y$ block, those neighbouring odd-length equalities fix the incoming and outgoing ranks. The total multiplicity across levels 4 through 28 is 675. No outgoing $d_2$ can start in this census. [Fortuity source and original census](https://github.com/davidsd/fortuity-ai/tree/4256f207986029540060fb44dd494746828f65ce)

## Replay

`python3 verify.py` checks the six pinned public tables, computes all 617 exact character dimensions, computes the Cartan multiplicity and verifies the recorded rank bounds. The full regeneration additionally checks tree closure of every individual graviton factor, constructs the finite-colour matrices from the rational products and recomputes the selected rows, columns and determinants. Use Python 3.13, NumPy, SciPy, SymPy, python-flint and a C++17 compiler. Set `REPO` to a checkout of [ee9e064](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7).

~~~sh
python3 verify.py
python3 verify.py --rebuild "$REPO" --workers 2
~~~
