# The second differential at Taylor excess two

In the relative jet complex of undeformed four-dimensional $\mathcal N=4$ Yang–Mills with simple gauge algebra, let $D$ count all Taylor slots and let $Y$ count letters. Over $\mathbb Q$ the second supercharge differential vanishes on every source with $D-Y\le 2$ and $D<12$. At $D=12$ the only possible source diagram in this range is $(3,3,3,3)$. The result settles the previously modular $G_2$ level-24 test over the rationals. The fan kernels and marked product are those of Li, Liu and Chang. The reverse Taylor correction and the excess-one argument were obtained previously in the public project. [Fan kernels and quantum ring](https://arxiv.org/abs/2611.24160) · [First-order reverse Ward identity](https://github.com/lbotao399-stack/quantum-lifts-results/blob/acf631f4f92c2c271bbf1d6076d5d576e5dd0b7b/algorithms/su3_j30_classical_ranks/bmn_return/README.md#the-nine-letter-target-survives-every-perturbative-page) · [Previous excess-one argument](https://github.com/lbotao399-stack/quantum-lifts-results/blob/acf631f4f92c2c271bbf1d6076d5d576e5dd0b7b/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#a-single-taylor-representation-can-hit-the-sixth-glueball)

## Second Ward identity

On a word with $D-Y\le2$, three selected letters have at most five Taylor slots. The fan $Q_2$ consumes six flavour slots, hence vanishes. The reverse odd shear $\bar P_0$ replaces one dotted derivative by one flavour derivative. Its two-input correction $\bar P_1$ consumes one slot of each other flavour and adds a derivative in the other dotted direction. Both preserve $D-Y$. The public first-order equations on gauge invariants are

$$
\{Q_0,\bar P_0\}=0,\qquad
\{Q_1,\bar P_0\}=-\{Q_0,\bar P_1\}.
$$

In $\{Q_1,\bar P_1\}$ the disjoint pairs cancel by odd parity. A pair used twice has at most four slots. For pairs overlapping on one position the three initial letters must carry exactly five flavour slots. Fix the restored flavour to be $i$. Its multiplicity is one, and the other two multiplicities are two each. All three inputs then have zero dotted degree. The four possible unordered patterns and both ordered compositions are enumerated by the verification program.

The connected colour diagram is a tree with three structure constants and five external legs. Invariance of the metric moves the external root to any specified leg. Its value is then the pairing of that leg with a multilinear Lie polynomial in the other four. Expand each bracket as an associative commutator. To extract a cyclic coefficient in the trace of five independent matrix variables, assign its five labels successively to $E_{12},E_{23},E_{34},E_{45},E_{51}$ in $\mathfrak{gl}_5$ with its trace metric. Precisely that cyclic word has trace one. Fixing the first input edge by cyclic rotation, assigning the next two input edges in all $4\cdot3$ orders and retaining every output coordinate tests all cyclic coefficients. Transposing the three input units to pass to the trace-dual coordinate frame reverses the orientation of the test cycle and leaves the family of tests equivalent by relabeling the five vertices. Vanishing of these coefficients implies the rooted polynomial vanishes in the free associative algebra, hence in the free Lie algebra. Invariance gives the same identity with every supplied invariant metric.

The exact first-order solve uses all pairs through total degree five, including the possible predecessors for an incoming $d_1$. The one-leg colour loop is the adjoint Casimir. For the connected four-leg first-order equations Jacobi leaves at most two independent pairings. Their evaluations in the cross-product algebra at $(0,0,1,1)$ and $(0,1,0,1)$ give the two coordinate vectors. In the second-order test all six choices of dotted direction and restored flavour give the following counts per choice.

| Verification | Count |
|---|---:|
| Independent pivot coefficients in the first Ward solve | 184 |
| Nonempty first Ward equations in the three-colour replay | 10,719 |
| Forward-shear colour-pair tests through total degree five | 1,899 |
| Three-input tests on a directed five-cycle | 48 |
| Nonzero output entries in the two separate compositions, before addition | 1,216 |
| Residual entries after addition | 0 |

Consequently $\{Q_1,\bar P_1\}=0$ on this full excess range. The forward shear replaces a flavour slot by a dotted slot with its divided-jet coefficient. Its graded commutator with $Q_0$ vanishes by Taylor covariance. The four-edge evaluations verify its graded commutator with $Q_1$ directly. The connected forward colour symbol has four external legs. The same colour-word argument applies on four vertices. [Exact tests](verify.py) · [Rational output](generated.json)

For a leading cycle $O_0$ and a first correction $Q_0 O_1=-Q_1 O_0$, the reverse action has first correction $\bar P_0O_1+\bar P_1O_0$. Substitution gives

$$
Q_1(\bar P_0O_1+\bar P_1O_0)+\bar P_0 Q_1O_1
=-Q_0\bar P_1O_1+\{Q_1,\bar P_1\}O_0
=-Q_0\bar P_1O_1.
$$

Thus $d_2$ intertwines the linear Taylor algebra, with the odd sign and even weight shift $(2,2\mid-2,-2,-2)$. Classical Taylor degree $D$ is obtained by polarizing $D$ labelled slots. Its cochains are subquotients of sums of covariant tensor powers of the two-even, three-odd Taylor space. Super Schur–Weyl duality makes these tensor powers semisimple in characteristic zero. The first Ward equations make both the kernel and the incoming image of $d_1$ stable under the classical linear shear action, so this description passes to $E_2=H(H(Q_0),d_1)$ at each degree. [Taylor character convention and public calculation](https://github.com/lbotao399-stack/quantum-lifts-results/blob/acf631f4f92c2c271bbf1d6076d5d576e5dd0b7b/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#a-single-taylor-representation-can-hit-the-sixth-glueball)

A covariant diagram $\lambda$ has odd highest entries $\max(\lambda'_j-2,0)$ for column labels $j=1,2,3$, where $\lambda'_j$ is the column height. The twist by the square of the Berezinian subtracts two from these three entries. A polynomial target therefore requires $\lambda'_3\ge4$. Its source contains the rectangle $(3,3,3,3)$, of area twelve. This proves the asserted bound.

## Two perturbative endpoints in G₂

Use the defining seven-dimensional trace. In the following tuples the flavour degrees precede the two dotted degrees. At $(5,2,2;1,1)$, $Y=9$, the existing complete rational bounds give $2872-1825-1046=1$ for the tree cohomology. Both neighbouring first-differential spaces vanish. Here $D=11$ and $D-Y=2$, hence $d_2=0$ over $\mathbb Q$. For $r\ge3$ the outgoing $d_r$ would require a negative second flavour degree. For $r\ge2$ the incoming $d_r$ would require a negative dotted degree. This source is one-dimensional on perturbative $E_\infty$. [Rational source ranks and previous two-prime primitives](https://github.com/lbotao399-stack/quantum-lifts-results/tree/e965e8015f72466a8f26fb7993af53baa0bcd16d/reviews/independent-g2-and-r007/g2-filtered)

The physical outgoing charge block with doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(1,1,12,6,6)$ has dotted degrees $(t,t)$, flavour degrees $(6-t,3-t,3-t)$ and length $11-t$. Here $t$ ranges from zero through three. The three rank columns below are cochain dimension and adjacent tree ranks or rank lower bounds. A filled row proves vanishing over $\mathbb Q$.

| $t$ | $Y$ | Dimension | Incoming | Outgoing | $\dim H(Q_0)$ |
|---:|---:|---:|---:|---:|---:|
| 0 | 11 | 226 | 201 | 25 | 0 |
| 1 | 10 | 1222 | 1046 | 176 | 0 |
| 2 | 9 | 944 | 796 | 148 | 0 |
| 3 | 8 | 115 | 101 | 13 | 1 |

For the last row, at lengths seven and nine the two saturation equations are $252=151+101$ and $13=13+0$. They fix the two rational ranks at length eight. Adding the eight specified rational graviton products to the 101 boundaries gives rank 102 at two primes. The first recipe spans the quotient and all eight recipes vanish identically on the Cartan. [Incoming physical charge convention and level-24 products](https://github.com/lbotao399-stack/quantum-lifts-results/tree/d8fbf518ce4448ece8c844f5291e02a98e3a0f4f/reviews/independent-g2-and-r007/g2-level24-charged-products) · [Two-prime target replay and rational Cartan calculation](target.py)

The incoming $d_1$ source for this last line is spanned by marked graviton products. Its incoming $d_2$ source is the one-dimensional line just proved permanent. Its incoming $d_3$ source consists of the two marked products at $(6,3,3;0,0)$, $Y=10$. For $r\ge4$ the predecessor dotted degrees are negative. The target is itself a closed marked product. Thus the last row also has dimension one at perturbative $E_\infty$. An instanton lift through this outgoing channel still requires a nonzero physical matrix element between these two perturbative cohomology lines.

## Replay

Use Python 3.13. The Ward calculation uses the public reverse-shear operator and its ordered-pair runtime at commit `acf631f4f92c2c271bbf1d6076d5d576e5dd0b7b`. The target calculation uses the earlier native engine at `ee9e06494a926fb88baeb6fd875cee36762f63a7`, the published G₂ filtered certificate, a C++17 compiler, NumPy, SciPy, SymPy and python-flint. Supply the corresponding checkout directories as `WARD` and `ENGINE`, and the extracted G₂ filtered certificate directory as `FILTERED`.

~~~sh
python3 -B verify.py --repo "$WARD" | diff -u generated.json -
python3 -B target.py --repo "$ENGINE" --filtered "$FILTERED" --out rebuilt-target.json
cmp target.json rebuilt-target.json
~~~
