# Where to search for a second differential into the sixth glueball

In the full relative holomorphic jet complex at $F_4$, $j=36$, every classical source cycle for the incoming $d_2$ on $S^6$ has a classical primitive wherever a linear combination of its four polarized scalar colours is regular. Regular means that the adjoint action has rank $48$. The primitive also exists near a tuple containing a scalar whose semisimple part has centralizer with semisimple factors of type $A_1$, $A_2$, $C_2$, or $A_2+A_1$, including products of rank-one factors. The remaining proper Levi types are $B_3$ and $C_3$. Tuples with only nilpotent linear combinations and no regular combination also remain. This is a support statement for the *classical source* of the differential. The global corrected pairing into $S^6$ is open. [Physical polarization](physical_flag.json) · [Earlier pages and corrected pairing](README.md#one-row-decides-the-second-differential)

### A polynomial primitive on every regular chart

Use the ordinary four-variable polarization of Taylor type $(4,4,4,2)$, the transpose of the physical diagram $(4,4,3,3)$. The variables $X_0,\ldots,X_3$ are even first jets. An alternating second jet $Y_{ij}$ satisfies $Q_0Y_{ij}=[X_i,X_j]$. A source has twelve first jets and one second jet. Its restriction to four Cartan elements with any marked colour vanishes. The Cartan torus first projects the marked colour to the Cartan. The longest Weyl element then changes the sign of all thirteen factors. [Native Weyl action and scalar diagrams](source_support.json)

Choose $X_0$ regular. Its centralizer is abelian. Over the regular open subset of the Lie algebra the kernels of $\operatorname{ad}_{X_0}$ form a rank-four vector bundle. Three copies over the same base are irreducible, and their restriction to regular semisimple $X_0$ is dense. Consequently the Cartan vanishing just proved extends to any three elements of this centralizer, with arbitrary marked colour.

For $A=\operatorname{ad}_{X_0}$ choose a nonzero minor $\Delta=\det A_{I,J}$ of order $48$. Let $D$ apply $(A_{I,J})^{-1}$ to the rows $I$ and insert the result in the columns $J$. Then $ADA=A$ and $DAD=D$. For $i>0$ use the following polynomial path over the ring with $\Delta$ inverted and its degree-one homotopy $h_t$:

$$
X_0(t)=X_0,\qquad h_tX_0=0,\qquad
X_i(t)=X_i-(1-t)DA X_i,\qquad h_tX_i=D Y_{0i}.
$$

Set $h_tY_{ij}=0$, extend $h_t$ by the graded Leibniz rule with spectators evaluated along the path, and take

$$
Y_{ij}(t)=Y_{ij}-\int_t^1
\bigl([h_uX_i,X_j(u)]+[X_i(u),h_uX_j]\bigr)\,du.
$$

In particular $Y_{0i}(t)=Y_{0i}-(1-t)ADY_{0i}$. Direct substitution gives $Q_0Y_{ij}(t)=[X_i(t),X_j(t)]$ and $\partial_t O(t)=Q_0h_tO+h_tQ_0O$. At $t=0$ all four first jets commute in the centralizer of $X_0$. Hence for every classical source cycle $O$,

$$
O(0)=0,\qquad
O=Q_0\int_0^1 h_t O\,dt .
$$

The integral is coefficientwise rational. If the chosen Taylor direction occupies $m$ boxes, its marked pair contains it either once or zero times. The numbers of other scalar factors are respectively $13-m$ and $12-m$. Each costs at most one denominator $\Delta$. The corresponding marked factors cost respectively one and two. The integral therefore gives a **polynomial** primitive after multiplication by $\Delta^{14-m}$. The four displayed weights give powers $10,10,10,12$. A basis-dependent minor gives an identity on the colour-coordinate chart.

There is also a gauge-invariant version on the regular semisimple locus. Write $\det(t-A)=t^4p(t)$ and take $\Delta=p(0)$. The polynomial quotient $(p(0)-p(t))/t$, evaluated at $A$ and divided by $\Delta$, supplies $D$. The spectrum is symmetric under sign and hence this polynomial vanishes at zero. It inverts $A$ on its image. This choice makes the denominator-cleared primitive gauge invariant with the same exponent.

At a regular semisimple test element both the selected minor and $p(0)$ are $582820744123750809600000000$. At the principal nilpotent they are $1152$ and $0$ respectively. All sixteen regular representatives pass the rational generator identities, including the six independent commutators in their centralizers. [Exact homotopies](regular_source_homotopy.py) · [Minors, rows and columns](regular_source_homotopy.json)

### Vanishing on the reduced commuting variety

For any commuting quadruple, take the simultaneous semisimple and nilpotent parts in the complex Lie algebra. Their common semisimple centralizer is a Levi algebra. Its centre contains all the semisimple parts, and its positive nilradical in a suitable Borel contains all the nilpotent parts. The central torus forces the marked colour to lie in the Levi. Removing the alternating marked pair from $(4,4,4,2)$ leaves the two scalar diagrams $(4,4,3,1)$ and $(4,3,3,2)$. Each column can use at most as many centre vectors as the dimension of the centre. The resulting minimum numbers of positive-root entries are:

| Semisimple rank of the Levi | Minimum for $(4,4,3,1)$ | Minimum for $(4,3,3,2)$ | Largest possible positive-root height |
|---:|---:|---:|---:|
| 2 | 4 | 5 | 3 |
| 3 | 8 | 8 | 5 |
| 4 | 12 | 12 | 11 |

Each such entry contributes height at least one. The marked adjoint colour can compensate at most the height in the last column. Thus the invariant vanishes. In rank one the commuting nilpotents are multiples of a single element of $\mathfrak{sl}_2$. Replacing that common element by a family of semisimple elements proves vanishing by specialization. Rank zero is the Cartan calculation. This proves pointwise vanishing for every source cochain, including on the other components of the reduced commuting variety. Detecting a global source requires information about the equations or their relations beyond evaluation at its ordinary points. [All sixteen root calculations](source_support.py)

### Reduction to a smaller gauge algebra

Write $\mathfrak l$ for the centralizer of the semisimple part of $X_0$ and $\mathfrak m$ for its complementary root spaces. Restrict $X_0$ to the open subset of $\mathfrak l$ where $\operatorname{ad}_{X_0}$ is invertible on $\mathfrak m$. In the same homotopy take $D$ to be this inverse on $\mathfrak m$ and zero on $\mathfrak l$. Now $X_i(0)$ is the projection of $X_i$ to $\mathfrak l$. The endpoint of the formula for $Q_0Y_{ij}(t)$ retains the Levi commutator. The central torus makes the original invariant linear in the mark depend at the endpoint only on the Levi projection of $Y_{ij}(0)$.

A classical Levi primitive may also contain a third jet $Z_{ijk}$, with $Q_0Z_{ijk}=[X_i,Y_{jk}]-[X_j,Y_{ik}]+[X_k,Y_{ij}]$. Set $Z_{ijk}(t)=Z_{ijk}-\int_t^1 h_uQ_0Z_{ijk}\,du$. Since $Q_0^2=0$, applying $Q_0$ gives the same three brackets evaluated on the path. At zero, Levi projection commutes with brackets with the three remaining Levi scalars. Thus substituting the endpoint into a Levi primitive gives an actual primitive on this slice. Conjugation from the group times the slice to the Lie algebra is smooth there: its derivative covers $\mathfrak l$ by variation on the slice and $\mathfrak m$ by $[X_0,\mathfrak m]=\mathfrak m$. Flat base change detects local exactness. The fourteen native slice tests use singular anchors, noncommuting projected scalars and all four third-jet equations. [Generator and integral checks](regular_source_homotopy.json)

For a marked $A_2$ or $C_2$ factor with two additional central scalar colours, removing two horizontal strips leaves 27 possible ordinary diagrams. For each diagram the highest component among three even and two odd physical polarizations occurs with multiplicity one. A8 jobs 5705974 and 5706051 tested the corresponding blocks in both groups at primes 1000003 and 1000033. Let $\rho$ denote Cartan evaluation. In all 54 blocks,

$$
\dim C-\operatorname{rank}Q_{0,\mathrm{in}}-
\operatorname{rank}\begin{pmatrix}Q_{0,\mathrm{out}}\\ \rho\end{pmatrix}=0 .
$$

The native cochain dimension is an integer character calculation. Each nonzero modular minor is a rational rank lower bound. Since the incoming image lies in the stacked kernel, equality of these bounds forces injectivity of Cartan restriction on rational cohomology. The largest block has $2971-2674-297=0$. If the mark occurs on an abelian factor, only terminal invariant functions remain on the semisimple factors. Their Cartan restrictions are injective for types $A$ and $C$. For an $A_1$ factor with the mark, use the all-degree one-second-jet theorem. With the mark on the other factor, an unmarked $A_1$ supplies one symmetric row of even degree by restriction to its one-dimensional Cartan. This is already included among the two horizontal strips in the 27-diagram test. Tensoring these injections accounts for $A_1+A_2$ and the products of rank-one factors. [All-degree inputs for types $A$, $C$ and $A_1$](../../../su2_one_second_jet/README.md) · [Full subgroup pivot indices and Slurm receipts](levi_minors.json.gz)

For $B_3$ and $C_3$ the product of their longest Weyl element with the longest $F_4$ element fixes their roots and reverses the central direction. On the semisimple algebra any residual action fixing its Cartan pointwise is inner. With a semisimple mark, only an even number of central scalars can occur. With a central mark, an odd number of central scalars is required, leaving an odd number of semisimple first jets. In the classical commuting quotient for these defining representations, odd traces vanish by transposition and reordering. The odd orthogonal adjoint invariants are generated by traces: its orthogonal and special-orthogonal conjugations coincide. Thus the central-mark part is zero, including scheme-theoretic terminal invariants.

Of the five semisimple-mark diagrams left by this parity, $(4,4,2)$, $(4,4,4)$ and $(4,4,3,1)$ have zero Cartan kernel. In the last diagram A8 job 5706421 gives $2999-2465-534=0$ for $SO(7)$ and $2957-2431-526=0$ for $USp(6)$. The remaining shapes are $(4,4,2,2)$ and $(4,4,4,2)$. [Subgroup verifier](verify_levi_minors.py) · [Physical polarization sizes](levi_polarizations.json)

### A smaller matrix for genuinely finite-colour classes

Let $I$ be the kernel of specializing stable traces to the indicated finite gauge group. Give the terminal cochains homological degree zero and the one-mark cochains degree one. The short exact sequence of the three classical complexes gives

$$
H_1(C_{\rm stable})\longrightarrow H_1(C_{\rm finite})
\longrightarrow I_0/Q_0 I_1\longrightarrow H_0(C_{\rm stable}).
$$

The quotient of finite cohomology by the image of stable cohomology is precisely the kernel of the last map. To compute this quotient it suffices to build $Q_0$ on finite trace identities and evaluate the remaining terminal identities in the stable commuting ring. The following are exact matrix sizes, with rows before columns:

| Gauge group | Ordinary diagram | $Q_0$ on the full single-diagram source | $Q_0$ on the trace identities |
|---|---|---:|---:|
| $SO(7)$ | $(4,4,2,2)$ | $305\times1568$ | $126\times360$ |
| $USp(6)$ | $(4,4,2,2)$ | $296\times1538$ | $135\times390$ |
| $SO(7)$ | $(4,4,4,2)$ | $843\times5720$ | $820\times3658$ |
| $USp(6)$ | $(4,4,4,2)$ | $809\times5525$ | $854\times3853$ |

For the orthogonal group, fundamental-vector pairings decompose into symmetric-group representations with even rows, cut off at seven rows. For the symplectic group the diagrams are transposed and are cut off at six rows. The permutation sign from transposition interchanges the two within-pair adjoint projectors. For the stable character, count the perfect matchings fixed by each permutation. Applying the Taylor character to these three colour characters gives the table. Summing the single-diagram dimensions with their branching multiplicities reproduces all three adjacent dimensions of the physical native blocks used as controls. In particular the alternative $SO(7)$ polarization $(4,4,0;2,2)$ has $16347,7266,1435$ cochains in those three lengths. [Exact comparison](levi_relation_spaces.json) · [Representation calculation](levi_relation_spaces.py)

For the exceptional group itself the 703-row classical interpolation minor gives determinants $1005\bmod1009$ and $734\bmod1013$. Against its specified seven quotient coordinates, 750 four-column adjoint probes have evaluation rank 89 and zero residual rank at both primes. This is a modular result on those probes. On the purely nilpotent part, there are 30 Weyl orbits of four distinct mutually commuting positive-root generators up to Weyl conjugacy. At these quadruples the first transverse Taylor order allowed by colour weight is four, five or six on respectively 12, 14 and 4 orbits. Weight allowance gives a lower bound on where a nonzero coefficient could occur. [A8 tables and coordinates](four_slot_results.json) · [All root quadruples](nilpotent_flag_orders.json) · [Independent classical minor](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank)

### Replay

Run in this directory with Python, NumPy, SymPy and python-flint. The final command regenerates one complete pair of pivot lists and each determinant from native trace evaluations.

~~~sh
python3 regular_source_homotopy.py | diff -u regular_source_homotopy.json -
python3 source_support.py | diff -u source_support.json -
python3 verify_levi_minors.py
python3 levi_relation_spaces.py | diff -u levi_relation_spaces.json -
python3 nilpotent_flag_orders.py --out /tmp/f4-root-orders.json
cmp nilpotent_flag_orders.json /tmp/f4-root-orders.json
c++ -O3 -std=c++17 -pthread -I ../../../taylor_hooks/native ../../../taylor_hooks/su2_cartan_kernel.cpp -o /tmp/f4-levi-minor
python3 verify_levi_minors.py --engine /tmp/f4-levi-minor --case SO7-4-4-3-1-1000003.json --threads 8
~~~
