# Operator mixing at the charges of the fifth F4 glueball

The glueball $S^5$ has $E=15$, $J_L=J_R=0$ and $R=(5,5,5)$. Work in its dotted and flavour singlet. There are exactly five independent Cartan-visible **classical target** directions at these charges. The one-second-jet source at $E=29/2$, $J_L=1/2$, $J_R=0$, $R=(9/2,9/2,9/2)$ has 19 directions in the free Weyl-invariant Cartan algebra. The two numbers concern different ends of a possible supersymmetry-current matrix element. The fifth glueball itself restricts to zero on the rank-four Cartan. [Earlier perturbative survival proof](README.md#the-fifth-glueball-survives-to-perturbative-infinity)

| Dotted jet pairs | First scalar jets of each flavour in the target | Cartan target dimension | Cartan source dimensions by second-jet type |
|---:|---:|---:|---:|
| 0 | 5 | 0 | 2, 0, 0 |
| 1 | 4 | 3 | 0, 0, 0 |
| 2 | 3 | 0 | 7, 6, 1 |
| 3 | 2 | 2 | 0, 0, 0 |
| 4 | 1 | 0 | 1, 1, 1 |
| 5 | 0 | 0 | 0, 0, 0 |

For $k$ dotted pairs set $p=5-k$, the number of first scalars of each flavour in the target. Its colour characters on the Cartan $\mathfrak h$ are the Schur characters of shapes $(p,p,p)$ and $(k,k)^{\mathsf T}$; transposition accounts for the odd dotted letters. For the marked source the three types are a two-scalar mark, a mixed mark and a two-dotted mark. Removing their boxes gives the respective scalar and dotted shapes before dotted transposition:

$$
((p,p-1,p-1),(k,k)),\qquad
((p,p,p-1),(k,k-1)),\qquad
((p,p,p),(k,k-2)).
$$

Omit a shape with negative entries. The mark contributes one additional copy of $\mathfrak h$. Averaging these characters over the 1152 Weyl elements gives the table. A separate fixed-weight calculation with spin and flavour subtraction agrees on each target and on the six source totals. [Native Weyl characters and the earlier source census](../../../f4/cartan.json) · [Both computations](fifth_cartan_mixing.py)

### Native representatives for all five Cartan target directions

Use the existing form $B=\operatorname{Tr}_{26}/6$. For a triple of scalar test vectors $X$ in the Cartan, let $G$ be their $3\times3$ Gram matrix for $B$. For a single dotted test vector $T$ write $v$ for the three entries $B(X_i,T)$, where $i$ enumerates the three scalar vectors. Two dotted inputs can be set to the same $T$ after extracting their spin singlet, since their two colour slots are symmetric. The three evaluated polynomials at one dotted pair are

$$
(\det G)^2 B(T,T),\qquad
\det G\,v^{\mathsf T}\operatorname{adj}(G)v,\qquad
\operatorname{Tr}_{26}\!\left(T^2
       \big(X^{\mathsf T}\operatorname{adj}(G)X\big)^2\right).
$$

In the last expression the matrices commute on the Cartan. Its native extension uses the symmetrized six-input defining trace on the two dotted and four scalar entries before contracting with the two adjugates. Thus it is a polynomial gauge invariant, including at singular $G$. Independent polarizations recover the physical fermionic cochains. For three dotted pairs use a second triple $Z$ of test vectors for the Schur shape $(2,2,2)$. The two polynomials are the product of the two Gram determinants and the square of the determinant of the mixed Gram matrix. All five native polynomials use only first jets and are $Q_0$-closed. All tree boundaries vanish on Cartan evaluation. The two integer minors, on the explicit test points in the certificate, are

$$
43133924351672320000,\qquad -667762688 .
$$

Their ranks $3$ and $2$ exhaust the Cartan target dimensions. They supply subtraction directions for a matched current calculation. The physical nonperturbative matrix element into the Cartan-invisible part is open. [Points, exact matrices and determinants](fifth_cartan_mixing.json)

### An entire target block vanishes before instantons

At four dotted pairs and one scalar per flavour there are 144 invariant terminal cochains. The four Taylor multiplicities are $13,40,55,36$ for diagrams $(7,4),(6,4,1),(5,4,2),(4,4,3)$. The separate same-species bracket constraints leave zero in the first two. The shortened Ward system eliminates the third. In the last one the ordered harmonic multiplicity is 22, exchange has trace 12 and the simultaneous symmetry leaves $(22+12)/2=17$ dual coordinates. The mixed classical brackets have rank 17. The selected rational determinant is

$$
-\frac{2833763966967988198259375}{1174865639966244864}.
$$

The full classical target in this block therefore has dimension zero. Local replay and A8 job 5705731 agree on the matrix, all pivots and the determinant. The pure-glueball block has five dotted pairs. [Shortened Ward equations](glueball_factors.md) · [All 17 rows](fifth_adjacent_ward.json)

~~~sh
python3 fifth_cartan_mixing.py | diff -u fifth_cartan_mixing.json -
python3 fifth_adjacent_ward.py --out /tmp/f4-fifth-ward.json
cmp fifth_adjacent_ward.json /tmp/f4-fifth-ward.json
~~~
