# The crossed Albert candidate and the first page of the F₄ glueball tower

At level $j=36$ in the relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, the two proposed order-704 determinants for the crossed eight-Albert source B both vanish when the gauginos are put in the order used to define the graph coefficients. The [subsequent global calculation](../focused_powers/sixth_survival.md) proves perturbative survival of $[S^6]$. Here $S$ is the glueball in the native convention $B=\operatorname{Tr}_{26}/6$, and $d_r$ is the order-$r$ supercharge differential. [Dot’s exact input archive](https://github.com/lbotao399-stack/quantum-lifts-results/blob/783cdacb1171fb21ef5bdaa155fe14202a8d026d/reviews/f4-j36-eight-albert-d2/f4-j36-eight-Albert-certificate.zip) · [Exchange with Dot](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1#issuecomment-7105180683)

Dot independently recovered the ordered values and withdrew this certificate and its seventh-power corollary. [Withdrawal and separate replay](https://github.com/lbotao399-stack/quantum-lifts-results/blob/7fcfb7175703dd1ee6dcf4a71f2885c0d7dc8c73/reviews/f4-j36-eight-albert-d2/README.md)

| Pure power | Level | Verified result in the relative complex |
| --- | ---: | --- |
| $S^5$ | 30 | Survives to perturbative $E_\infty$. |
| $S^6$ | 36 | Nonzero on $E_2$. The crossed-graph minors vanish at both stated primes. |
| $S^7$ | 42 | Nonzero on $E_2$ by a separate rational Ward certificate. An incoming $d_2$ is open. |
| $S^8$ | 48 | Nonzero classically and hit by $d_1$. The incoming image on the one-dimensional glueball target has rank one. |

## Two orderings of the twelve gauginos

The two incidence determinants refer to the ordered edge lists `cp` and `cm`, six for each dotted species. In the supplied evaluator a leaf creates exterior variables in increasing incident-edge order. An internal node takes the exterior product of its first and second child for each species. Consequently the actual order is the sequence of owned edges obtained by visiting the leaves from left to right. Converting the contraction back to the declared edge order multiplies its value by the product of the two permutation signs. This follows directly by interchanging adjacent odd insertions of a given species. The graph weights in the supplied bank already use the declared edge lists. [Source packet and code](https://github.com/lbotao399-stack/quantum-lifts-results/tree/783cdacb1171fb21ef5bdaa155fe14202a8d026d/reviews/f4-j36-eight-albert-d2)

For left-marked graph 0 these lists are `[8,10,13,14,16,17]` and `[10,13,14,15,16,17]`. The evaluator visits the union in the order `[15,8,10,13,16,14,17]`. Restricting this sequence to each species gives one and four inversions respectively. Swapping the two children at the root’s right child transposes degrees $(0,1)$ and $(6,5)$. Its sign is $(-1)^{0\cdot6+1\cdot5}=-1$. The two exact extremal evaluations before the order conversion are $(-5750,24200)$ on the archived tree and $(5750,-24200)$ on the swapped tree. Both give $(5750,-24200)$ in the declared edge order.

| Prime | First mark | Archived score | Ordered score | Score of the 703-pivot section | Ordered difference | Archived determinant | Ordered determinant |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1009 | Left | 663 | 44 | 44 | 0 | 551 | 0 |
| 1009 | Right | 663 | 44 | 44 | 0 | 551 | 0 |
| 1013 | Left | 457 | 92 | 92 | 0 | 478 | 0 |

The two marks also give the same ordered graph sum over $\mathbb Q$. In both cases the sum before the packet’s normalization is $-92613596121446400$, hence the **unsubtracted** score is

$$
\frac{-92613596121446400}{2(28/3)^4}
=-\frac{299015516814300}{49}.
$$

Twelve disjoint-fan graphs are identical on the two marks. The left bank has 18 graphs and the right bank has 30, giving 36 distinct contractions with their declared edge occupations. The rational graph sum reduces to $44$ and $92$. For the section, the replay verifies the row identity on all 703 evaluation points and builds the last column from the four trace scores and the nested score in the packet’s normalization. The integer-only replay independently matches the eight-cubic outgoing row at all 703 points at both primes. Independent FLINT determinants for the classical block are $1005\bmod1009$ and $734\bmod1013$. Appending the physical column gives the last two columns of the table. The nonzero determinants of the original classical block and the rank-703 upper bound justify subtracting its section over $\mathbb Q$. The subsequent integer reconstruction gives the same section value $-299015516814300/49$ over $\mathbb Q$, so the corrected score for this source is exactly zero. The [global sixth-power proof](../focused_powers/sixth_survival.md) covers the remaining sources. [Exact rational subtraction on 720 points](../focused_powers/rational_section/README.md) [Classical rank](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank) · [Independent determinant replay](verify_crossed.py) · [Exact poles for every graph](graph_poles.json)

## The seventh and eighth powers

The zero-flavour target is strictly closed. The classical invariant quotient is $\mathbb C[S]/(S^9)$. For each power $p=7,8$, let $\lambda$ denote a covector annihilating the incoming classical image, normalized by the relatively prime native branch weights in the table. On $E_2=H(E_1,d_1)$ the class $[S^p]$ survives precisely when $\lambda Q_1$ annihilates the incoming $\ker Q_0$. By finite-dimensional duality this is the linear equation $\lambda Q_1=\eta Q_0$ for an intermediate covector $\eta$. [Classical chiral-ring theorem](https://doi.org/10.1002/prop.202400027) · [Native Ward equations and conventions](../first_page/README.md#the-complete-first-page-test)

| Power | Weights of $\lambda$ | Dimensions of the corresponding exterior representations | $\lambda(S^p)$ |
| ---: | --- | --- | ---: |
| 7 | $(874,1105,1035)$ | $(1264120,15997696,19214624)$ | $-194893931232000$ |
| 8 | $(3478750,6020625,5173597)$ | $(15611882,13530946,65609375)$ | $19160509614246000000$ |

The same-dot classical bracket quotient in each exterior degree is the sum of the abelian-ideal modules. The program grows the ideals in the native positive-root system, tests the mixed bracket on each pair of extremal vectors and obtains the unique line of branch weights. The value on $S^p$ is $(-1)^{p(p-1)/2}p!$ times the sum of the weighted dimensions above.

For the matter-fermion and scalar input, the particular solution uses the native symmetric scalar intertwiner. The 2,704 ordered colour pairs and 1,378 symmetric pairs satisfy its identities. The free part of $\eta$ is an invariant pairing of two degree-$(p-1)$ ghost harmonics with the harmonic scalar three-form. Root-string bounds on its extremal vector, the bracket on the three-form, and exchange of the dotted species give a finite rational kernel. Independently, the Weyl character formula for its tensor multiplicity, with the symmetric square on equal branches, gives the same dimensions. Highest and lowest vectors generate each ghost harmonic. Moving raising and lowering operators onto the other arguments with the invariant form recovers the entire functional from these extremal values. For the mixed-jet test the implementation moves the root action onto the coefficient three-form with the transposed minus sign and tests the whole resulting tensor. [Equations and exact linear algebra](first_page_tower.py) · [Probe coordinates and the seventh-power rational solution](tower_probes.json)

- For $p=7$ the free dimensions are $43+16+9=68$. The 68 independent equations fix the stored rational solution. All 1,117 differentiated-scalar components vanish, including 1,074 in which the particular term is nonzero. The four extremal pairs for the mixed second ghost give $6920+2554+2554+6920=18948$ vanishing components. Same-dot second jets vanish in the bracket quotient. Together these exhaust the incoming profiles, so $[S^7]\ne0$ on $E_2$.
- For $p=8$ the free dimensions are $2+1+22+7+22+37=91$. Append the particular column on 92 differentiated-scalar equations. The $91\times91$ homogeneous minor and the $92\times92$ augmented minor are $(156,375)$ modulo 1009 and $(382,29)$ modulo 1013. They are reductions of the same rational equations. The augmented rank is 92 and the equation has no solution over $\mathbb Q$. Reductivity allows invariant projection of any putative solution. Thus some classical cycle has $Q_1$ pairing nontrivially with $\lambda$, and the one-dimensional target is hit by $d_1$.

For the fifth power the only admissible second source has dimension $141-80-61=0$. The Taylor selection rule permits incoming orders at most two, and the survival of $S^6$ excludes an incoming $d_1$ on $S^5$ by multiplication by $S$. For the sixth power the independent target pairing is $-15313094596800$. The seventh-power product statement is conditional on finding a source. If $d_2[X]=[S^6]$, the product theorem gives $d_2([X][S])=[S^7]$ since $S$ is closed. The untwisted interpretation requires the quantum comparison in Appendix B. The audit contributes one independently confirmed $d_1$ target, zero verified $d_2$ lifts, and zero new fortuitous generators. [Fifth-power proof and sixth-power target](../first_page/README.md#the-fifth-glueball-survives-to-perturbative-infinity) · [Appendix B](https://arxiv.org/abs/2611.24160)

## Replay

Install NumPy, SciPy, SymPy and python-flint. From this directory, run the following against the archive at commit `783cdac`. Its SHA-256 is `b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972`. The first command reads the archived individual pole evaluations and reconstructs the determinants. The second repeats all distinct native contractions over $\mathbb Q$ and the tree-swap test. The third independently rebuilds all 703 outgoing evaluations in bounded integer arithmetic modulo the product $1009\cdot1013$. The last rebuilds the seventh- and eighth-power Ward tests from the repository’s native structure constants. [Expected output](results.json)

~~~sh
python3 verify_crossed.py /path/to/f4-j36-eight-Albert-certificate.zip
python3 verify_crossed.py /path/to/f4-j36-eight-Albert-certificate.zip --native
python3 check_outgoing.py /path/to/f4-j36-eight-Albert-certificate.zip --workers 4
python3 first_page_tower.py
~~~
