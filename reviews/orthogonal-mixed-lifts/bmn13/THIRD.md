# A second independent third-page image in the mixed SO(13) terminal space

The full relative Taylor-jet complex of undeformed four-dimensional $\mathcal N=4$ SYM uses the defining orthogonal trace and the ordered-pair kernel $Q_1$ with $K=1$. The two derivative-free product primitives $O_0^{(0)},O_0^{(1)}$ have the seed order and highest-component projection of the [lower-page certificate](https://github.com/lbotao399-stack/quantum-lifts-results/tree/e9c48f4842a0d68af7fe62f9ba7085a922af112b/reviews/orthogonal-mixed-lifts/bmn13). Their source has flavour $(7,5,5)$, length sixteen and level $j=34$. Denote by $D$ the ordered product of the six distinct substitutions of a flavour derivative by a dotted derivative on divided jets. Apply $(\dot+,1),(\dot+,2),(\dot+,3),(\dot-,1),(\dot-,2),(\dot-,3)$ in that order. The resulting descendant has level $j=40$. [Jet conventions and the supercharge](https://arxiv.org/abs/2611.24160)

The terminal dotted and flavour degrees of the descendant's third target are $(6,6;2,0,0)$. In the three ordered trace coordinates of this sector, let $\lambda=(6,-1,0)$. Write $S=\operatorname{Tr}_{13}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$. The two substitutions of a dotted derivative of $S^7$ by flavour 1 have coordinates $(7,42,0)$, which $\lambda$ annihilates. For the class surviving the first obstruction the exact result is

$$
\lambda\,d_3\!\left(D\left(358475625[O_0^{(0)}]-832[O_0^{(1)}]\right)\right)
=-50677394964480000.
$$

The orthogonal $S^7$ differential and its source are due to Dot. The displayed calculation detects an image independent of the vector $(7,42,0)$. It also detects a third differential on the $j=34$ parent since the six substitutions commute with the differential on its correction equations. The full terminal space before finite trace identities has dimension three. Cartan restriction has rank one, with coordinates $(-92160,15360,-7680)$ in the same order. To transport Dot's differential from level 42 down to level 40 one needs the filtered superconformal Ward action for the two substitutions in the opposite direction. Under that hypothesis the two third images span the Cartan kernel and this terminal space gives no excess perturbative survivor for the proposed comparison to $USp(12)$. The rational pairing above and the exclusion of its own line from a nonperturbative mismatch use only the directly verified upward action. [Dot's pure-power result](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md) · [Filtered Ward hypothesis](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex) · [Cartan target criterion and hypotheses](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e9c48f4842a0d68af7fe62f9ba7085a922af112b/reviews/orthogonal-mixed-lifts/bmn13/README.md#cartan-visible-targets-at-all-perturbative-orders)

## Descent to the third page

The previous two-prime minors give two independent source classes modulo stable cycles. The first differential has rank one, with the kernel displayed above. No positive-order differential enters the zero-derivative source. The second target has odd total length fifteen. Reversing an odd orthogonal trace modulo adjacent graded commutators makes every terminal monomial at this length zero. All three source and correction grades have one second jet and otherwise first jets. For $r\ge2$, the $r+1$ contracted inputs of $Q_r$ contain at most $r+2<3r$ flavour slots. Thus only $Q_0$ and $Q_1$ enter these equations.

Write $L_1,L_2$ for the first dual and corrected second dual. Their defining equations are $\lambda Q_1=L_1Q_0$ and $-L_1Q_1=L_2Q_0$. The stable identities cover 314 and 16,978 potentially contributing native columns. The first relation space at the preceding weight vanishes in the required flavour representation. On the second dual subtract the commuting distinct-plane functional with coefficients $(15/128,33/32,4/3,27/64,7/12)$ on zero-based full columns $(0,1,2,10,11)$. The residual vanishes on the 76 determinant polarizations of the lower-page certificate after projection to the highest flavour component. The path reduction there proves completeness modulo tree commutators within the relations. The second dual now evaluates finite-colour cochains. For any lower corrections $Q_0O_1=-Q_1O_0$ and $Q_0O_2=-Q_1O_1$ in this quotient,

$$
\lambda Q_1O_2
=L_1 Q_0 O_2
=-L_1 Q_1 O_1
=L_2 Q_0 O_1
=-L_2 Q_1 O_0 .
$$

After the six substitutions and highest-component projection the two values of $L_2Q_1$ are $455081984/3$ and $4448563200000$. Only their specified linear combination has a third differential. The residues of its signed third pairing modulo $1000003$ and $1000033$ are $248799$ and $367978$.

## Odd substitutions on the correction equations

The native flavour-to-dotted substitution removes its chosen flavour from a jet, with the exterior sign for passing the larger flavours and with the divided-jet factor equal to one plus the old degree in that dotted direction. It is an odd derivation on trace words. All twelve second jets and all five first jets obey the anticommutator with $Q_0$. In this sector a $Q_1$ pair can only be a matter second jet and its complementary scalar. It consumes exactly three distinct flavours and produces two gaugino first jets. Substituting a contracted input destroys the pair. Substituting an output vanishes. For a spectator the two orders of the two odd operations have opposite signs. The trace reconnections are even tensor permutations with their Koszul signs. This proves the anticommutator on every trace word of the three grades. The six substitutions anticommute among themselves. Their product is even and is a determinant character in flavour and dotted space, so flavour highest components and dotted singlets are preserved. The code checks 102 classical letter identities, 612 substitution anticommutators, 870 ordered pair tests and 1,770 complete native trace equations.

## Independent arithmetic

The modular pass takes the transpose of each trace reconnection. A split rejoins two oriented output traces. A join cuts one output trace into two oriented pieces, each containing one ghost. The orthogonal term reverses one arc of a trace. A contracted empty arc supplies the matrix dimension thirteen. Both identifications of the ghosts with the matter and scalar inputs and all three missing flavours are enumerated. An active trace in the chosen primitive has even length, and a separate scalar bilinear dresses it. These conditions give 42,298 candidate parents, including parents with zero coefficient or zero contraction.

The rational replay obtains the flavour projector independently. On the eight labelled scalar slots of weight $(4,2,2)$, the sum of all transpositions has eigenvalue equal to the sum of box contents (column minus row) on each symmetric-group shape. The target shape has eigenvalue 2. The other shapes dominating this weight have distinct eigenvalues from the set $\{4,7,8,10,12,14,20,28\}$. Spectral interpolation therefore projects to the required highest component. The native odd substitutions applied to this projected covector reproduce all 9,053 rational word coefficients modulo both primes. The rational covector has denominator three in the ordered word basis.

Independently, the forward source expansion enumerates every raw three-flavour scalar word of lengths $3,5,7,9,11$. For each word, multiply its actual auxiliary two-by-two matrices, take the trace, reverse the colour word by successive adjacent transpositions, and accumulate its elementary commutator primitives. At odd length $l$ the denominator is $4l$, from the half trace logarithm and the half reversal. In a product average the primitive over its odd trace factors, expand the exponential of the remaining connected traces with multiplicity factorials, and multiply by the specified scalar bilinear. This gives rational coefficients on 11,994 distinct active traces. The pinned native forward $Q_1$ independently recomputes every parent contraction, including the zeros. The two source coefficients and the contraction agree at every parent at both primes. Exactly 34,966 parents contribute to at least one seed.

## Replay

Use Python 3.13, NumPy, SciPy, python-flint and a C++17 compiler. Supply a checkout of the [trace runtime at ee9e064](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7). From the directory containing the scripts run

~~~sh
python3 replay_bmn_third.py --repo "$REPO" --workers 2
python3 replay_bmn_third.py --repo "$REPO" --workers 2 --rebuild-adjoints
~~~

The second command also regenerates the stable adjoints, the finite determinant rows, and the native adjoint checks. The rank and first-differential certificates have their own replay in the lower-page directory. The powers of $g_{\rm YM}$ for a physical supercharge matrix element and its square are respectively four and eight under Appendix B's comparison and a nonsingular radial metric. The calculated coefficient is a supercharge obstruction in the stated trace and source conventions.
