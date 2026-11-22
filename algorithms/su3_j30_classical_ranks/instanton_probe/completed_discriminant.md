# Discriminant clearing for the completed SU(3) R-007 class

Let $F$ be the charge-conjugation-odd ten-letter R-007 source of flavour $(7,4,4)$, $u$ the perturbative parameter and $\widehat F$ any homogeneous odd-under-charge-conjugation completion closed under the perturbative supercharge $Q$. Write $\Delta$ for the characteristic discriminant of $\phi^1$, of flavour $(6,0,0)$, and use the marked product for multiplication on $H(Q)$. Write $\mathcal G$ for the subalgebra generated over $\mathbb Q[u]$ by the closed single-graviton representatives, using the coefficient field of the marked product if an extension is required. Then for some finite $K\ge6$ and graviton polynomials $G_1,G_3$ there is an identity in perturbative cohomology
$$
[\Delta]^K[\widehat F]=u[G_1]+u^3[G_3].
$$
Here $G_1,G_3$ denote polynomials in the *completed* generators, including subsequent powers of $\Delta$. This proves that the completed class lies in the graviton algebra after localization. Its bare leading term has the explicit classical annihilator $\Delta^6$. The argument below gives existence of $K$ for the completed class. [Classical $\Delta^6$ identity](discriminant_ward.md) · [All-order completion of $F$](https://github.com/lbotao399-stack/quantum-lifts-results/blob/2e22109bc140b10e6a9bf70892dcf6b586510475/nonperturbative-lift/contributions/dot-r007-all-order/ALL_ORDER_PROOF.md)

## Five classical reductions

For every classical trace cocycle $A$, classical localization gives a nonnegative integer $k$, a polynomial $G$ in classical gravitons, and a polynomial cochain $H$ with $\Delta^k A=G+Q_0H$. When $A$ vanishes on the Cartan one can take $G=0$. Use the homogeneous components of the identity and project to the charge-conjugation parity of $A$. Each $G$ has a closed marked-product lift. The same construction gives $\Delta$ as a closed scalar. Only the first flavour occurs in $\Delta$, so its representative and all its powers have no positive-order corrections. [Classical localization, §6.2, (109)–(112)](https://arxiv.org/html/2610.22332v3#S6.SS2) · [Perturbative product and degree rules](https://arxiv.org/abs/2611.24160)

Start with the closed remainder $\widehat F$. At a reduction write $A$ for its order-zero term after removing the power of $u$ already extracted. Choose $k,G,H$ as above. Multiply the remainder by $\Delta^k$ using the marked product, subtract the closed lift of $G$ and subtract $QH$. The result is closed and its order-zero coefficient vanishes. Dividing this **cochain** by $u$ gives the next closed remainder since the cochain module is free in $u$. This argument permits $u$-torsion in cohomology. At the initial step use $k=6$ and $G=0$.

After $r$ divisions, let $K_r$ count the discriminants already used. The leading term has flavour and dotted degrees and letter number
$$
(7+6K_r-r,\,4-r,\,4-r;\,r,r),\qquad Y=10+6K_r-r .
$$
Indeed an order-$s$ marked-product graph has $s$ cubic vertices and $2s$ physical edges. Each vertex supplies and each contracted pair consumes one copy of every flavour, giving a loss of $s$ in each. Its remaining letter number is reduced by $s$ since $3s-2(2s)=-s$. The physical spin gradings then insert $s$ dotted derivatives of each sign. A correction $QH$ has the same homogeneity with $H$ chosen one letter shorter at its leading order. Multiplying by any additional $\Delta$ adds six first-flavour slots and six letters. The total number of transverse flavour and dotted indices stays $(4-r)+(4-r)+r+r=8$. [Marked graph expansion and quantum charges](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex)

| $r$ | Second and third flavour degrees | Dotted degrees | Sign of charge conjugation on a Cartan monomial of this length | Possible graviton subtraction |
|---:|---|---|---:|---|
| 0 | $(4,4)$ | $(0,0)$ | $+1$ | zero |
| 1 | $(3,3)$ | $(1,1)$ | $-1$ | $G_1$ |
| 2 | $(2,2)$ | $(2,2)$ | $+1$ | zero |
| 3 | $(1,1)$ | $(3,3)$ | $-1$ | $G_3$ |
| 4 | $(0,0)$ | $(4,4)$ | $+1$ | zero |

On diagonal traceless jets negative transpose negates each letter. Thus the Cartan component of a charge-conjugation-odd remainder vanishes whenever $Y$ is even. This eliminates $G$ at $r=0,2,4$. After the fifth division both transverse flavour degrees would be $-1$, with all subsequent coefficients having still smaller degrees. The remainder is identically zero. Unwinding the five equalities gives the displayed identity with $K$ the sum of the five clearing exponents. In subsequent multiplications the graded Leibniz rule preserves exact terms, and associativity holds in cohomology.

## Consequence for an instanton action

Suppose the first matched instanton differential acts as a $u$-linear derivation on perturbative cohomology, or on it after inverting $u$, and annihilates the graviton lifts used above. Denote this differential by $D$. Applying it to the identity gives
$$
[\Delta]^K D[\widehat F]=0 .
$$
Any image with Cartan-visible leading class survives multiplication by every power of $[\Delta]$. Indeed its leading Cartan polynomial is multiplied by the nonzero scalar discriminant, which acts injectively on each exterior coefficient. A putative bounding equation, including one after inverting $u$, would make this product an incoming perturbative boundary on a finite page. Cartan evaluation annihilates each such boundary. This contradiction excludes a Cartan-visible leading image under the stated physical hypotheses. The argument uses the completed source and requires no finite bound on the number of inputs in a contact. [Perturbative Cartan-image theorem](../../perturbative_cartan_targets.md)

For Dot's outgoing nine-letter representative the exact three-jet value is $-14112/5$ at $\Delta=400$. It therefore cannot by itself be the leading image with nonzero coefficient under these hypotheses. In the same $(5,2,2;2,2)$ block the classical Cartan kernel is one-dimensional. Four successive classical Taylor lowerings take its generator to $240$ times a nonzero nine-letter BMN class. Thus if the matched primary Ward identity also requires all these lowerings to annihilate the leading image in classical cohomology, the image in this block is zero. A mixed leading image is constrained by its total Cartan polynomial. Establishing the physical product matching for the instanton insertion, or finding a channel with a nonzero physical matrix element, remains open. [Three-jet evaluation](discriminant_ward.json) · [One-dimensional kernel and nonabelian return](../bmn_return/README.md)
