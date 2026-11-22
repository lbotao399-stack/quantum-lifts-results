# Cartan restriction with one second jet in SU(2)

In the full relative holomorphic jet complex of undeformed $\mathcal N=4$ SYM with gauge group $SU(2)$, Cartan restriction is injective on the entire classical cohomology containing one second jet and otherwise first jets, at every Taylor degree. Every class is a sum of products of components of $\operatorname{Tr}(dC)^2$ with one Taylor derivative of such a component. At left spin $J_L=0$ or $1/2$ every classical class survives every perturbative page. A nonzero supercharge differential at any positive order in this group requires target $J_L\ge 1$ and source $J_L\ge 3/2$. The proof uses the classical tree complex and the perturbative Cartan-image and product theorems. The physical identification of the formal complex uses Appendix B of the all-loop manuscript. [Perturbative image theorem](../../perturbative_cartan_targets.md) · [Supercharge grading and product theorem](https://arxiv.org/abs/2611.24160) · [Classical relative complex](https://arxiv.org/abs/2209.06728)

Write $Y$ for letter number and $n$ for total Taylor degree, counting all dotted and flavour derivatives of the relative superfield $C$. Then $2J_L=n-Y$. At $n-Y=1$ there is exactly one second jet. Every $Q_r$ changes $(n,Y)$ by $(-r,1-r)$, so its induced obstruction lowers $2J_L$ by one. At $n-Y=0$ the tree cohomology is the ring of invariant functions of the supercommuting first jets. Its Cartan restriction is injective for all $SU(N)$ and $USp(2N)$. For $GL_N$ and $Sp_{2N}$ this follows from the ordinary commuting-scheme restriction theorem applied to arbitrarily many labelled slots, then polarization with the physical signs. For $SL_N$ split every matrix into its scalar and traceless parts. Both ordinary restriction rings factor by the same polynomial ring in scalar parts, and setting those parts to zero gives the traceless isomorphism. This already removes all terminal targets at all ranks in those two series. [Vaccarino, Theorem 3](https://arxiv.org/abs/math/0602660) · [Chen–Ngô, Theorem 4.1](https://arxiv.org/abs/2102.01849) · [Project's multilinear transfer, section 5](../../../notes/Gram_determinants_commuting_schemes.tex)

## A bound on generators

Start with an arbitrary space $V$ of ordinary labels for even first jets. In the complexified adjoint of $SU(2)$ the Lie bracket is the three-dimensional cross product in the standard $\mathfrak{so}_3$ realization. Write $X_i$ for a first jet, $Y_{ij}$ for its alternating second jet and $Z_{ijk}$ for the alternating third jet. Subscripts label Taylor directions. With the alternating wedge order the relevant tree maps are

$$
Q_0Y_{ij}=[X_i,X_j],\qquad
Q_0Z_{ijk}=[X_i,Y_{jk}]-[X_j,Y_{ik}]+[X_k,Y_{ij}].
$$

The first-jet polynomial ring is the polynomial ring on a three-row matrix of colour components. The entries of $[X_i,X_j]$ are its signed two-by-two minors. The cycles with one $Y$ are precisely the polynomial relations among these minors. The boundaries come from two $Y$'s or one $Z$. The two-$Y$ boundaries give the Koszul relations, obtained by multiplying one minor by another and subtracting the reverse product.

Modulo these Koszul relations the relations are generated in polynomial degree three, counting a $Y$ as degree two. Indeed the minors form a quadratic Gröbner basis for the rank-one matrix ideal. Order each monomial by its row labels and sort column labels by straightening inverted pairs with a minor. The sorted monomial is unique for each pair of row and column multisets, as is also seen by substituting a column vector times a row vector. For two leading minor monomials with a common factor the overlap has degree three, so reducing its difference to zero gives a cubic relation. For relatively prime leading monomials the Koszul relation already has the required leading cancellation. These relations generate by leading-term descent on an arbitrary relation, equivalently by the Gröbner syzygy theorem. This argument also covers one or two label directions. [Segre syzygies and prior linear-presentation results](https://arxiv.org/abs/math/0404417)

Each minor acts by zero on the one-$Y$ cohomology, since its product with a cycle is the tree variation of the corresponding $Y$ times that cycle. In the quotient by the minor ideal, degree $k$ scalar coefficients have colour and label space $\operatorname{Sym}^k(\mathbb C^3)\otimes\operatorname{Sym}^k V$, by the same rank-one factorization. Every cubic relation lies inside the tensor product of two adjoint colour vectors, one from $X$ and one from $Y$. Its possible colour spins are $0,1,2$. In $\operatorname{Sym}^k(\mathbb C^3)$, the spin $l$ part is a traceless tensor of degree $l$ multiplied by powers of the colour quadratic form, with $k-l$ even and nonnegative. Multiplication of symmetric label tensors onto $\operatorname{Sym}^k V$ is surjective. Taking colour invariants is exact in characteristic zero. Consequently the invariant one-$Y$ cohomology, including the quotient by the $Z$ boundaries, is generated in Taylor degrees $3,4,5$ over the scalar invariant ring. On the Cartan this ring is the even polynomial algebra $R=\bigoplus_{t\ge0}\operatorname{Sym}^{2t} V$, where its degree-one generators are scalar bilinears.

## Four finite integer calculations

Use multilinear Taylor labels. In degree $n$ choose the two labels of $Y$ and contract its colour and the remaining first-jet colours by invariant tensors of $SO(3)$. Products of metrics and at most one volume tensor span these tensors. Through six even vector slots all metric matchings are independent, as independently colouring each matched pair with its own Cartesian axis extracts its coefficient. The three-vector invariant is the volume. In the preceding cochain space the program uses every triple of labels on $Z$ and every disjoint pair of label pairs on two $Y$'s. The signed tree maps use the displayed bracket and the odd product rule. The outgoing rank is evaluated on every assignment of three coordinate axes to the scalar colours. All reported ranks are integer-matrix ranks over $\mathbb Q$. The program checks both $Q_0^2=0$ and zero Cartan value on every incoming column. [Integer matrices and replay](verify.py) · [Rank records](ranks.json)

| Taylor slots | Middle dimension | Incoming rank | Outgoing rank | Classical dimension | Cartan rank of displayed cycles |
|---:|---:|---:|---:|---:|---:|
| 3 | 3 | 0 | 1 | 2 | 2 |
| 4 | 6 | 4 | 2 | 0 | 0 |
| 5 | 30 | 14 | 6 | 10 | 10 |
| 7 | 315 | 258 | 36 | 21 | 21 |

For three distinct labels $i<j<k$ use the invariant form $B$ on the adjoint. The three terms $B(X_i,Y_{jk}),B(X_j,Y_{ik}),B(X_k,Y_{ij})$ have tree variations with signs $+,-,+$ times the same invariant triple. The first plus the second and the second plus the third are closed. These are signed Taylor derivatives of $B(X_i,X_j)$. The Cartan ranks in the last column use just these two cycles times all pairwise contractions of the remaining first jets. The checks in degrees three through five therefore prove that the cubic cycles alone generate over $R$.

## Why degree seven suffices for injectivity in every degree

Write $x_i,y_{jk}$ for the Cartan entries of the first and second jets. The Weyl reflection negates both. With one $y$ the Cartan module is
$\bigoplus_{t\ge0}\operatorname{Sym}^{2t+1}V\otimes\wedge^2 V$.
Give its cubic Taylor part module degree zero. The cubic cycles restrict to the kernel of alternation from $V\otimes\wedge^2V$ to $\wedge^3V$. Call the submodule they generate $M$. Modulo the cubic relations the tensor $x_i y_{jk}$ is fully alternating in its three labels. Temporarily allowing arbitrary scalar polynomial coefficients, the three interchanges give

$$
x_a x_b y_{cd}=-x_a x_c y_{bd}
= x_b x_c y_{ad}=-x_b x_a y_{cd}.
$$

Thus a Cartan monomial with at least two $x$'s is a sum of multiples of the cubic cycles. With an odd number at least three of $x$'s the multiplier has even degree and belongs to $R$. The quotient of the Cartan module by $M$ is therefore $\wedge^3 V$ in module degree zero, with every positive-degree scalar acting by zero.

The full Cartan module has a presentation over $R$ with generators $x_i y_{jk}$ and relations in module degree one. Indeed fixing $y_{jk}$ leaves the odd polynomials in the $x$'s. Any two encodings of the same monomial by scalar pairs and a singleton agree after swapping the singleton with an entry of a pair and re-pairing the other entries. The even ring $R$ itself is defined by quadratic re-pairing relations, since these swaps connect all pairings of any fixed multiset.

For completeness, $\operatorname{Tor}_i^R(-,\mathbb Q)$ records the degrees of the generators of step $i$ of a minimal graded free resolution over $R$, with $\mathbb Q=R/R_{>0}$. The linear presentation just proved places the first Tor of the full Cartan module in degree one. The second Tor of $\mathbb Q$ over $R$ is in degree two. For this latter assertion, present $R$ by a polynomial ring on the scalar pairs with a quadratic ideal. The relations on its augmentation generators are generated by the usual degree-two commutativity relations and lifts of the defining quadrics, by lifting a relation to the polynomial ring and subtracting these quadrics. In the exact Tor sequence for $M$ inside the full Cartan module, its first Tor lies between the second Tor of the quotient and the first Tor of the full module. Hence the defining relations among the cubic generators of $M$ occur only in module degrees one and two, Taylor degrees five and seven.

Cartan restriction gives a surjection from the actual one-$Y$ cohomology to $M$, with the same cubic generators. The degree-five and degree-seven equality of the dimensions and Cartan ranks proves that every generating relation of $M$ already holds in the actual cohomology. This proves injectivity. Multilinear tests cover arbitrary repeated labels by polarization in characteristic zero. For the physical three even and two odd first jets, tensor the multilinear Taylor modules with the signed tensor powers of their label space and take permutation coinvariants. This operation is exact in characteristic zero and commutes with splitting a jet. It restores the physical symmetric dotted and alternating flavour Taylor derivatives. The argument needs Taylor covariance of $Q_0$ only.

Thus at odd $n\ge5$ the full one-second-jet module on first-jet labels is $\operatorname{Sym}^{n-2} V\otimes\wedge^2 V$, using graded powers for physical labels and the intrinsic parity of the second jet. At $n=3$ it is the kernel of the alternation above. At even degree it is zero. Four dotted slots require both different dotted labels in the graded symmetric power and both slots of the graded exterior square. Their dotted representation is the product of the alternating square of the doublet with its symmetric square, which has spin one. In particular the dotted singlet with four dotted slots vanishes identically in this sector, for every scalar flavour weight. The incoming BMN second differential has precisely two dotted slots of each sign.

The perturbative Cartan-image theorem rules out every incoming differential on the one-second-jet sector. Its outgoing differentials land in the Cartan-injective terminal sector. Equivalently, the identified products of bilinear gravitons have closed quantum completions by the product theorem. This proves perturbative survival. A physical instanton matrix element requires additional input. [Cartan-image proof](../../perturbative_cartan_targets.md) · [Graviton generators and quantum product](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L4695)

## Direct test in the first remaining target spin

For a BMN source the second supercharge obstruction adds two derivatives of each dotted sign. In the dotted singlet, the full trace engine tested every dominant flavour partition of scalar degree zero through nine, at target $2J_L=2$ and at $2J_L=1$ as a control. These are 53 flavour partitions, covering levels $12,14,\ldots,30$. At $2J_L=2$ the first prime gives total Cartan image ranks by level $0,0,3,0,4,0,4,0,4,0$. At $2J_L=1$ all classical blocks vanish. Every individual rank is reproduced at a second prime. In each block the classical cochain dimension is the incoming rank plus the rank of the outgoing map stacked with Cartan restriction; the latter rank equals the outgoing rank plus the restriction rank on its kernel. These integer-character dimensions and modular lower bounds force zero Cartan kernel over the rationals in each block. Thus no incoming second obstruction reaches these 53 targets at $2J_L=2$. The argument proves the kernel bound even if separate modular cohomology dimensions are only upper bounds. [Rank records](bmn_targets.json) · [Native engine](../taylor_hooks/native/engine.cpp)

Run from this directory with NumPy and python-flint installed.

~~~sh
python3 verify.py | diff -u ranks.json -
c++ -O3 -std=c++17 -pthread ../taylor_hooks/native/engine.cpp -o /tmp/su2-terminal-engine
python3 bmn_targets.py /tmp/su2-terminal-engine | diff -u bmn_targets.json -
~~~
