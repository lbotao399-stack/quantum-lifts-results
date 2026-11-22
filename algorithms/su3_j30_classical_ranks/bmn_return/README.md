# The rational target line and its incoming differential at level 30

In the $SU(3)$ relative polynomial Taylor complex, the nine-letter highest-weight block with dotted degrees $(2,2)$ and flavour degrees $(5,2,2)$ has exactly one classical class restricting to zero on the Cartan. The incoming $d_2$ has rank zero over $\mathbb Q$ on the entire source block. The first-order reverse Taylor identity below also excludes incoming $d_1$. The Cartan-kernel line survives to perturbative $E_\infty$.

At dotted degree zero, start with flavour $(5,4,4)$ and nine letters. Separate the cochains by charge conjugation $C$, which commutes with the classical differential $Q_0$. Independent integer characters give the cochain dimensions. Nonzero modular minors of the incoming and outgoing evaluations give rational lower bounds on their ranks. The resulting upper bounds on cohomology are two at even length eight for $C=+1$, one at odd length nine for $C=-1$, and zero at every other length. The two Euler characteristics, with length sign, are $2$ and $-1$. Both upper bounds are attained. In particular the odd cohomology at length nine has dimension $708-356-351=1$. Its Cartan highest-weight dimension is zero. [Character and rank checks](verify.py)

The four defining-trace factors in [the product recipe](product.json) give 477 terms after highest-weight projection. Appending this product to the classical boundaries raises the total rank from 658 to 659. The two augmented determinants are 563838 modulo 1000003 and 242547 modulo 1000033. The Euler saturation fixes the unaugmented rational rank. Four signed Taylor substitutions take the product to 34492 terms in the required nine-letter block. The four reverse substitutions recover 240 times the original polynomial coefficient by coefficient. Each substitution anticommutes with $Q_0$ on generators and preserves finite-colour identities. Restriction to diagonal matrices commutes with these substitutions. The raised polynomial is therefore Cartan-zero and represents a nonzero rational class. The existing upper bound $63983-31838-32144=1$ proves that it spans the whole kernel. [Polynomial reconstruction](product.py) · [Large-block rank data](../README.md)

The only incoming $d_2$ source at these degrees has no dotted jets, flavour $(7,4,4)$ and ten letters. At lengths nine and ten, respectively, the minors for $Q_0$ give the following upper bounds in the two charge-conjugation sectors.

$$
2768-1113-1651=4,\qquad 2896-1073-1822=1.
$$

Four even products of graviton components have Cartan determinant $3528$ in ordinary traceless diagonal coordinates. The cyclic colour basis gives $2812385772$ with its own pivot columns. Both are determinants after clearing the denominators $7,7,14,14$. The four classes admit closed marked products by Li–Liu–Chang, Theorem 9.1. Dot's odd class R-007 has a nonzero classical block pairing and admits an all-order formal completion. These five independent classes exhaust the rational source. Thus $d_2=0$ on the whole block. The proof uses no quantum Taylor-lowering identity. [Graviton products](source-products.json) · [Dot's proof and input links](https://github.com/lbotao399-stack/quantum-lifts-results/blob/2e22109bc140b10e6a9bf70892dcf6b586510475/nonperturbative-lift/contributions/dot-r007-all-order/ALL_ORDER_PROOF.md) · [Marked-product theorem](https://arxiv.org/pdf/2611.24160v1#page=86)

With Python 3.13, NumPy, SciPy and python-flint, run from the repository root

~~~sh
python3 algorithms/su3_j30_classical_ranks/bmn_return/verify.py
python3 algorithms/su3_j30_classical_ranks/verify_chains.py
~~~

Append `--rebuild-native` to the first command to recompile the C++17 evaluators and regenerate all twelve local rank records, including their ordered minors. The default checks their coordinates and pins, recalculates the integer characters and both source product determinants, reconstructs the two classical representatives and checks the return. The large-block command replays the pre-existing target rank records.

## A quantum representative of the outgoing candidate

Write $T_0$ for the classical generator of the nine-letter Cartan kernel above. The finite first correction below solves $Q_0T_1=-Q_1T_0$ in the native trace convention. The marked product gives an all-order closed choice of representative of $[T_0]$. A physical instanton matrix element from R-007 to this line remains open. [Exact pair primitives](first_correction.json) · [Four-step transport](raised_first.json)

The four factors of the BMN product have lengths $2,2,2,3$ and parities $0,1,1,0$. Each obeys $Q_0F_i=Q_1F_i=0$ as an ordered trace polynomial. Since $Q_1$ acts on two inputs, its action on the four-factor product is the sum of the six pair images with their graded spectator signs. For each pair the replay merges adjacent letters in the support and solves the resulting rational tree equation. The pair sizes in lexicographic order are

| Factors | 12 | 13 | 14 | 23 | 24 | 34 |
|---|---:|---:|---:|---:|---:|---:|
| Terms in the first image | 3 | 6 | 24 | 18 | 150 | 68 |
| Parent columns tested | 5 | 12 | 42 | 18 | 144 | 62 |
| Terms in the primitive | 2 | 4 | 10 | 6 | 84 | 20 |

Reinserting spectators gives 192 terms. Projecting onto the original highest weight and odd charge conjugation gives 7,762 terms. Direct application of $Q_1$ to the 477-term BMN polynomial and $Q_0$ to this correction cancels all 21,176 ordered-trace coefficients. This identity already holds before imposing the finite-colour trace identities.

The forward Taylor substitution $P$ removes a specified flavour from a letter and increases its specified dotted degree by one, with the divided-jet coefficient and exterior sign in `product.py`. It is the odd shear that translates the dotted coordinate by the corresponding odd coordinate times an odd parameter. In the all-loop fan formula the functional derivatives and both output fields have a common $\theta$. Every output position is an affine combination of input positions with coefficient sum one. In the two displayed output points, the inner sums of even-indexed parameters telescope to one, followed by the outer sum of odd-indexed parameters. A common translation depending only on $\theta$ therefore shifts both endpoints by exactly that translation. The dotted derivatives commute with it. Equivalently, on independent input odd coordinates, every difference $\theta_I^{(j)}-\theta_I^{(k)}$ annihilates the product of all adjacent odd differences. Extracting the shear parameter gives $Q_nP+PQ_n=0$ for every $n$. The shear preserves the relative condition at the origin. [All-loop supercharge and endpoint positions](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L189-L233)

For each of the four forward substitutions we checked this identity on the ordered two-input $Q_1$ kernel for every pair of the 37 possible letters through total Taylor degree three. In each run 43 ordered pairs have a nonzero output action which cancels the input action. We also checked $Q_0P+PQ_0=0$ on all letters of the source and correction at each step. The four substitutions give a 34,492-term $T_0$ with the original hash and a **70,141-term** $T_1$ with common denominator $180$. Both orders use the existing source scale. [Operator replay](raised_first.py)

For all-order existence, the only factor with three letters has five flavour slots. A bare $Q_2$ requires six, and $Q_n$ for $n\ge3$ needs at least four letters. Thus the four factors are strictly closed. Apply the marked product, the equivariant highest-weight projection, and the four forward substitutions. The leading term is $T_0$. Its flavour degrees are $(5,2,2)$, so the homogeneous completion stops after its second correction. This existence argument allows the choice of first correction to vary by a tree cycle. It establishes outgoing perturbative closure. The reverse Taylor identity below settles the incoming first differential. The physical instanton coefficient remains open. [Marked-product theorem](https://arxiv.org/pdf/2611.24160v1#page=86)

~~~sh
cd algorithms/su3_j30_classical_ranks/bmn_return
python3 first_correction.py --check first_correction.json
python3 raised_first.py --check raised_first.json
~~~

## The nine-letter target survives every perturbative page

Let $\bar P_0$ lower one specified dotted derivative and restore the specified scalar flavour on its letter, with the divided-jet exterior sign of the reverse operation in [the product replay](product.py). For each of the four reverse operations in that replay we construct an odd two-input, one-output correction $\bar P_1$. On gauge-invariant polynomials the exact identity is

$$
\{Q_0,\bar P_0\}=0,\qquad
\{Q_1,\bar P_0\}+\{Q_0,\bar P_1\}=0.
$$

The second identity holds here on words for which every letter and every pair of distinct letter positions has total Taylor degree at most seven. Nine letters with fourteen slots satisfy this bound because the other seven letters of a pair each use at least one slot. All four reverse operations preserve the total and the length. The gauge contraction for $\bar P_1$ is the invariant colour bracket with both inputs raised by the native trace metric. For the shear adding flavour $i$ and lowering dotted direction $\alpha$, the correction consumes one copy of each of the other two flavours and adds one derivative in the other dotted direction.

In unrestricted colour coordinates the sum of anticommutators is a sum of compositions whose rightmost operation is a simultaneous infinitesimal gauge rotation. The remaining first-order operator changes the jet of one letter and rotates its colour. This term vanishes on every invariant polynomial. The sum of anticommutators has differential order at most two and annihilates constants. Values on individual coordinates and on their pair products therefore determine it on the full sector by the graded Leibniz rule for its bilinear symbol. Each output on each of these inputs is compared as a polynomial, with its exterior signs.

For each shear, 24,736 nonzero rational linear equations determine 846 pivot coefficients in the bracket and gauge terms. A8 job 5702721 used the ordinary defining-trace metric on the six matrix units and two traceless diagonals of $\mathfrak{sl}_3$. It substituted the supplied coefficients into every single-colour input and all 64 ordered two-colour choices, including equal and nonorthogonal colours. Each shear passes 1,786,154 component equations, giving 7,144,616 checks. A separate direct trace expansion, A8 job 5702502, cancels all 119,334 nonzero coefficients in $Q_1T_0+Q_0T_1$. [Rational coefficients](reverse_ward.json) · [Operator replay](reverse_ward.py) · [A8 outputs and hashes](reverse_ward_a8.json)

Suppose $d_1[X]=[T_0]$, where $d_1$ is induced by $Q_1$ on classical cohomology. The nine-letter input $X$ has dotted degrees $(1,1)$ and flavour degrees $(6,3,3)$. The second Ward identity lets each reverse operation pass through $d_1$ with a minus sign. Lowering twice in each dotted direction annihilates $X$. On $T_0$ the same ordered operations give 240 times the nonzero 477-term BMN class. Its two augmented determinants are 563838 and 242547, and Euler saturation fixes the rational boundary rank. This contradiction excludes every incoming $d_1$ to the line. The entire incoming $d_2$ was excluded above. An incoming order $r\ge3$ would require dotted degrees $(2-r,2-r)$. Together with the marked-product completion this proves $[T_0]\ne0$ on $E_\infty$ of the perturbative supercharge. The one-instanton supercurrent matrix element from R-007 to this line remains undetermined.

~~~sh
python3 reverse_ward.py solve /tmp/reverse-ward.json
cmp /tmp/reverse-ward.json reverse_ward.json
python3 reverse_ward.py full-sl3 reverse_ward.json --workers 8
~~~
