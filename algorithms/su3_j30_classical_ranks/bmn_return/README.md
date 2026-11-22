# The rational target line and its incoming differential at level 30

In the $SU(3)$ relative polynomial Taylor complex, the nine-letter highest-weight block with dotted degrees $(2,2)$ and flavour degrees $(5,2,2)$ has exactly one classical class restricting to zero on the Cartan. The incoming $d_2$ has rank zero over $\mathbb Q$ on the entire source block.

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
