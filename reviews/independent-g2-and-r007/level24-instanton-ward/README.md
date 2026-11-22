# A superconformal Ward test for the charged $G_2$ instanton channel

In undeformed four-dimensional $\mathcal N=4$ Yang–Mills theory, the proposed level-24 instanton transition in $G_2$ has rank zero if the physical Taylor Ward generators act on perturbative cohomology with their stated leading terms and the graviton generators admit the exact primary product with its classical leading term. The matched microscopic instanton integral remains unevaluated. [Perturbative endpoints](https://github.com/lbotao399-stack/quantum-lifts-results/blob/3843a42d05f064b9fb855571c03917707dee9a20/reviews/independent-g2-and-r007/excess-two-ward/README.md) · [Physical primary product](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

Each tuple lists the three flavour degrees followed by the two dotted Taylor degrees. The candidate source at $(5,2,2;1,1)$ and length $Y=9$ survives every perturbative page. Its physical outgoing block has the unique perturbative line at $(3,0,0;3,3)$, $Y=8$. The source has flavour Dynkin labels $[3,0]$ and doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(2,2,12,6,6)$. The target has $(1,1,12,6,6)$. The degrees and letter length of separate summands may change at fixed physical charges. [Source ranks and first corrections](https://github.com/lbotao399-stack/quantum-lifts-results/tree/e965e8015f72466a8f26fb7993af53baa0bcd16d/reviews/independent-g2-and-r007/g2-filtered)

On the target replace one minus-dotted derivative by flavour 3 and a second minus-dotted derivative by flavour 2. Write $L$ for these two successive substitutions and $R$ for their two reverses in the opposite order. On divided jets the lowering has the exterior sign for moving the restored flavour past the larger flavours. The reverse replacement has this sign times the incremented minus-dotted multiplicity. Both act on products by the graded Leibniz rule.

For the first archived product $P$, the numbers of nonzero defining-trace words at the five successive steps are $4,11,21,11,4$. Direct addition over the rationals gives

$$
RL P=6P .
$$

The equality holds on all eight archived recipes, including a zero polynomial, before imposing finite-colour identities. The endpoint calculation proves $[P]\ne0$ on perturbative cohomology. If the physical Ward actions induce these leading Taylor actions, $L$ is injective on the unique perturbative target line. [Product and return](return_product.py) · [All eight computations](return.json)

On the source the two same-spin substitutions are alternating in their flavour labels. The allowed flavour components follow from

$$
[3,0]\otimes[0,1]=[3,1]\oplus[2,0].
$$

Let $t$ count plus-dotted derivatives in the lowered weight. At minus degree $t-2$ the possible flavour triples are $(7-t,4-t,3-t)$ and $(6-t,4-t,4-t)$, with $Y=10-t$. Nonnegativity leaves exactly the five rows below. The last two columns give rational rank lower bounds. The joined columns are tree boundaries and rational products of graviton components in the defining seven-dimensional trace.

| Lowered flavour and dotted degrees | $Y$ | Dimension | Joined rank | Outgoing tree rank |
|---|---:|---:|---:|---:|
| $(5,2,1;2,0)$ | 8 | 1206 | 702 | 504 |
| $(4,2,2;2,0)$ | 8 | 1397 | 911 | 486 |
| $(4,1,0;3,1)$ | 7 | 674 | 388 | 286 |
| $(3,1,1;3,1)$ | 7 | 1431 | 895 | 536 |
| $(2,0,0;4,2)$ | 6 | 183 | 110 | 73 |

All products are tree closed. In each row the joined rank plus the outgoing rank equals the exact cochain dimension. The products therefore span the tree cohomology over the rationals. At $t=3$ both cohomologies are zero. The replays use primes $2147483647$ and $1000033$ and recompute determinants for the product-augmented ranks. [Verifier](verify.py) · [Ranks and determinants](rebuilt.json)

The marked perturbative product supplies closed lifts of the product columns. Under the stated physical product hypothesis, their survivors admit representatives closed under the full supercharge. They span the lowered physical source after taking its successive perturbative grades. For a perturbative source $[X]$ at the original charge, the leading instanton differential $\delta$ satisfies

$$
L\delta[X]=\delta L[X]=0.
$$

The Ward identity on cohomology drops commutators with the perturbative supercharge arising from instanton corrections to the generators. The composition of the two odd substitutions is even. Injectivity on the permanent target line forces $\delta[X]=0$ in this channel. [Perturbative product and fan kernels](https://arxiv.org/abs/2611.24160) · [Prior Ward strategy in the symplectic theory](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/sp3_omega)

## Replay

Use Python 3.13, NumPy, SciPy, SymPy, python-flint and a C++17 compiler. Set REPO to the public trace runtime at [ee9e064](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7), CHARGED to the [charged product verifier](https://github.com/lbotao399-stack/quantum-lifts-results/tree/d8fbf518ce4448ece8c844f5291e02a98e3a0f4f/reviews/independent-g2-and-r007/g2-level24-charged-products), and ENDPOINT to the [rational endpoint packet](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3843a42d05f064b9fb855571c03917707dee9a20/reviews/independent-g2-and-r007/excess-two-ward).

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify.py --repo "$REPO" --charged "$CHARGED" --endpoint "$ENDPOINT" --rebuild | diff -u rebuilt.json -
~~~
