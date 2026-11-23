# The crossed eight-Albert source for the sixth F₄ glueball

The proposed certificate for $[S^6]\in\operatorname{im}d_2$ at $j=36$ fails. Here $S$ is the glueball with invariant form $\operatorname{Tr}_{26}/6$, and $d_2$ is the second differential of the supercharge filtration in the relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM. Both corrected augmented determinants vanish at the specified primes. Whether a source hits $[S^6]$ over $\mathbb Q$ remains open. Dot supplied the crossed source B and withdrew this certificate after the tree-order error was identified. The earlier source with a classical primitive is a separate case. [Crossed-source withdrawal](https://github.com/lbotao399-stack/quantum-lifts-results/blob/7fcfb7175703dd1ee6dcf4a71f2885c0d7dc8c73/reviews/f4-j36-eight-albert-d2/README.md) · [Earlier source](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1#issuecomment-7074126654)

## The rational graph calculation

The two lists `cp` and `cm` fix the exterior order of the six gauginos of each dotted species. The archived contraction visits the edges in a tree-dependent order. Restoring their declared order multiplies each term by the product of the two permutation signs. On left-marked graph 0, the edge visit order is $(15,8,10,13,16,14,17)$. The two lists have one and four inversions, giving $(-1)^{1+4}=-1$. Interchanging the two children of the same internal vertex gives the following exact pole evaluations. [Original graph-order report](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1#issuecomment-7105180683)

| Graph 0 | Archived tree | Tree with the children exchanged |
|---|---:|---:|
| Tree-order poles | $(-5750,24200)$ | $(5750,-24200)$ |
| Declared edge-order poles | $(5750,-24200)$ | $(5750,-24200)$ |

This replay contracts all $18+30-12=36$ distinct graphs over $\mathbb Q$, carrying the external-label and edge-order permutation signs separately. The first twelve graphs are identical in the two first-mark choices. For **each** mark the corrected score before subtracting the classical section is

$$
\frac{-92613596121446400}{2(28/3)^4}=-\frac{299015516814300}{49}.
$$

Using the tree order gives $-11362589638943400/49$. Swapping the two children in graph 0 alone changes that score by $7774403437171800/49$ with the source, external states, and classical section fixed. The edge-ordered score is invariant under this swap. [Exact pole evaluations](exact-poles.jsonl) · [Independent contraction](contract.py)

## Subtracting the classical section

Let $M$ be the $703\times703$ matrix of classical pivot images and $y$ the outgoing row of B. The full outgoing rank over $\mathbb Q$ is 703, so $c=yM^{-1}$ subtracts the entire $Q_0$ image. The first correction exists because the intermediate classical map surjects onto all 633 terminal cochains. For this source the bare $Q_2$ would need six flavour slots on three letters with at most four. For the resulting cycle $X$ after the six Ward substitutions, a first correction $X_1$ solves $Q_0X_1=-Q_1X$. The second-page class is $d_2[X]=[Q_2X+Q_1X_1]=[Q_1X_1]$. [Classical rank and Cartan upper bound](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank) · [First correction and Ward return](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential)

Write $s$ for the column of scores of the pivot preimages. Regenerating the three six-cubic pivot recipes and evaluating their eleven distinct terminal graphs over $\mathbb Q$ gives, in zero-based row order and for either first mark,

$$
(s_{700},s_{701},s_{702})=\left(-\frac{598031033628600}{7},0,0\right).
$$

For the unclosed graph score $z$, the order-704 determinant is $\det(M)(z-cs)$. Both complete matrices were eliminated directly. Each entry in the table applies to either first mark. The original archive reported only the left mark at 1013. The rational replay also evaluates the right mark there. [Six-cubic poles](section-poles.jsonl) · [Determinant replay](verify.py)

| Prime | $\det M$ | Tree-order $z$ | Edge-ordered $z$ | Section $cs$ | Archived minor | Edge-ordered minor |
|---:|---:|---:|---:|---:|---:|---:|
| 1009 | 1005 | 663 | 44 | 44 | 551 | 0 |
| 1013 | 734 | 457 | 92 | 92 | 478 | 0 |

The rational section is regular at both primes. Thus the reduced numerator of $z-cs$ is divisible by $1009\cdot1013=1022117$. Its exact value over $\mathbb Q$ remains undetermined.

## Pure glueball controls

On $E_2=H(E_1,d_1)$, where $E_1$ is the classical cohomology, an incoming first differential is excluded by a Ward covector $\lambda$ with $\lambda Q_1=\eta Q_0$ and $\lambda(S^p)\ne0$. The following are previously established controls. The relative complex and the paper's native normalization are used throughout. [Sixth-power Ward proof](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#the-complete-first-page-test) · [Seventh- and eighth-power programs and results](https://github.com/lbotao399-stack/quantum-lifts-results/tree/53462d3102e77eea433207ab1d0fd7e09a8303f6/algorithms/dual_singlets_20261118/f4/j36_source/eight_albert_review)

| Pure power | Level $j$ | Perturbative result |
|---|---:|---|
| $S^5$ | 30 | Survives to $E_\infty$. An incoming $d_1$ would hit $S^6$ after multiplication by $S$. The only possible $d_2$ source has dimension $141-80-61=0$, and the allowed incoming order is at most two. |
| $S^6$ | 36 | Nonzero on $E_2$, with $\lambda(S^6)=-15313094596800$. Its Cartan restriction vanishes by exterior degree six on a rank-four Cartan. An incoming $d_2$ remains open. |
| $S^7$ | 42 | Nonzero on $E_2$, with $\lambda(S^7)=-194893931232000$. The independent Ward test solves 68 coefficients over $\mathbb Q$ and checks 1,117 differentiated-scalar and 18,948 mixed-jet components. An incoming $d_2$ remains open. |
| $S^8$ | 48 | Hit by $d_1$. The homogeneous and augmented Ward minors of orders 91 and 92 are $(156,375)$ modulo 1009 and $(382,29)$ modulo 1013. The one-dimensional classical target is in the image. |

The quantum product theorem gives $d_2([X][S])=[S^7]$ **if** $d_2[X]=[S^6]$, since $S$ is closed. This implication still needs its antecedent. These filtration orders count supercharge corrections. A determination of an anomalous dimension requires the dilatation operator. The identification with untwisted physical local operators uses the quantum comparison in Appendix B. This calculation tests perturbation theory only. [Quantum product and induced derivation](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/f927e3492c1c8b730026041478fcc8ab62741a4a/main.tex#L4769-L4793) · [Appendix B comparison](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/f927e3492c1c8b730026041478fcc8ab62741a4a/main.tex#L5124)

## Reproduction

Use Python 3.10 or later with NumPy. From this directory in the research branch or the standalone bundle, run

```sh
shasum -a 256 -c SHA256SUMS
python3 -B verify.py
python3 -B exact.py --graph-zero
python3 -B exact.py --workers 2 > crossed-rebuilt.jsonl
diff -u exact-poles.jsonl crossed-rebuilt.jsonl
python3 -B section.py --workers 2 > section-rebuilt.jsonl
diff -u section-poles.jsonl section-rebuilt.jsonl
```

The optional provenance check uses Node.js 18 or later and the original ZIP. It verifies the ZIP hash, native tensors, individual archived pole records, the outgoing row at both primes, and the graph construction for each six-cubic pivot. `verify.py --classical DIRECTORY` accepts the two `selected-*.u16.gz` matrices in a different location. [Pinned original archive](https://github.com/lbotao399-stack/quantum-lifts-results/blob/783cdacb1171fb21ef5bdaa155fe14202a8d026d/reviews/f4-j36-eight-albert-d2/f4-j36-eight-Albert-certificate.zip) · [Full independent classical row rebuild](https://github.com/lbotao399-stack/quantum-lifts-results/blob/53462d3102e77eea433207ab1d0fd7e09a8303f6/algorithms/dual_singlets_20261118/f4/j36_source/eight_albert_review/check_outgoing.py)

```sh
python3 -B check_inputs.py /path/to/f4-j36-eight-Albert-certificate.zip
```
