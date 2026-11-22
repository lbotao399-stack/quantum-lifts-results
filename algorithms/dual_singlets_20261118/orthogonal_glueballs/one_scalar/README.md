# Earlier-page survival with one scalar

In the full relative holomorphic jet complex, let $S=\operatorname{Tr}_N(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$ use the defining trace of $SO(N)$. With $\phi^1$ the first scalar flavour, the class of $S^{k-1}\operatorname{Tr}_N(\bar\lambda_{\dot+}\phi^1)$ is nonzero on the third page in both rows.

| $N$ | $k$ | $j$ | Finite-colour relation dimensions at orders 0, 1, 2 | Ordered-word minor modulo $1000003$ |
|---:|---:|---:|---|---:|
| 11 | 6 | 35 | 0, 0, 123 | 797069 |
| 13 | 7 | 41 | 0, 0, 17 | 376015 |

The terminal weight has $k$ positive dotted gauginos, $k-1$ negative dotted gauginos and one scalar. In the stable classical quotient, graded commutation and orthogonal transposition leave only products of bilinears at this weight. An even trace can contain each odd letter once at most and there is only one scalar. The displayed product is the unique surviving monomial. The zero finite-colour relation dimension in the first column of each triple preserves its class. All positive-order supercharges annihilate this representative by flavour degree.

Let $\lambda$ extract the coefficient of this product. In the stable trace algebra the verifier constructs $\mathrm{first}$ and $\mathrm{second}$ from graded commutator primitives and checks

$$
\lambda Q_1=\mathrm{first}\,Q_0,\qquad
-\mathrm{first}\,Q_1=\mathrm{second}\,Q_0.
$$

Inverse reconnection and classical trace splitting enumerate the potentially nonzero columns. The two checks cover 55 and 2,178 columns respectively, and the two rows have 17 and 383 nonzero entries. A bare $Q_2$ reaching the target would take a source of length $2k+1$ with $2k+2$ Taylor slots. It therefore has one second jet and otherwise first jets. Three active letters supply at most four of the six required flavour slots. [Transposed equations](../../orthogonal_trace_dual.py)

At the first dual weight there are zero finite-colour identities in the full weight space. At the second dual weight flavour equivariance selects the highest weight $(3,2,2)$. The invariant-tensor character gives the dimension in the table at fixed dotted degrees $(k-2,k-3)$. The polarized Pfaffians use two scalar determinant projectors and one derivative of flavour 1. Each exceeds degree $N$ and hence vanishes after specialization. The displayed minors prove rational independence, so these rows span the entire required relation space. The centre element $-1$ acts trivially on adjoints and interchanges the two components of $O(N)$ for both odd $N$. [Characters and polarization](../../orthogonal_one_scalar.py)

On the stable commuting quotient, let $C_i$ extract coefficient $i$ in the zero-based distinct-plane basis stored in the replay. On every finite relation the second row equals

$$
\begin{aligned}
SO(11):\quad &-\tfrac34 C_1-\tfrac5{12} C_2+\tfrac38 C_5-\tfrac12 C_6,\\
SO(13):\quad &-\tfrac1{64} C_0 .
\end{aligned}
$$

Subtracting that functional preserves the second transposed equation because the $C_i$ kill classical boundaries. The difference descends through every finite-colour relation in the target flavour representation. This excludes incoming $d_1$ and $d_2$ on the two specified lines. Later incoming differentials are open.

With the repository's Python trace dependencies, run from this directory

~~~sh
python3 verify.py --rebuild --workers 3
~~~

The replay regenerates both adjoints, the ordered-word minors, the indicated commuting-plane coefficients and the rational second scores for all 140 independent Pfaffians. [Verifier](verify.py) · [SO(11) inputs](one-scalar-relations-11-6.json.gz) · [SO(13) inputs](one-scalar-relations-13-7.json.gz)
