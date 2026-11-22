# Quadratic scalar targets in fourteen, fifteen and sixteen colours

Use $S=\operatorname{Tr}_{N}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$ and $A=\phi^1$ in the full relative holomorphic jet complex. At $p=N-6$ the target with $p$ gauginos of each dotted sign and two copies of $A$ has the following ordered classical trace coordinates,

$$
S^p\operatorname{Tr}(A^2),\qquad
S^{p-1}\operatorname{Tr}(\bar\lambda_{\dot-}A)\operatorname{Tr}(\bar\lambda_{\dot+}A),\qquad
S^{p-1}\operatorname{Tr}(\bar\lambda_{\dot+}\bar\lambda_{\dot-}A^2).
$$

The finite-colour relations at tree level have zero projection in this sector. The coordinates are dotted singlets, highest in the symmetric two-scalar flavour representation. In even dimension the calculation is reflection even. Each group supplies a rank-one incoming $d_1$ and two independent classes surviving the incoming $d_2$. Positive-order outgoing terms vanish on these representatives by their single scalar flavour. [Trace coordinates, grading and tree primitive](../scalar_multipliers/README.md)

| Gauge group | $p$ | $j=6p+4$ | Ordered first image |
|---|---:|---:|---|
| $SO(14)$ | 8 | 52 | $(160,-320,-3200)$ |
| $SO(15)$ | 9 | 58 | $(-400,800,8800)$ |
| $SO(16)$ | 10 | 64 | $(960,-1920,-23040)$ |

At the first transpose the only determinant type consists of identical open paths with both gauginos, three positive and three negative single gaugino axes, and singleton scalars. The two exterior triples of auxiliary axes each transform by the determinant of the adjoint of $SL(2)$, hence trivially. The two residual scalar axes contribute one invariant symmetric quadratic form. Evaluating them on the two Cartan axes or on the raising and lowering axes gives rows in the ratio $2:1$ in all three groups. The finite relation has odd trace length and admits a stable classical primitive by graded trace reversal. Each nonzero row therefore represents an actual first image.

In each case the first image is proportional to $(1,-2,-(2p+4))$. Its annihilator is spanned by $(2,1,0)$ and $(2p+4,0,1)$. The corresponding second transpose equations cover 4,678 and 15,928 native columns. A bare second-order term requires six flavour slots on three inputs, while the predecessor has only one second jet and supplies at most four. The remaining contribution is the iterated first-order pairing. [Construction of the transpose and path reduction](README.md)

The finite-colour determinant identities at the second step have 46 path types. Reduction by the auxiliary $SL(2)$ action leaves 208 representatives of 1,848 weight-zero components per group. On these 208 rows the selected twenty commuting-plane columns have rank 20. Adjoining both second scores leaves rank 20 in each group. The replay solves over the rationals and checks both solutions against every row. Subtracting these commuting functionals makes both transposes descend to finite colour. [All three rational solutions and minors](adjacent_results.json) · [Rows for fourteen](rows-14-8.json.gz) · [fifteen](rows-15-9.json.gz) · [sixteen colours](rows-16-10.json.gz)

The original rows combine identical even paths and use a trace with a marked starting point. The verifier defaults to checking the first two independent rows in each group and the full rational residuals. The native option evaluates these rows with individually labelled paths and the original trace logarithm. Both returned the same receipt. The rebuild option regenerates every row.

~~~sh
python3 adjacent_colours.py | diff -u adjacent_results.json -
python3 adjacent_colours.py --native | diff -u adjacent_results.json -
~~~
