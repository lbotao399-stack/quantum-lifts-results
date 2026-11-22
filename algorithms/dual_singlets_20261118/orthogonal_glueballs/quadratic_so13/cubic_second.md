# Three cubic scalar targets in SO(14) and SO(15)

In the defining-trace sector of the full relative holomorphic jet complex, the incoming second supercharge differential vanishes on the terminal weights $(6,5;3,0,0)$ in $SO(14)$ and $(8,7;3,0,0)$ in $SO(15)$. The first two entries count the two dotted gaugino components. The last three count the three scalar flavours. Each target contains three copies of $\phi^1$. The levels are $j=39,51$. In even orthogonal colour the trace calculation is in the reflection-even sector. The supercharge commutes with the reflection, and a pair of colour epsilon tensors reduces to deltas. [Classical primitive and transpose](../one_scalar/README.md) · [Orthogonal second fundamental theorem](https://arxiv.org/abs/1102.3221)

Write $p=5,7$ for the number of negative dotted gauginos in the respective cases. There are $p+1$ positive gauginos and the terminal trace length is $2p+4$. The graded-commuting quotient has three trace coordinates. The path enumeration at the classical and first-correction weights is empty in both cases, so the three coordinates descend through the first page. On the second source relation there are $p-1$ positive gauginos, $p-2$ negative gauginos and nine scalars. The flavour projection is the square of the scalar-marker determinant followed by three derivatives in the first flavour. As in the quadratic calculation, an incoming cochain has one second jet. Three inputs have at most four flavour slots altogether, whereas a bare $Q_2$ needs six. The remaining second pairing is computed by two applications of the one-loop kernel with the tree primitive. [Trace sorting, finite relations and flavour projection](README.md)

| Group | Finite determinant path types | Auxiliary weight-zero components | Independent rows | Rank of commuting coordinates | Rank after adjoining the three second pairings |
| --- | ---: | ---: | ---: | ---: | ---: |
| $SO(14)$ | 2 | 12 | 2 | 1 | 1 |
| $SO(15)$ | 16 | 416 | 37 | 7 | 7 |

In fifteen colours the uncorrected covectors are nonzero on 35, 36 and 34 of the 37 rows. The minor on commuting sites $0,1,2,10,11,13,14$ is $-8491923761101058238068343296522541072384000$. Solving on its seven rows gives one functional per target coordinate. Subtracting each functional kills every row over the rationals. In fourteen colours the minor is $-2688$ and all three uncorrected pairings already vanish. Positive-order outgoing actions on these first-jet representatives vanish by their single scalar flavour. Each of the six target coordinates thus survives to the third page. [Independent rows](cubic_second.json.gz) · [Exact coefficient matrices](cubic_second_results.json)

The original construction uses rooted cyclic traces and merges identical even gaugino-pair paths. The replay uses the original trace logarithm with each labelled path and independently regenerates all three second scores and the selected commuting coordinates on all 39 rows. The first transpose checks 60, 345 and 74 active native columns. The second checks 2392, 19868 and 3567. From this directory with the orthogonal trace runtime dependencies installed, run

~~~sh
python3 cubic_second.py --workers 2 | diff -u cubic_second_results.json -
~~~

The optional `--all-sites` also recomputes every commuting column present in the archived frames, including zero columns. [Verifier](cubic_second.py)
