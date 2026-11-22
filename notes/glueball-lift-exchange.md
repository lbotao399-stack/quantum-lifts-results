# Glueball lift correspondence

Dot and Codex are comparing the calculations in [GitHub issue #1](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1).

## November 21, 2026 · Results of the exchange

The calculations use the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM. Here $S$ is the glueball, $Q_0$ the classical supercharge and $d_r$ the differential on page $r$. Physical identification of the complex with untwisted local operators uses the comparison hypothesis of the all-loop manuscript.

| Target | Dot's source | Codex cross-check |
| --- | --- | --- |
| $F_4$, $j=36$, $S^6$ | Dot withdrew $O_0=P^6(Na)$ and the $S^7$ consequence derived from it. | $Q_0X=O_0$ in all 10,401 trace words. Restoring the two exterior metric factors gives zero on all four archived recipes. A different source for $S^6$ is open. |
| $SO(11)$, $j=36$, $S^6$ | Gram cycle minus five times the path, plus $500/3$ times the alternate path. | All 122 first-correction coordinates cancel. The third coefficient is $1572864000/512=3072000$. A separate source gives $8640$. |
| $SO(12)$, $j=42$, $S^7$ | Reflection-even promoted three-edge path. | The 89 nonzero source coordinates equal $374400$ times tail 14 minus $187200$ times tail 146 on the full union of 121 coordinates. The third coefficient is $-33177600$. |
| $SO(13)$, $j=42$, $S^7$ | Promoted cycle plus $5/2$ times the promoted path. | The 122-coordinate residual equals $2400$ times the finite Gram relation. The third coefficient is $-5529600+(5/2)(1843200)=-921600$. Thus $[S^7]$ is zero on page four. |

Dot reports regenerating all three signed source circuits. Codex verified the three archive hashes, all 90 member hashes, the three fast replays and the earlier-page certificates for eleven, twelve and thirteen colours. A separately written rational verifier reads both primitive choices in the archived moment tables and tests every coordinate and normalized third coefficient. This comparison uses Dot's upstream contraction tables. [Dot's withdrawal](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1#issuecomment-7074126654) · [Dot's pinned sources](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md) · [Independent coordinate verifier](../algorithms/dual_singlets_20261118/orthogonal_glueballs/verify_dot_coordinates.py)

For the independent $SO(11)$ source, `PfaffianCoefficient(3,0)` specifies the leading degree-15 relation. The leading cochain is its graded reversal primitive. Its first correction uses $135$ times determinant row $1$ minus $90$ times each of rows $9$ and $18$ (zero-based). In the archived one-loop kernel the common prefactor is $K=1$. The target uses the defining trace $S=\operatorname{Tr}_{11}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$. The pairing is $-17280+135(288)-90(192)-90(-48)=8640$. [Source, corrections and target normalization](../algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36/README.md) · [Coefficient extraction](../algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36/replay_forward.py)

[The classical $F_4$ preimage and four corrected pairings](../algorithms/dual_singlets_20261118/f4/j36_source/README.md) · [Earlier pages for thirteen colours](../algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36/README.md) · [Earlier pages for twelve colours](../algorithms/dual_singlets_20261118/orthogonal_glueballs/even_second/README.md)
