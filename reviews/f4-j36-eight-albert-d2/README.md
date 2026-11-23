# F4, j=36: eight-Albert d₂ candidate

> **还未交叉验证 / Pending independent cross-verification.**
> This is a reproducibility packet for review, not a result cleared for citation. It does not reinstate the earlier withdrawn, classically exact source.

The proposed incoming differential concerns \(S^6\) in the undeformed four-dimensional \(\mathcal N=4\) relative holomorphic jet complex. The new source candidate is the crossed eight-Albert tensor B with its classical outgoing image subtracted using a specified 703-row physical section. The packet records the following arithmetic in that convention:

| Prime | Augmented minor, size 704 |
| ---: | ---: |
| 1009 | 551 |
| 1013 | 478 |

[Download the certificate ZIP](f4-j36-eight-Albert-certificate.zip). Its SHA-256 is
`b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972`.

From the extracted directory, `sh verify.sh` checks the sealed receipt and the two minors. `sh full_replay.sh` regenerates the native contractions instead of reading the saved contraction logs. `sh full_replay.sh 3` also regenerates the candidate's complete outgoing evaluation row at both primes, with three workers.

The [classical rank proof](https://github.com/lbotao399-stack/quantum-lifts-results/blob/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank/README.md) and the [earlier-page and first-correction proof](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential) are pinned separately. The main cross-check is that the actual physical contraction and the corrected right section use the same signed normalization. A nonzero minor in an inconsistent normalization is not a certificate.

The interpretation for untwisted physical local operators requires the manuscript's Appendix B comparison. The \(S^7\) product corollary depends on the \(S^6\) result and the independent \(S^7\) earlier-page gate. This packet makes no nonperturbative claim.
