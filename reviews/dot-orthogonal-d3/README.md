# Cross-check of the orthogonal \(d_3\) sources and the \(F_4\) correction

Discussion: [issue #1](https://github.com/lbotao399-stack/quantum-lifts-results/issues/1). The independent comparison starts from commit [ee9e064](https://github.com/lbotao399-stack/quantum-lifts-results/commit/ee9e06494a926fb88baeb6fd875cee36762f63a7). This branch is for reproducible cross-checks; it does not modify the default branch.

## \(F_4\), \(j=36\): the previously proposed source is withdrawn

Independent exact replays agree that the 6,417-term \(X\) from the pinned comparison obeys \(Q_0X=O_0=P^6(Na)\) on all 10,401 collected trace words, with the *same* \(O_0\) as the original package. The old exporter omitted the two branch-dependent exterior metric factors \(-4,+2\); restoring them gives opposite contributions of magnitude \(283511008535040000/7\) and exact total zero. Those are not an overall loop rescaling, nor a change to \(Q_2\). The leading class is already zero on \(E_1\). The old \(S^6\) witness and the \(S^7\) consequence obtained from it are withdrawn. This does not prove that no different source can hit \(S^6\).

**Proposal for the next comparison.** Keep the pinned metric, Grassmann order and \(Q_0\) conventions; audit the *target* earlier-page quotient independently of any source score. In the \(j=36\) singlet, compute the complete classical predecessor image \(Q_0\colon W=2,Y=12\longrightarrow W=1,Y=13\) at \((\theta;n)=(222;44)\), or give an independently checked classical separating row for a candidate. Use only the resulting nonzero source classes to search for \(d_1\)-closed classes and evaluate \(d_2\). A new representative of the withdrawn class, a different first correction, or a change in overall loop units does not repair it. Please post any new candidate here with its sign conventions and independent source gate.

## Orthogonal source comparison

These claims concern the fixed **relative holomorphic jet complex** and defining-vector trace convention, not unconditional identification with the untwisted physical local-operator complex. \(s\) denotes the named target after setting \(\lambda(s)=1\). The table states scores for **these specific sources** with the common one-loop prefactor stripped. It does not require distinct sources into a one-dimensional target to give the same coefficient.

| Gauge algebra / target | Leading Gram relation before strict return | Finite first-correction Gram relation | Cubic coefficient |
| --- | --- | --- | ---: |
| \(\mathfrak{so}(11)\), \(j=36,\ S^6\) | \(G_{\mathrm{cycle}}-5G_{\mathrm{path}}+(500/3)G_{\mathrm{alt}}\) | \(0\): all 122 stable Cartan coordinates cancel | \(+3072000\) |
| \(\mathfrak{so}(12)\), \(j=42,\ S^7\), reflection-even | Single promoted Gram13 source from the 3-edge path | \(374400\,\Pi\mathrm{Tail14}-187200\,\Pi\mathrm{Tail146}\), using the endpoint conventions in the corrected full replay | \(-33177600\) |
| \(\mathfrak{so}(13)\), \(j=42,\ S^7\) | \(G_{\mathrm{promoted\,cycle}}+(5/2)G_{\mathrm{promoted\,path}}\) | \(2400\,G_{\mathrm{sat}}\) with the signed saturated Gram14 normalization | \(-921600\) |

For each row, the leading cochain is the strict sixfold \(P\)-return of the symmetry-projected, graded trace-reversal primitive of the displayed Gram relation. To construct the subsequent cochains, subtract the displayed finite relation from the stable first residue, then take its \(Q_0\) primitive by graded adjacent swaps and trace reversal. The second stable residue is odd in total trace length, so the same reversal algorithm gives its \(Q_0\) primitive. Choose signs so \(Q_0x_1=-Q_1x_0\) and \(Q_0x_2=-Q_1x_1\) after specialization. Bare \(Q_{r\ge 2}\) vanish on each of these \(W=1\) blocks. The last number is the actual \(\lambda(Q_1x_2)\), using a finite dual that kills all earlier incoming images. The replay files include finite graded cochain constructions and two classical primitive choices, rather than a dense expansion of every jet monomial.

In particular, the **SO(12) replay uses all 89 nonzero source moments, not just the initial 54 high-factor moments**. The two actual tail vectors and source match on all 121 union-support coordinates after removal of the historical output factor guard. The \(SO(13)\) \(j=42\) example is distinct from any mixed \(j=49\) computation.

Replay bundles and command descriptions are added under this directory. For the independent \(SO(11)\) score \(8640\), please give the exact source vector and units before making a direct comparison to \(3072000\).
