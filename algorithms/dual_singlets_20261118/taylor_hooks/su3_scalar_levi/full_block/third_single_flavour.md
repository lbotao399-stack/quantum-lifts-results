# Level-36 third targets on the two-colour Levi

For the $SU(3)$ BMN weight $(8,5,5)$, filtered Taylor equivariance puts the third supercharge obstruction in the diagram $(6,6,1,1,1)$. The three tests concern source lengths $10,11,12$ with colour charge-conjugation parities $-,+,-$, respectively. The target length is two less than the source length. On the full $\mathfrak{sl}_2\oplus\mathfrak{gl}_1$ block the multiplicity of this diagram in classical cohomology is zero in each case.

Retain the two even Taylor directions and one odd flavour direction. In the Schur character the boxes filled by that one odd direction form a vertical strip. Removing it from $(6,6,1,1,1)$ must leave at most two rows. The three choices give the components

$$
((6,6),3),\qquad ((6,5),4),\qquad ((5,5),5).
$$

The pair in each entry gives the two even degrees and the final integer gives the odd degree. Each has multiplicity one by the one-variable skew Schur rule. Restriction in colour commutes with the Taylor substitutions. A copy of the specified diagram in the restricted cohomology would therefore contribute to every column below.

| Source length | Target parity | Cochains at $((6,6),3)$ | Cochains at $((6,5),4)$ | Cochains at $((5,5),5)$ | Rational cohomology in the three components |
|---:|:---:|---:|---:|---:|:---:|
| 10 | $-$ | 32064 | 59470 | 20009 | $0,0,0$ |
| 11 | $+$ | 21825 | 51622 | 23972 | $0,0,0$ |
| 12 | $-$ | 7980 | 24583 | 15634 | $0,0,0$ |

The abelian spectator monomials have zero classical variation. For each active $SU(2)$ block with positive branching multiplicity, the adjacent minors at prime 2147483647 saturate its characteristic-zero cochain dimension. The rational rank lower bounds sum to that dimension and $Q_0^2=0$. An independent two-colour character gives the nine cochain dimensions displayed above. A8 job 5687955 supplied 249 previously missing sector records in 124 seconds with 24 allocated CPUs. The three earlier physical-weight Cartan-kernel upper bounds $4,7,3$ thus reduce to $0,0,0$ for the selected diagram on this block. The unrestricted three-colour target remains open. [Bounds](third_single_flavour_bounds.json) · [Additional minors](third_single_flavour_minors.json.gz) · [Slurm receipt](third_single_flavour_job.json)

In the environment of the other Levi calculations, run from the repository root

```sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_scalar_levi/full_block"
python3 "$d/su3_third_levi_weights.py" --records "$p/minors_2147483647.json.gz" --records "$p/third_remaining_minors.json.gz" --records "$p/third_single_flavour_minors.json.gz" --out /tmp/su3-third-single-flavour.json
cmp /tmp/su3-third-single-flavour.json "$p/third_single_flavour_bounds.json"
```
