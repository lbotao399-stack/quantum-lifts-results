# Paired Albert covariants in the second obstruction to the sixth glueball

In the four-scalar polarization of the $F_4$ source at $j=36$, consider two separate bilinear defining traces and two open words paired through the $26$ in $\operatorname{End}(26)$. In Taylor type $(4,4,3,3)$, every classical cycle made from this family has a classical primitive. Consequently this family supplies zero incoming second differential to $S^6$, where $S$ is the defining-trace glueball. The full selected source contains 3,720 cochains and its classical split has rank 703. A corrected augmented minor of order 704 anywhere in that source would establish a lift. [Full classical rank](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3b30b8ff242977debd90581c12f4a767d7875c6c/reviews/independent-f4/sixth-glueball-classical-rank) · [Supercharge and corrected-row criterion](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential)

## Reduction to quartics

The second jet with scalar labels $i,j$ has classical variation $[X_i,X_j]$. All first scalar jets $X_i$ are closed. The Albert covariant projects an ordered word onto the traceless Jordan multiplication operators using the trace form. It kills words of length at most three, and adjacent exchanges therefore leave its four-input value unchanged. Call this symmetric value $J$. At thirteen letters the only possible four lengths are $2,2,4,5$. Metric transposition sends the projected five-input word to the negative of its reversal. Reversing by adjacent exchanges gives twice the five-input word as the sum over its ten ordered pairs, replacing each pair by its Lie bracket in $J$ and retaining the other three inputs. These are ordinary colour identities before super-polarization. [Quartic symmetry and native tensors](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/mixed_colours.md)

If the marked jet occurs inside the bracket with a scalar, invariance moves the scalar adjoint action onto the seven other inputs of the paired quartics, three on the quartic containing the mark and four on the opposite quartic. This eliminates the nested case already on source cochains. The spectator traces multiply that identity.

## Kernel and primitives

For the four-dimensional scalar-label space $V$, the ordinary multiplicity to extract is $(4,4,4,2)$. The two bilinear traces contribute $\operatorname{Sym}^2(\operatorname{Sym}^2 V)$ when neither is marked. A marked pair and a scalar bracket both have label space $\bigwedge^2 V$. In opposite quartics the source distinguishes the two copies of $\bigwedge^2 V\otimes\operatorname{Sym}^3 V$. Replacing the mark by its bracket symmetrizes the copies. Within one quartic the same argument symmetrizes the two copies of $\bigwedge^2 V$. When the mark lies in a bilinear, the split is $\operatorname{Tr}_{26}([X_i,X_j]X_k)$, totally alternating in its three labels. The exact Weyl character calculation gives the following multiplicities.

| Marked jet | Source before additional colour identities | Upper bound for split image | Dimension of formal kernel |
|---|---:|---:|---:|
| Opposite quartic | 9 | 7 | 2 |
| Same quartic | 5 | 4 | 1 |
| Bilinear trace | 5 | 1 | 4 |

At the first twelve integer colour assignments of the supplied generator, the same twelve classical split rows give determinants $413$ modulo $1009$ and $948$ modulo $1013$. Hence the combined rational image has dimension exactly $7+4+1=12$. Any quotient by additional native identities can remove only formal kernel vectors. In each of the first two cases the antisymmetrized difference is $Q_0$ of the invariant with both pairs replaced by odd second jets, using the minus sign on the second variation. In the third case the kernel is the kernel of the alternation of the three trace labels, tensored with the remaining factors and projected to the specified diagram. Its marked trace component is closed. The other quartic contains a scalar bracket. Replacing that bracket by its second jet supplies a primitive for the product, with the sign from the odd marked trace. The characteristic-zero Young projection commutes with each split and with these primitive constructions.

## Native quantum check

The nine opposite-quartic sources give the following values of $\eta Q_1$ modulo $1009$, with $\lambda Q_1=\eta Q_0$ in the supplied native intermediate frame.

$$
(281,905,130,182,468,645,827,130,307).
$$

The classical functional in `classical_scores.json` gives the same nine values on their splits. These residues test signs and normalization on individual cochains. The rational kernel statement above follows from the primitives and the minor. The native replay resolves the five-letter identity before moving the ghosts. It uses 2,747,440 weighted colour assignments and checks 25,000 contractions against literal five-letter matrix products. The two-trace intermediate frame and the contracted-pair normalization come from the project's earlier trace calculation. [Frame, ordered action and reduction identities](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/trace_obstruction.md)

## Replay

Use Python 3.10 or later with NumPy and a C++17 compiler. The first command checks every supplied file. The second rebuilds both character bounds and both twelve-by-twelve matrices. The third also rebuilds the native pair contractions and all 706 scalar evaluation points.

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify.py | diff -u verified.json -
python3 -B quantum.py --native --all-points | diff -u quantum_verified.json -
~~~
