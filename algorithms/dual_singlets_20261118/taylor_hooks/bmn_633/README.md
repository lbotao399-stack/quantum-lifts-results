# A graviton product in the rank-three Cartan kernel

At BMN level $j=24$, flavour highest weight $(6,3,3)$ and letter length ten, the unique $SO(7)$ classical class that vanishes on the Cartan is a linear combination of graviton products. The same polynomial is a classical boundary for $USp(6)$.

## Rational dimensions

The classical supercharge $Q_0$ raises length by one. Character dimensions at lengths eight through twelve and certified lower bounds on its four ranks are

| Group | Cochain dimensions | Ranks out of lengths 8, 9, 10, 11 |
|---|---|---|
| $SO(7)$ | 929, 1863, 1837, 772, 107 | 712, 1151, 665, 107 |
| $USp(6)$ | 929, 1861, 1832, 770, 107 | 712, 1149, 663, 107 |

The two adjacent rank sums attain the length-nine and length-eleven dimensions in each group. Since $Q_0^2=0$, all four bounds are exact over $\mathbb Q$. The length-ten cohomology dimensions are $1837-1151-665=21$ and $1832-1149-663=20$. [Character calculation and eight minors](resolved.json)

The Weyl projection on the three signed Cartan planes gives $732-583-453+230+262-168=20$. Twenty highest-weight combinations of gravitons restrict to independent polynomials there. Their integer minor has determinant $3152750479073944222118163439163729086669528693714452480$. Every classical boundary restricts to zero, so the two kernels on cohomology have dimensions one and zero. [Product recipes](products.json.gz) · [Cartan minor](replay.json)

## The polynomial

Write $P_0$ for the polynomial reconstructed from the seed recipe by subtracting the unique linear combination of those twenty product vectors with the same Cartan coordinates. Its 1101 signed graviton monomials expand to 1212 trace monomials with common denominator 540. The three raising operators and $Q_0$ annihilate the polynomial before finite-colour identities. Its rank-three Cartan expansion vanishes. The rank-four expansion has 370 nonzero orbit coefficients, including $6144/5$. [Seed](product_recipe.json.gz) · [Coefficients and checks](product.json) · [Expanded integer numerator in the native trace alphabet](product.txt.gz)

At prime 1000033, evaluation on a full $SO(7)$ target frame gives

$$
\operatorname{rank}(Q_0)=1151,\qquad
\operatorname{rank}(Q_0\mid 540P_0)=1152,\qquad
\det(\text{augmented minor})=840083\pmod{1000033}.
$$

The first rank is also exact over $\mathbb Q$. Thus $[P_0]$ spans the $SO(7)$ Cartan kernel. All 21 classes are spanned by products and have perturbatively closed completions by the ring theorem. Injectivity of the $USp(6)$ Cartan restriction gives $[P_0]=0$ there. The comparison alone supplies no instanton matrix element. [Augmented minor and probe coordinates](product_augmented.json) · [Perturbative product theorem](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex)

An independent signed-minor construction uses 665 classical equations and 20 Cartan equations on 686 trace generators. Its order-685 minor has residues 2052382176 modulo 2147483647 and 163439 modulo 1000033. The earlier 2570-term descendant has 295 nonzero Cartan coefficients. Subtracting its product projection gives augmented rank 1151 at the tested prime. [Integer generators and probes](rational.json.gz) · [Two-prime replay](frozen.json) · [Descendant](descendant.json)

For each of eleven lower $SO(7)$ sectors checked by A6 jobs 7657344 and 7657349, the incoming-rank lower bound plus the rank lower bound for the stacked outgoing and Cartan maps equals the full cochain dimension. Their rational Cartan kernels on cohomology vanish. [Sector ranks and file hashes](parents/receipts.json)

## Replay

From the repository root, use Python with NumPy, SciPy and python-flint, and a C++17 compiler.

```sh
d=algorithms/dual_singlets_20261118
t="$d/taylor_hooks"
b="$t/bmn_633"
export PYTHONPATH="$PWD/$d:$PWD/algorithms/su3_j34_20261118${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p /tmp/bmn633-replay
python3 "$d/so7_bmn_633.py" --folder "$b" --out /tmp/bmn633-replay/dimensions.json
python3 "$b/product_replay.py" --export /tmp/bmn633-replay/product.txt --rank "$b/product_augmented.json" --out /tmp/bmn633-replay/product.json
c++ -O3 -std=c++17 -pthread -I "$t/native" "$t/su2_levi_image.cpp" -o /tmp/bmn633-replay/image
PAR_THREADS=4 /tmp/bmn633-replay/image /tmp/bmn633-replay/product.txt 6 3 3 0 0 10 1000033 778117 SO7 > /tmp/bmn633-replay/minor.json
cmp /tmp/bmn633-replay/minor.json "$b/product_augmented.json"
c++ -O3 -std=c++17 -pthread -I "$t/native" "$b/frozen_rows.cpp" -o /tmp/bmn633-replay/frozen
python3 "$b/frozen_replay.py" --recipe "$b/rational.json.gz" --check "$b/cofactor.json" --engine /tmp/bmn633-replay/frozen --out /tmp/bmn633-replay/frozen.json
```
