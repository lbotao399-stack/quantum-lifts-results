# Rational BMN source bounds for third-order tests in $SU(3)$

At BMN level $j=2(q_1+q_2+q_3)$, let $Y$ be the number of letters and $C=\pm1$ the colour charge-conjugation eigenvalue. Each BMN graviton and every product of gravitons has $C=(-1)^Y$. The public Fortuity calculation gave respective lower bounds 4, 8 and 11 on the number of opposite-parity highest weights in the three columns below, conditional on its diagonal multiplicities. Here the bound follows from explicit rational products and the two separate Euler characteristics. [Public index data](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/kernel_parity_j40.json) · [Parity lemma and scope](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex)

| Level | Flavour weight | Euler characteristics at $C=+,-$ with sign $(-1)^Y$ | Independent products at $C=+,-$ | Opposite-parity lower bounds at $C=+,-$ |
|---:|---|---|---|---|
| 36 | $(8,5,5)$ | $(0,-2)$ | $(2,4)$ | $(2,2)$ |
| 38 | $(9,5,5)$ | $(4,-2)$ | $(7,7)$ | $(3,5)$ |
| 40 | $(10,5,5)$ | $(6,-7)$ | $(12,12)$ | $(6,5)$ |

For $C=+$, odd-length cohomology equals even-length cohomology minus its Euler characteristic. For $C=-$, even-length cohomology equals odd-length cohomology plus its Euler characteristic. At level 36 these inequalities give $2-0=2$ and $4+(-2)=2$. At levels 38 and 40 they give $7-4=3$, $7-2=5$, $12-6=6$ and $12-7=5$. The rows refer to distinct charges. Some of the classes can be descendants or products involving lower fortuitous classes.

The single-trace factors are descendants of the quadratic and cubic scalar traces. For each chosen factor the program checks classical closure on oriented trace words, and checks its charge conjugate coefficient by coefficient. The highest-weight projection of their product commutes with both $Q_0$ and Cartan restriction. Classical boundaries restrict to zero. On the ordinary traceless Cartan the diagonal entries of each letter are $x,y,-x-y$, with exterior variables for fermionic letters. The signed integer product determinants at level 36 are $-8,-6552,3$ at lengths $10,11,12$, of orders $1,4,1$. At level 40 they are $32,-49152,10030613004288,2305926750528,40$ at lengths $10,11,12,13,14$, of orders $1,3,9,9,2$. All denominators and selected monomials are saved. [Level 36](replay_weight_8_5_5.json) · [Level 38](replay_weight_9_5_5.json) · [Level 40](replay_weight_10_5_5.json)

The first calculation uses the cyclic Cartan basis with diagonal entries $(1,z,z^2)$ and $(1,z^2,z)$, where $z$ is a primitive cube root of unity. A product of these diagonal matrices has trace three when the two multiplicities agree modulo three and zero otherwise. A finite-field flavour quotient selects the trial products. The saved determinants are obtained by exact highest-weight projection. The independent calculation uses $x,y,-x-y$ and reconstructs the cochain dimensions by a Weyl-weight subtraction at each length. [Algorithm](../../su3_bmn_parity_products.py) · [Rational replay](../../su3_bmn_parity_replay.py) · [A8 accounting](jobs.json)

## Four inherited classes at level 36

At weight $(8,5,5)$ the two colour parities have rational classical dimensions $(0,0)$ at lengths $8,9,17,18$. At length 10 the even and odd dimensions are

$$
4287-704-3582=1,\qquad 4537-849-3687=1 .
$$

The modular minors bound the rational dimensions from above. A quadratic/cubic trace product and the triply raised level-30 operator give the matching lower bounds. The latter pairs to $20736$ on the $SU(2)$ block. At length 12 the flavour projection of the level-30 operator times the traceless scalar–fermion bilinear pairs to $31104$ in a different block. Both carry odd charge conjugation. The level-30 input has 6,139 terms. Its 19,990 formal classical-variation terms cancel at three colours in all 14 letter contents, with 456,667,072 integer contributions. [Source and image bounds](source_lengths.json) · [Classical closure](descendants/classical.json) · [Odd descendants and boundary functionals](descendants/descendants.json) · [Published source and construction](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/dmk_o330_hw.meta.json)

At length 11 there are two independent even descendants. In the published level-32 operator, four abelian scalars on restriction to $SU(2)\times U(1)$ multiply $-1020600$ times the $SU(2)$ threshold class. In the active block the central weight has rational dimension $153-47-105=1$, with threshold pairing $-20$. The neighbouring root weight has dimension $96-28-68=0$. In the highest-weight projection of the twice-raised operator onto $[3,0]$, only the $1111$ parent component contributes to the coefficient of $(\hat\phi^1)^4$ in cohomology. The other two components have active weights obtained by permuting $(5,4,3)$. For the retained component, applying $Q^3_+$ then $Q^2_+$ to the threshold pairs to $160$. [Published quartic scalar coefficient](https://arxiv.org/html/2412.08695v2#S2.SS1) · [Integer and rational block calculations](quartic.json)

For the once-raised level-34 operator the same pure scalar coefficient vanishes identically. Its coefficient with three abelian scalars and one abelian fermion pairs to $400982400$. Taking these two functionals as rows and the level-32 and level-34 descendants as columns gives

$$
(-1020600)(160)(400982400)=-65478821990400000
$$

for the determinant. The different lengths separate the two odd classes. The two even and the two odd classes account for the *lower bound* of four at this weight. A further primitive or any higher obstruction requires a separate calculation. [Level-34 calculation](descendants/level34.json) · [Quartic scalar replay](../../su3_level36_scalar_channel.py) · [Level-30 and level-34 replay](../../su3_bmn_descendants.py)

A8 job 5683898 gives outgoing eleven-letter image minors of orders 8,977 and 8,520 at prime 2,147,483,647. The two independent even descendants and the four independent odd Cartan products attain the rational upper bounds

$$
12561-3582-8977=2,\qquad 12211-3687-8520=4.
$$

These six classes exhaust the eleven-letter BMN cohomology of weight $(8,5,5)$. [Image minor coordinates](length11_image/image_11.json) · [A8 receipt](length11_image/receipt.json) · [Rational bounds by length](source_lengths.json)

Under filtered superconformal equivariance the third-order Taylor diagrams have rows $(6,6,1^3)$, $(6,6,1^4)$ and $(6,6,1^5)$ in increasing order of level. Here $1^k$ means $k$ rows of one box. Determining a nonzero third obstruction still requires a representative surviving the two lower orders and a nonzero pairing after quotienting both lower correction freedoms.

From the repository root, the calculation at level 36 can be replayed with

```sh
d=algorithms/dual_singlets_20261118
python3 "$d/su3_bmn_parity_replay.py" "$d/taylor_hooks/su3_bmn_sources/weight_8_5_5.json" --out /tmp/su3-855-replay.json
cmp /tmp/su3-855-replay.json "$d/taylor_hooks/su3_bmn_sources/replay_weight_8_5_5.json"
```
