# All-order formal survival of the SU(3) F1 tower

In the same normalized $\mathcal N=4$ holomorphic Taylor-jet complex as the [level-30 proof](../dot-r007-all-order/ALL_ORDER_PROOF.md), the corrected $C=+1$ 8640-term $j=32$ seed has a finite closed lift over $\mathbb Q[[u]]$ and is not $u$-torsion. Consequently the two seed operators and their marked scalar products give at least one non-torsion formal class at every $j=24+2n$, $n\ge3$, in the symmetric flavour $[n,0]$ of the $F_1$ tower. This is a formal perturbative existence result, not a uniqueness result, a physical all-coupling BPS statement, or a nonperturbative lift.

## Exact classical input and quantum syzygy

Use the 47 ordered terms of [dMK et al. Eq. (2.26)](https://arxiv.org/html/2412.08695v2#S2.SS2), with the projected tensors and $C$-even representative made explicit in the [frozen level-32 source](https://github.com/lbotao399-stack/quantum-lifts-results/blob/9361ca1bef07dbb79ea43e5afb2c2beff9640e5f/algorithms/dual_singlets_20261118/taylor_hooks/su3_level32/O32_C_even.json.gz). Its decompressed JSON SHA-256 is $897fbf27bced5562684f9537ed1bd9e1b8a352f7149cc3b99ddb22c89b6c823e$. The factorized replay reconstructs this entire source byte for byte.

Write the highest-weight seed as $\sum_a \alpha_a f_a r_a$ with the printed ordered coefficients. The local relations $Q_0 r_a=R_a$ are the same five already certified for level 30 and two additional $r_{14}$ tensors of types $[1,0]$ and $[0,2]$. The companion generic traceless $3\times3$ coefficient calculation checks their highest components, including an extra mixed-index check for $[0,2]$. The projected equivariant tensor identities follow by flavour lowering. Their $R_a$ are polynomials in the closed $u_2,u_3,v_2,v_3$ graviton components. Replacing these by free graded-commutative component symbols with only the stated linear projections yields
$$
 \sum_a \alpha_a (-1)^{|f_a|} f_a R_a=0
$$
in the free algebra. The ungraded negative control has 716 terms, so the differential Koszul sign is essential. Linearity of the split inclusion into the raw symbols lets the quantum reordering use only a basis of the **closed** projected components.

The same positive-order homotopy argument as for level 30, from [LLC Eqs. (6.13), (6.26) and Theorem 9.1](https://arxiv.org/pdf/2611.24160v1), therefore gives
$$
 S=\sum_a \alpha_a(-1)^{|f_a|} f_a^\star\star R_a^\star=QH,\qquad H=O(u).
$$
All starred objects and $Q$ use the same scheme. We project each $r$ to $C=(-1)^{Y_r+1}$, which leaves the relation $R$ unchanged and puts each term of the reconstructed sum in the $C=+1$ sector, because the total source length is 11. This is just the projection already made in the frozen input.

## Small local repairs

The initial weights $\theta$ of $12_{03},12_{11},12_{22},14_{10},14_{02},14_{13},16_{23}$ are $(3,3,0),(3,2,1),(4,2,0),(3,2,2),(3,3,1),(4,3,0),(5,3,0)$. Loop order one removes one of each $\theta_i$ and adds $(1,1)$ dotted jets. Thus the only $Q_0$-closed first residuals $Q_1r_a-[u]R_a^\star$ that need a check are $12_{11},14_{10},14_{02}$. They are all Cartan-null from their opposite-to-letter-number colour parity.

These are exact integer highest-weight chain dimensions and rational lower ranks of the genuine incoming and stacked outgoing maps, computed by two-prime replay on invariant/HW trace inputs. The outgoing rank is that of $A=(Q_0,\rho)$ with the **same** evaluated columns, not an unrelated Cartan rank.

| Slot $(\theta;n_+,n_-)_{Y}$ | $n$ | lower $\operatorname{rank} I$ | lower $\operatorname{rank} A$ |
|---|---:|---:|---:|
| $(2,1,0;1,1)_4$ | 12 | 4 | 8 |
| $(2,1,1;1,1)_5$ | 32 | 21 | 11 |
| $(2,2,0;1,1)_5$ | 18 | 9 | 9 |
| $(1,0,0;2,2)_4$ | 7 | 3 | 4 |

The inclusions $\mathrm{im}\,I\subseteq\ker A$ and saturation of each $n$ give zero Cartan-null cohomology over $\mathbb Q$. Repair the first order equivariantly to obtain $T_a=r_a+u b_a$, $Q T_a=R_a^\star\bmod u^2$. Only $14_{10}$ has a possible second local residual $e^i=[u^2](QT_{14_{10}}^i-R_{14_{10}}^{\star i})$. Its highest member is $(1,0,0;2,2)_4$.

Set $W=\sum_a\alpha_a f_a^\star\star T_a-H$. It is first-order closed by the strict Leibniz rule. Every first defect has been canceled, so its second residual in the highest component is **only**
$$
 31500\,\epsilon^{abc} v^{1}_a v^{1}_b v^{1}_c\,e^1 .
$$
The sign is positive because the single printed coefficient is $-31500$ and the prefactor contains three odd $v$'s. The [all-page Cartan criterion](../dot-cartan-image/PROOF.md) makes the restriction of this actual second residual zero. The exact Cartan prefactor is nonzero: at $\phi^1=\mathrm{diag}(1,-1,0)$, other scalars zero, $\psi_a=\eta_a\mathrm{diag}(1,-1,0)$, it equals $1008000\,\eta_1\eta_2\eta_3$, with independent odd $\eta_a$. Its odd variables are only the $\beta$ letters from $\psi$, whereas the highest $e^1$ slot contains only scalar and relative $c$ jets, no $\beta$. Hence multiplication is injective on that subalgebra: take a nonzero scalar-polynomial coefficient of a $\beta$-monomial; it is a non-zero-divisor over the free $c$-exterior basis. No claim of injectivity on the entire Cartan superalgebra is needed. It follows that $\rho(e^1)=0$, and by covariance all components restrict to zero. The last rank row provides the local second primitive. No later local terms are allowed by the theta weights; after this repair $W$ is exactly closed. For a homogeneous finite $W$, the rational fan differential gives rational descent exactly as at level 30.

## Split pairing for the two seeds and the tower

Use the block restriction $\phi^1=\mathrm{diag}(X+h\mathbf1_2,-2h)$ and all other letters in the upper two-colour block. The coefficient of $h^4$ in the level-32 seed is literally $-8164800\,i$ times the frozen 117-term $SU(2)$ representative in the unscaled Pauli convention, with separating dual pairing $244944000\,i$. This is independent of the convention used to print a scalar coefficient in the cited paper. The zero-jet projection intertwines the formal differential with $Q_0$, so the dual annihilates every formal boundary and its nonzero pairing excludes $u$-power torsion, just as for the level-30 seed.

For $n\ge5$, use $n-3=2a+3b$ with $a,b\ge0$ and the class of the fixed-bracket marked product $u_2^a\star u_3^b\star O_{30}$ in the same highest scalar flavour. Its leading polynomial is the usual one in the $F_1$ tower. The leading block coefficients of $u_2,u_3$ are $6h^2,-6h^3$; their next powers are absent by $\mathrm{tr}\,X=0$. No product correction is needed for the leading $h^n$ pairing. The seed's highest possible power $h^5$ would multiply a pure $\psi^5$ coefficient. It vanishes: coefficientwise classical closure of a derivative-free pure-$\psi$ polynomial forces every $\partial_{\psi_m}$ to vanish, by the independent $\phi_j^B\phi_k^C$ coefficients of its commutator variation and $[\mathfrak{su}(2),\mathfrak{su}(2)]=\mathfrak{su}(2)$. The preceding $h^4$ coefficient cannot combine with a missing next-to-leading scalar-product power. Thus the $h^n$ coefficient pairs to $6^a(-6)^b$ times the nonzero level-30 dual. A second use of the zero-jet projection proves non-torsion for these representatives. At each level the chosen leading class has $Y=n+7$ and $C=(-1)^n$, opposite to same-length multigravitons. No assertion about other representatives, multiplicities, or arbitrary products is made.

## Replay

See REPLAY.md. The checks reconstruct the classical source and free syzygy, the small invariant rank lower bounds, and the block-separating dual. They do not list the full higher-loop coefficient vectors, compute the exact physical radial operator dictionary, or measure an instanton obstruction.

