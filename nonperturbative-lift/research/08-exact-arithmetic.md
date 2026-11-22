<a id="np-section-06"></a>

## 8. Exact arithmetic and the scope of the results

The finite-field census of the source manuscript supplies candidate classes in additional sectors.

Rational certification combines explicit chain dimensions, nonzero minors, and full polynomial identities to determine characteristic-zero cohomology.

> **Rank calculation.**
>
> For adjacent $N_b$ sectors, the source two-step differential is
>
> $$
> D=\begin{pmatrix}Q_0&0\\Q_1&Q_0\end{pmatrix}.
> $$
>
> Project $\operatorname{im}D$ onto its first summand.
>
> The kernel is $\operatorname{im}Q_0+Q_1\ker Q_0$ in the second summand.
>
> Consequently,
>
> $$
> \operatorname{rank}d_1
> =\operatorname{rank}D
> -\operatorname{rank}(Q_0\text{ on first summand})
> -\operatorname{rank}(Q_0\text{ on second summand}).
> $$
>
> A nonzero minor modulo a prime at which the denominators are invertible gives a lower bound on the rational rank.
>
> When the exact character dimension and $Q_0^2=0$ or $d_1^2=0$ make the upper and lower bounds coincide, they determine the rank and cohomology dimension exactly.
>
> Explicit primitives are solved over the rational field and substituted into every unprojected coefficient.
>
> $\square$

| Gauge group and sector | Exact result |
|---|---|
| SU(2), leading two-ghost scalars, all $k$ | Dimension $1$ on $E_2$ for even $k$, and $0$ for odd $k$ |
| SU(2), $C_2$ and $C_4$ product | Full rational identity; nilpotency index $3$ for $[C_2]$ |
| SU(2), first partner and weight $(7,7)$ partner class | Full rational cocycles; one-dimensional specified scalar blocks |
| SU(2), weight $(9,9)$, leading $(5,5)$ | Dimension $1$ on $E_\infty$; explicit representative through the third correction |
| SU(3), leading two-ghost scalars, all $k$ | $E_2=0$ for $k\equiv2\pmod3$, and dimension $1$ otherwise |
| SU(3), $C_1$ partner | Full rational cocycle; vanishing $C_0$ product and $C_0^2$ indeterminacy image |
| SU(3), $C_1^2,C_3,C_4$ | Nonzero full classes; one-dimensional specified homogeneous scalar blocks for $C_3,C_4$ |
| SU(3), leading $k=6$ | Dimension $1$ on $E_3$ |
| SU(4), leading tree scalars, all $k$ | Rank-six free module over the Casimir ring |
| SU(4), Cartan $Q_1$ kernel, all $k$ | Zero for $k\equiv3\pmod4$, and dimension $1$ otherwise |
| SU(4), $k=2,4,5$ | Full cocycle at $k=2$; exact first lifts at $k=4,5$ |
| SU(4), $k=6$ | The first-lift primitive remains open; the all-degree full-$E_2$ theorem is conditional on this step |

The saved rational engine uses $Q_r(\mathrm{physical})=(64\pi^2)^{-r}Q_r(\mathrm{stored})$.

The coefficient arrays retain the stored normalization; $\hbar$ follows the formal loop bookkeeping of each source.

<a id="np-a84ccb2fe764f5cd-1"></a>

### The finite-field census and higher product conjectures

The source census uses the prime $p=2^{31}-1$.

Each tuple denotes $(N_b,N_c)$, with multiplicity one.

Sectors established over the rationals use the preceding proofs.

The remaining entries record modular evidence.

| Weight | SU(2) census | SU(3) census |
|---|---|---|
| $(0,0)$ | $(0,0)$ | $(0,0)$ |
| $(1,1)$ | $(0,2)$ | $(0,2)$ |
| $(2,2)$ | Empty | $(1,2),(0,4)$ |
| $(3,3)$ | $(2,2)$ | Empty |
| $(4,4)$ | Empty | $(3,2),(2,4)$ |
| $(5,5)$ | $(4,2),(1,5)$ | $(4,2),(3,4),(0,7)$ |
| $(6,6)$ | $(2,4)$ | Uncomputed in the source table |
| $(7,7)$ | $(6,2),(3,5)$ | Uncomputed in the source table |
| $(8,8)$ | $(4,4),(1,7)$ | Uncomputed in the source table |
| $(9,9)$, partial census for $N_b\ge4$ | $(5,5)$ | Uncomputed in the source table |

For $SU(2)$, the weight-$(8,8)$ census is complete, while the weight-$(9,9)$ census is partial.

In the leading-component normalization $C_0a_{2i+2j}=a_{2i}a_{2j}$, the source tower conjecture uses $D_0=-C_0$ and $D_j=(2j-1)C_{2j}$:

$$
D_iD_j=D_0D_{i+j},\qquad i,j\ge1.
$$

Here $i,j$ label tower elements.

The completed $i=j=1$ relation is the full rational product above; full products in higher degrees test the remaining cases.

The source also conjectures the Cartan partner ratio $-2(2j+1)(2j+3)$.

The first three saved modular/partial-chain ratios are $-30,-70,-126$.

These conjectures organize the computed products and partner pattern into a common statement.
