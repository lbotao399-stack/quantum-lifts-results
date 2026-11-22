# The one-loop bracket of trace Casimirs is a tree boundary

The fields and one-loop normalization are those of the preserved product engine. For integers $m,n\ge2$ and the source rotation generators $J_+,J_-$, which raise and lower the two holomorphic derivative labels, respectively,

$$
\begin{aligned}
&Q_1\bigl(\operatorname{tr}(b^m)\operatorname{tr}(b^n)\bigr)
-\operatorname{tr}(b^m)Q_1\operatorname{tr}(b^n)
-\operatorname{tr}(b^n)Q_1\operatorname{tr}(b^m)\\
&\qquad=4mn\,Q_0\left[
\left(1-\frac12J_-J_+\right)
\operatorname{tr}\bigl((\partial_1b^{m-1})b(\partial_2b^{n-1})\bigr)
\right].
\end{aligned}
$$

The parentheses $\partial_\alpha b^r$ mean the derivative of the entire matrix power. The scalar projector here antisymmetrizes its two first derivative labels. The identity holds for the source $\mathfrak{sl}_N$ matrix theory; no finite-degree Cayley–Hamilton relation is used in its proof.

> **Proof.** Use the engine's two independent even matrix probes $p,q$, held fixed when computing the two-input symbol. The one-loop two-$b$ symbol is
>
> $$
> \begin{aligned}
> Q_1\bigl(\operatorname{tr}(pb)\operatorname{tr}(qb)\bigr)
> =-\operatorname{tr}\Bigl(&[p,\partial_1c][q,\partial_2b]
> +[q,\partial_1c][p,\partial_2b]\\
> &-[p,\partial_2c][q,\partial_1b]
> -[q,\partial_2c][p,\partial_1b]\Bigr).
> \end{aligned}
> $$
>
> This is the source two-structure-constant color contraction written with the invariant matrix trace. The trace-dual completeness relation replaces its internal color contraction by a trace of commutators; its scalar-matrix subtraction vanishes because the commutators are traceless. Thus the displayed coefficient is independent of $N$.
>
> $Q_1$ has differential order two. The cross term between the two traces therefore substitutes their matrix gradients $mb^{m-1}$ and $nb^{n-1}$ into this symbol. Their scalar-matrix components drop from every commutator, so the same calculation applies to the trace-dual gradients on $\mathfrak{sl}_N$. Since the two matrix powers commute, cyclicity and the Jacobi identity identify the two symmetric probe contributions. The result is
>
> $$
> -2mn\,\operatorname{tr}\left(
> [b^{m-1},\partial_1c][b^{n-1},\partial_2b]
> -[b^{m-1},\partial_2c][b^{n-1},\partial_1b]\right).
> $$
>
> In a gauge-invariant trace, bare-ghost gauge terms cancel. The remaining tree variation replaces $\partial_\alpha b^r$ by $[\partial_\alpha c,b^r]$. Differentiate $[b^{m-1},b^{n-1}]=0$ and $[b,b^{m-1}]=0$; cyclicity then gives
>
> $$
> \begin{aligned}
> &Q_0\left[
> \left(1-\frac12J_-J_+\right)
> \operatorname{tr}\bigl((\partial_1b^{m-1})b(\partial_2b^{n-1})\bigr)
> \right]\\
> &\qquad=-\frac12\operatorname{tr}\left(
> [b^{m-1},\partial_1c][b^{n-1},\partial_2b]
> -[b^{m-1},\partial_2c][b^{n-1},\partial_1b]\right).
> \end{aligned}
> $$
>
> Multiplication by $4mn$ proves the identity. $\square$

For $(m,n)=(2,2),(2,3),(3,3)$ the primitive reduces exactly to the three trace formulas used in the SU(3) universal first lift. The source two-probe symbol was also checked as a full component polynomial for $SU(3)$ and $SU(4)$, with $564$ and $2664$ monomials. Additional full SU(3) checks at $(2,4),(3,4),(4,4),(2,5)$ have zero residual. These checks fix the source sign and coefficient; the matrix calculation supplies the general identity.
