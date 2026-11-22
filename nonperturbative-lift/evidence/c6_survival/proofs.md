# SU(2) partner at (9,9)

Here $V_{N_b}^{N_c}(n_1,n_2)$ is the source's space of gauge-invariant words with $N_b$ copies of $b$, $N_c$ differentiated copies of $c$, and weights $(n_1,n_2)$. The following calculations restrict to rotation scalars. The perturbative page differential has

$$
d_r:E_r^{(N_b,N_c)}\longrightarrow E_r^{(N_b-r,N_c+1)}.
$$

**Proposition 1.** For the two-step complex $V_{N_b}^{k}\oplus V_{N_b-1}^{k}$, with differential $D=Q_0+Q_1$ truncated to these two $N_b$ values,

$$
\operatorname{rank}d_1^{(N_b,k)}
=\operatorname{rank}D^{(N_b,k)}
-\operatorname{rank}Q_0^{(N_b,k)}
-\operatorname{rank}Q_0^{(N_b-1,k)}.
$$

> **Proof.** Relative to the displayed direct sums,
>
> $$
> D^{(N_b,k)}=
> \begin{pmatrix}
> Q_0^{(N_b,k)}&0\\
> Q_1^{(N_b,k)}&Q_0^{(N_b-1,k)}
> \end{pmatrix}.
> $$
>
> Projection of $\operatorname{im}D$ onto its first summand has image $\operatorname{im}Q_0^{(N_b,k)}$ and kernel
>
> $$
> \{0\}\oplus\left(
> \operatorname{im}Q_0^{(N_b-1,k)}
> +Q_1\ker Q_0^{(N_b,k)}
> \right).
> $$
>
> The anticommutator relation gives
>
> $$
> Q_0Q_1x=-Q_1Q_0x=0\quad(x\in\ker Q_0),
> \qquad Q_1Q_0u=-Q_0Q_1u.
> $$
>
> Therefore
>
> $$
> \begin{aligned}
> \operatorname{rank}D^{(N_b,k)}
> &=\operatorname{rank}Q_0^{(N_b,k)}
> +\dim\left(\operatorname{im}Q_0^{(N_b-1,k)}+Q_1\ker Q_0^{(N_b,k)}\right)\\
> &=\operatorname{rank}Q_0^{(N_b,k)}
> +\operatorname{rank}Q_0^{(N_b-1,k)}
> +\dim\frac{\operatorname{im}Q_0^{(N_b-1,k)}+Q_1\ker Q_0^{(N_b,k)}}{\operatorname{im}Q_0^{(N_b-1,k)}}\\
> &=\operatorname{rank}Q_0^{(N_b,k)}
> +\operatorname{rank}Q_0^{(N_b-1,k)}
> +\operatorname{rank}d_1^{(N_b,k)}.
> \end{aligned}
> $$
>
> $\square$

**Proposition 2.** Every incoming differential $d_r$, $r\ge2$, to the rotation-scalar class at $(N_b,N_c)=(5,5)$ and $(n_1,n_2)=(9,9)$ has zero source.

> **Proof.** Its source has $(N_b,N_c)=(5+r,4)$. Each $b$ has weights $(1,1)$, while each differentiated $c$ needs at least one derivative. Thus
>
> $$
> 2\bigl(9-(5+r)\bigr)\ge4
> \quad\Longrightarrow\quad r\le2.
> $$
>
> For $r=2$, all seven $b$'s are undifferentiated and each of the four $c$'s has exactly one derivative. Bold numerals below denote the dimensions of the irreducible SU(2) representations. The adjoint is $\mathbf3$ and the rotation spinor is $\mathbf2$. Their rotation-scalar fermionic component is
>
> $$
> \left[\Lambda^4(\mathbf3\otimes\mathbf2)\right]_{\text{rotation spin }0}
> =\operatorname{Sym}^2\mathbf3
> =\mathbf5\oplus\mathbf1,
> $$
>
> while
>
> $$
> \operatorname{Sym}^7\mathbf3
> =\mathbf{15}\oplus\mathbf{11}\oplus\mathbf7\oplus\mathbf3.
> $$
>
> The two lists have no common gauge representation. The gauge-singlet multiplicity in their tensor product is zero. Hence $V_7^4(9,9)$ has zero rotation-scalar component. $\square$

**Proposition 3.** If the rotation-scalar groups

$$
E_2^{(3,6)}=E_2^{(2,6)}=E_2^{(1,6)}=E_2^{(0,6)}=0
$$

at $(9,9)$ vanish over the coefficient field used, a nonzero class in $E_2^{(5,5)}$ survives to $E_\infty$ over that field.

> **Proof.** The outgoing targets for $r=2,3,4,5$ are precisely the four displayed groups. Every later page is a subquotient of $E_2$, so these targets remain zero. For $r\ge6$, the target has $N_b<0$. Proposition 2 removes all incoming differentials. $\square$

The rank arguments below use finite-field computations as lower bounds for rational maps. Each subtraction of ranks in Proposition 1 is made after the subtracted ranks have been fixed over $\mathbb Q$.

**Proposition 4.** The rotation-scalar component of $B_{5,2}\,\partial_2C_0\,\partial_1^2C_0$ gives a nonzero rational $Q_0$ class at $(N_b,N_c)=(5,5)$ and $(9,9)$.

> **Proof.** The classical operators $B_{5,2}$ and $C_0$ are those of BGKWWY (3.18); $Q_0$ commutes with derivatives and rotations and is a graded derivation. Thus $Q_0$ annihilates the product and its rotation-scalar component. In the source's Chevalley component convention (the coefficient of the Cartan matrix in $b$ uses the trace-dual basis), exact rational evaluation gives
>
> $$
> \left.\left[B_{5,2}\,\partial_2C_0\,\partial_1^2C_0\right]_{\text{rotation spin }0}\right|_{\mathfrak t}
> =\frac12 b^5J,
> $$
>
> where $\mathfrak t$ is the Cartan subalgebra and the ledger's ordered fermion product is
>
> $$
> J=\partial_2c\,\partial_2^2c\,\partial_1c\,\partial_1\partial_2c\,\partial_1^2c.
> $$
>
> All five fermionic jets are distinct, so $b^5J\ne0$. Every term in $Q_0$ contains a Lie bracket, hence $(Q_0u)|_{\mathfrak t}=0$. The displayed operator represents a nonzero class. The independent component calculation contains 885 monomials and verifies $Q_0=0$ over $\mathbb Q$. $\square$

**Proposition 5.** The rotation-scalar source has $\dim H(Q_0)^{(5,5)}=1$ over $\mathbb Q$ and survives the incoming and outgoing $d_1$ maps.

> **Proof.** Exact character dimensions and finite-field rank lower bounds give
>
> $$
> \dim V_5^5=85,
> \qquad\operatorname{rank}Q_0^{(5,4)}\ge74,
> \qquad\operatorname{rank}Q_0^{(5,5)}\ge10.
> $$
>
> The chain relation gives $\dim H(Q_0)^{(5,5)}\le85-74-10=1$. Proposition 4 gives the reverse inequality. For the adjacent source and target,
>
> $$
> \begin{aligned}
> \dim H(Q_0)^{(6,4)}&\le14-12-2=0,\\
> \dim H(Q_0)^{(4,6)}&\le246-226-20=0.
> \end{aligned}
> $$
>
> These are precisely the source and target of $d_1$ touching $(5,5)$. $\square$

Two more target groups vanish by

$$
\begin{aligned}
\dim H(Q_0)^{(3,6)}&\le1851-1461-390=0,\\
\dim H(Q_0)^{(1,6)}&\le5231-2597-2634=0.
\end{aligned}
$$

**Lemma 6.** The primary-product witness descriptors define rational $Q_0$-closed operators.

> **Proof.** In the component engine's matrix convention, $Q_0c=c^2$ and $Q_0b=[c,b]$. The primary generators used by every descriptor are $A_2=\operatorname{tr}b^2$, $B_{1,\alpha}=\operatorname{tr}(b\partial_\alpha c)$, and $C_0=\operatorname{tr}(\partial_1c\partial_2c)$. Graded cyclicity gives
>
> $$
> \begin{aligned}
> Q_0A_2
> &=\operatorname{tr}\bigl([c,b]b+b[c,b]\bigr)
> =\operatorname{tr}(cb^2-b^2c)=0,\\
> Q_0B_{1,\alpha}
> &=\operatorname{tr}\bigl([c,b]\partial_\alpha c
> +b(\partial_\alpha c)c+bc(\partial_\alpha c)\bigr)\\
> &=\operatorname{tr}\bigl(cb\partial_\alpha c+b(\partial_\alpha c)c\bigr)=0,\\
> Q_0C_0
> &=\operatorname{tr}\bigl((\partial_1c)c(\partial_2c)
> +c(\partial_1c)(\partial_2c)
> -(\partial_1c)(\partial_2c)c
> -(\partial_1c)c(\partial_2c)\bigr)=0.
> \end{aligned}
> $$
>
> Derivatives commute with $Q_0$, and its graded Leibniz rule extends these identities to every listed product. For a weight-zero word of maximal rotation spin $9-N_b$, its scalar projector is
>
> $$
> \prod_{\substack{j=1,\ldots,9-N_b\\\text{rotation spins}}}
> \left(1-\frac{J_-J_+}{j(j+1)}\right).
> $$
>
> Here $J_+$ converts $\partial_2$ to $\partial_1$ with its multiplicity, and $J_-$ converts $\partial_1$ to $\partial_2$. On the rotation-weight-zero component of spin $j$, $J_-J_+$ has eigenvalue $j(j+1)$. The projector commutes with $Q_0$. Its denominators, the trace-dual basis denominator $2$, and the factorial denominators of the source $Q_1$ jet coefficients are invertible modulo $2147483647$. Every stored descriptor therefore specifies a rational cocycle and its valid reduction. $\square$

**Lemma 7.** Independent Cartan restrictions of $Q_1$ images of the witness operators certify the same rank lower bound for $d_1$ over $\mathbb Q$.

> **Proof.** For each witness $Q_0u=0$,
>
> $$
> Q_0Q_1u=-Q_1Q_0u=0.
> $$
>
> Cartan restriction kills $Q_0$ images. A nonzero minor of the matrix of Cartan coefficients therefore proves independence of the corresponding $[Q_1u]$ in $H(Q_0)$. The matrix is rational, with denominators invertible modulo the stated prime. A nonzero minor after reduction proves that its rational determinant is nonzero. $\square$

For the target $(N_b,N_c)=(2,6)$, the recomputed certificates give

$$
\begin{aligned}
\operatorname{rank}d_1^{(3,5)}&\ge12,
&\det\equiv1844352828\pmod{2147483647},\\
\operatorname{rank}d_1^{(2,6)}&\ge3,
&\det\equiv1717992652\pmod{2147483647},\\
\dim H(Q_0)^{(2,6)}&\le5408-3532-1861=15.
\end{aligned}
$$

**Proposition 8.** $E_2^{(2,6)}=0$ over $\mathbb Q$.

> **Proof.** The quantum identity $Q_1^2=-Q_0Q_2-Q_2Q_0$ induces $d_1^2=0$. Hence
>
> $$
> 0\le\dim E_2^{(2,6)}
> =\dim H(Q_0)^{(2,6)}
> -\operatorname{rank}d_1^{(3,5)}
> -\operatorname{rank}d_1^{(2,6)}
> \le15-12-3=0.
> $$
>
> $\square$

**Proposition 9.** $\dim H(Q_0)^{(0,6)}=6$ over $\mathbb Q$.

> **Proof.** Consider the rotation-scalar parts of products of three derivatives of $C_0$. The three derivative pairs $(r,s)$, meaning $\partial_1^r\partial_2^sC_0$, are
>
> $$
> \begin{array}{c|ccc}
> &\text{first factor}&\text{second factor}&\text{third factor}\\\hline
> 1&(5,0)&(1,2)&(0,4)\\
> 2&(2,4)&(0,1)&(4,1)\\
> 3&(2,2)&(4,2)&(0,2)\\
> 4&(1,2)&(3,3)&(2,1)\\
> 5&(1,1)&(5,0)&(0,5)\\
> 6&(1,3)&(2,2)&(3,1)
> \end{array}
> $$
>
> Lemma 6 gives rational $Q_0$-closed lifts. In the engine trace convention, restriction to the Cartan replaces $C_0$ by $2\partial_1c\,\partial_2c$. Restriction commutes with derivatives and with the rotation projector. Exact rational exterior-polynomial arithmetic gives a coefficient minor
>
> $$
> \begin{pmatrix}
> 48/7&0&0&0&-16&-16/7\\
> 16/7&16/21&-16/21&0&0&0\\
> 0&-32/105&-32/105&32/35&0&32/35\\
> 36/35&0&0&36/35&0&36/35\\
> -160/7&0&160/21&0&0&-160/7\\
> 0&0&0&0&0&-144/35
> \end{pmatrix},
> \qquad
> \det=-\frac{301989888}{420175}\ne0.
> $$
>
> Its column monomials, ordered fermions and full rational polynomials are recorded in `cartan_target_certificate.json`. The independent verifier `verify_target_certificate.py` also checks that each rational polynomial reduces to its previously stored finite-field polynomial. Because Cartan restriction annihilates $Q_0$ images, these six classes are independent. The rank lower bounds give the reverse bound:
>
> $$
> 6\le\dim H(Q_0)^{(0,6)}
> \le1241-418-817=6.
> $$
>
> The two lower bounds must therefore be attained individually:
>
> $$
> \operatorname{rank}_{\mathbb Q}Q_0^{(0,5)}=418,
> \qquad
> \operatorname{rank}_{\mathbb Q}Q_0^{(0,6)}=817.
> $$
>
> $\square$

**Proposition 10.** $E_2^{(0,6)}=0$ over $\mathbb Q$.

> **Proof.** The vanishing at $(1,6)$ gives exact rational ranks, since
>
> $$
> 5231=2597+2634,
> \qquad
> \operatorname{rank}_{\mathbb Q}Q_0^{(1,5)}=2597.
> $$
>
> The two-step rank run `RD_1_5_s83.json` used 3045 projected source samples and 3045 sketch coordinates at prime $p=2147483647$, with seed 83. It found rank 3021 and zero failed gauge/rotation invariance tests. Each source sample is the reduction of a rational gauge-and-rotation projector applied to an integer vector. Every projector denominator is invertible modulo $p$. The sketch is a linear map, so its rank is a lower bound for the rank of the unsketched differential. Thus
>
> $$
> \operatorname{rank}_{\mathbb Q}D^{(1,5)}\ge3021.
> $$
>
> Both ranks subtracted in Proposition 1 are now known over $\mathbb Q$. Therefore
>
> $$
> \operatorname{rank}_{\mathbb Q}d_1^{(1,5)}
> \ge3021-2597-418=6
> =\dim H(Q_0)^{(0,6)}.
> $$
>
> The incoming map is surjective. Its target has zero outgoing map because $d_1$ would require $N_b=-1$. Consequently
>
> $$
> \operatorname{rank}_{\mathbb Q}D^{(1,5)}=3021,
> \qquad
> \operatorname{rank}_{\mathbb Q}d_1^{(1,5)}=6,
> \qquad
> \dim E_2^{(0,6)}=6-6=0.
> $$
>
> $\square$

**Theorem 11.** For the source's nilpotent perturbative quantum differential, preserving gauge and rotation invariance and satisfying $\Delta(N_b,N_c)=(-r,+1)$ at order $r$, the rotation-scalar sector at $(n_1,n_2)=(9,9)$ has

$$
\dim_{\mathbb Q}E_\infty^{(5,5)}=1.
$$

> **Proof.** Proposition 5 gives $\dim E_2^{(5,5)}=1$. Every outgoing target at a higher page is a subquotient of one of
>
> $$
> E_2^{(3,6)}=E_2^{(2,6)}=E_2^{(1,6)}=E_2^{(0,6)}=0.
> $$
>
> The first and third vanish by the displayed $Q_0$ ranks, the second by Proposition 8, and the fourth by Proposition 10. Proposition 2 removes every incoming higher differential. The filtration is finite at the fixed weight, so the surviving class lifts to the full perturbative cohomology. $\square$

The theorem establishes existence and nontriviality of an all-order perturbative class with the stated leading component. An explicit full representative, the full quantum product syzygy $C_2\star C_4=-(5/3)C_0\star C_6$, and the nonperturbative differential on this class are separate calculations.
