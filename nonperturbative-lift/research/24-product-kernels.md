<a id="np-section-22"></a>

## 24. Appendix: one-loop and two-loop product kernels

This appendix uses the orthonormal component convention of the source product-kernel subsection.

The local normalization bridges in the body relate it to the trace-dual engine convention.

The source generators are $t_A$.

Here $A,B,C=1,\ldots,N^2-1$ label adjoint basis elements, with supplied metric $\delta_{AB}$ and structure constants $[t_A,t_B]=f_{AB}{}^Ct_C$, where $f_{ABC}=\delta_{CD}f_{AB}{}^D$.

The generator normalization is

$$
\sum_{\substack{C,D=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
f_{ACD}f_{BCD}=N\delta_{AB}.
$$

The jet directions are $\alpha,\beta=1,2$.

The source antisymmetric tensor satisfies $\epsilon^{12}=1$, with $\partial^\alpha=\epsilon^{\alpha\beta}\partial_\beta$.

The following action fixes the source loop coefficient; $A,B$ are fixed color labels of the two input letters:

$$
Q_1(b_Ac_B)=\frac12
\sum_{\substack{C,D,E=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
f_{ACD}f_{BCE}\,
\partial_\alpha c_D\,\partial^\alpha c_E,
$$

$$
Q_1(b_Ab_B)=\frac12
\sum_{\substack{C,D,E=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
f_{ACD}f_{BCE}
\bigl(\partial_\alpha c_D\,\partial^\alpha b_E
-\partial_\alpha b_D\,\partial^\alpha c_E\bigr).
$$

These equations fix the source $\hbar$ in this appendix.

The normalization entry for the electronic certificates gives the physical-loop factor of the coefficients stored in the body.

<a id="np-74b4ed7f5ce47c5f-1"></a>

### One-loop triangle

The source takes $A=R X_A(z_1)$ and $B=Y_B(z_2)R'$ as the two operators.

In the complete expressions, $A,B$ denote operators; as subscripts of $X_A,Y_B$, they denote the fixed adjoint labels of the selected letters.

The ordered remainders $R,R'$ result from removing the selected letters; $z_1,z_2$ are their displacements.

Apply $\partial_{z_1},\partial_{z_2}$ and then set the displacements to zero to obtain derivative letters.

Each selected pair contributes

$$
\begin{aligned}
\sum_{\substack{C=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
f_{ABC}
\int_0^1dr_1\int_0^{1-r_1}dr_2\,
\partial_\alpha R\bigl(r_2(z_2-z_1)\bigr)
\partial^\alpha Z_C\bigl(r_1z_1+(1-r_1)z_2\bigr)R'.
\end{aligned}
$$

Here $r_1,r_2$ are simplex coordinates.

For $(X,Y)=(b,c),(c,b)$, take $Z=c$; for $(X,Y)=(b,b)$, take $Z=b$.

The $(c,c)$ channel vanishes.

Sum over all selected pairs with the Koszul sign for moving $X$ to the right end of $A$ and $Y$ to the left end of $B$.

This gives $[\hbar](A\star B)$.

> **Calculation.**
>
> For selected letters with zero derivative order, the simplex area is
>
> $$
> \int_0^1dr_1\int_0^{1-r_1}dr_2=\frac12.
> $$
>
> so the pair contributes
>
> $$
> \frac12
> \sum_{\substack{C=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
> f_{ABC}(\partial_\alpha R)(\partial^\alpha Z_C)R'.
> $$
>
> The derivative acts on the complete remainder $R$.
>
> In $C_0\star{}$, one additional derivative on the selected ghost uses
>
> $$
> \int_0^1dr_1\int_0^{1-r_1}dr_2\,r_1
> =\int_0^1dr_1\int_0^{1-r_1}dr_2\,r_2=\frac16.
> $$
>
> Adding the two ghost insertions gives the source one-loop multiplication operator $\theta$:
>
> $$
> \theta(b_A)=-\frac13
> \sum_{\substack{B,C=1,\ldots,N^2-1\\\delta\text{-orthonormal adjoint contraction}}}
> f_{ABC}(\partial_\alpha\partial_\beta c_B)
> (\partial^\alpha\partial^\beta c_C),
> \qquad\theta(\partial^m c_A)=0.
> $$
>
> In the last equation, $m$ denotes a holomorphic derivative multi-index.
>
> Input-jet derivatives act on every displacement dependence in the simplex formula.
>
> In general, $\theta(\partial_\alpha b_A)\ne\partial_\alpha\theta(b_A)$.
>
> $\square$

The complete one-loop Leibniz equation is

$$
\begin{aligned}
&Q_0[\hbar](A\star B)-[\hbar]((Q_0A)\star B)
-(-1)^{|A|}[\hbar](A\star(Q_0B))\\
&\qquad=-Q_1(AB)+(Q_1A)B+(-1)^{|A|}A(Q_1B).
\end{aligned}
$$

Here $|A|$ is Grassmann parity.

<a id="np-1f6cef9d03fcf705-1"></a>

### Two-loop forests

In a product forest, internal vertices are interaction vertices, leaves attach to input letters, and each root carries one output letter.

At two loops, the three forest topologies are two independent triangles and two ordered single-tree channels.

Let $p_\alpha$ act on the remainder of $A$ and $q_{1\alpha},q_{2\alpha}$ on the two output letters.

The source antisymmetric pairing is

$$
\langle p,q\rangle=\epsilon^{\alpha\beta}p_\alpha q_\beta.
$$

For zero input-jet derivatives, the two-triangle kernel is

$$
\begin{aligned}
\frac1{24}\langle q_1,q_2+2p\rangle\langle q_2,q_1+2p\rangle
={}&\frac16\langle p,q_1\rangle\langle p,q_2\rangle\\
&+\frac1{12}\langle q_1,q_2\rangle
\bigl(\langle p,q_1\rangle-\langle p,q_2\rangle\bigr)
-\frac1{24}\langle q_1,q_2\rangle^2.
\end{aligned}
$$

The two single-tree channels consume two letters from $A$ and one from $B$, or one from $A$ and two from $B$, respectively.

Their kernels are

$$
-\frac1{24}\langle q,p\rangle^2,
\qquad
+\frac1{24}\langle q,p\rangle^2,
$$

Here $q_\alpha$ acts on the tree's single output letter.

Taylor coefficients of the original displacement kernels give arbitrary input jets.

The displayed kernels have zero input-jet order.

Every leaf assignment retains odd-letter order.

The forest weight is the reciprocal of its automorphism-group order, and the overall normalization is the square of the source one-loop normalization.

The source forest-topology counts at one through six loops are

$$
1,\quad3,\quad8,\quad23,\quad65,\quad194.
$$

<a id="np-1f7d097fe2a68a47-1"></a>

### The source hypothesis for the ring laws

We use [Theorem 5.4 of the all-loop manuscript](../sources/all-loop/theorems.md#b-00513) as an external input: if the interaction satisfies the quantum master equation, its induced product on $H(Q_{\mathrm{pert}})$ is associative, graded commutative, and unital.

The source [Theorem 5.2](../sources/all-loop/theorems.md#b-00465) reduces this condition to chiral gauge-anomaly cancellation.

Pure SYM has empty matter representation and satisfies the condition.

These statements hold in the perturbative Ward/holomorphic model.

For quantum cocycles $A,B,D$, the theorem gives

$$
[(A\star B)\star D]=[A\star(B\star D)],
\qquad [A\star B]=(-1)^{|A||B|}[B\star A].
$$

The nilpotency calculation in the body combines this ring theorem with two independent exact boundary certificates.

Its hypotheses are independent of the strict endomorphism identities $[L,h]=0$, $h^2=0$ and of the nonperturbative exact-product assumption.

Finite common-root calculations verify low-order instances.

The source all-order proof supplies the external framework used here.
