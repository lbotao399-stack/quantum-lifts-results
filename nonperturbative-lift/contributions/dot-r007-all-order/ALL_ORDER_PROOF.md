# Characteristic-zero all-order perturbative survival of SU(3) R-007

The corrected 6139-term $j=30$, $C=-1$, dotted-singlet representative $O_0$ in the [frozen input replay](../dot-r007-exact-input) admits an all-order closed lift over $\mathbb Q[[u]]$ in the normalized relative holomorphic Taylor-jet complex. The completed class has no formal-coupling torsion. This is an existence proof, not a list of its coefficient polynomials or a physical instanton calculation.

## 1. Lift the free syzygy, not a finite-rank trace relation

Use the corrected version of [dMK et al., Eqs. (2.20), (2.24) and App. A](https://arxiv.org/html/2412.08695):
$$
O_0=\sum_a c_a(-1)^{|f_a|}f_a r_a,\qquad
Q_0r_a=R_a,\qquad \sum_a c_a f_a R_a=0.
$$
The first local identities are the six finite-$SU(3)$ identities replayed in the frozen packet. The last is stronger than a finite-colour trace relation: it is zero as a free graded-commutative polynomial in the independent graviton component symbols. The replay uses independent $u_2,u_3,v_{2,\mathrm{raw}},v_{3,\mathrm{raw}}$ with their linear symmetry/traceless projections: 22 summands with 6079 separate monomials cancel exactly in the required highest component. The algebra on the projected $v$ subspace embeds in the algebra on the raw symbols by a split linear inclusion. In the quantum rewrite we use a basis of the **projected closed** graviton components, not any potentially nonclosed raw trace.

Use one common $Q(u)$ and matching marked product throughout, from [Li–Liu–Chang, Eqs. (6.13), (6.26), (8.3), Theorem 9.1](https://arxiv.org/pdf/2611.24160v1). Replace polynomials in closed gravitons by fixed-bracketing marked lifts $f_a^\star,R_a^\star$. The difference between any original tree and its graded sorted right comb is $Q$-exact with primitive starting at order $u$: apply the double-marked homotopy of Eq. (6.26), graded commutativity (put its root equal to the strict unit), and Leibniz to carry a local rewrite through its parent tree. A repeated odd generator is treated by graded commutativity and division by two. The ordered free words then cancel as above. Hence
$$
 S=\sum_a c_a f_a^\star\star R_a^\star=QH,\qquad H=O(u).
$$
This is an equality for the chosen trees, not an assumption that a general classical algebraic relation persists quantum mechanically.

## 2. The first two local repairs

The first local residual $e_a^{(1)}=Q_1r_a-[u]R_a^\star$ is $Q_0$-closed by nilpotency and Leibniz. Project $r_a$ to $C=(-1)^{Y_r+1}$, with $Y_r$ its letter number. This preserves $Q_0r_a=R_a$ since $R_a$ is a product of $Y_r+1$ letters, and preserves the weighted sum defining $O_0$, since $C O_0=-O_0$ and each $f_a$ is of ordinary graviton parity.

The highest theta weights of $12_{03},12_{11},12_{22},14_{13},16_{01},16_{23}$ are $(3,3,0),(3,2,1),(4,2,0),(4,3,0),(3,3,2),(5,3,0)$. At order $u^k$ the undeformed cubic $Q_k$ and marked-product coefficients remove $k$ from each theta weight and insert $(k,k)$ dotted jets. Thus four highest-weight residuals, and their whole flavour tensors, vanish. The remaining two lie at $(\theta;n_+,n_-)_{Y}=(2,1,0;1,1)_4$, $C=-1$, and $(2,2,1;1,1)_5$, $C=+1$. They restrict to zero on Cartan because their parity is opposite to $(-1)^Y$.

Here $I$ is the classical incoming map, and $A$ is the stacked map $(Q_0,\rho)$ out of the indicated complete invariant highest-weight slot. No translation quotient is taken. The last row will be used in step 3.

| Slot | Exact chain dimension | Rational lower rank of $I$ | Rational lower rank of $A$ |
|---|---:|---:|---:|
| $(2,1,0;1,1)_4$ | 12 | 4 | 8 |
| $(2,2,1;1,1)_5$ | 101 | 46 | 55 |
| $(1,1,0;2,2)_4$ | 28 | 8 | 20 |

Exact integer Weyl coefficients determine the dimensions independently. Modular evaluated ranks give lower ranks of the rational maps on true invariant/HW trace polynomials. The joint ranks use the evaluated outgoing ranks $4,49,15$ and the additional Cartan ranks $4,6,5$ **on the kernels of those evaluated outgoing matrices, in the same columns**, not ranks of an unrestricted Cartan map. The two lower bounds admit no rational Cartan-null cohomology, because $\mathrm{im}\,I\subseteq\ker A$ and the sums saturate the exact dimensions. A positive modular Betti number is not being lifted as a lower bound.

Thus choose homogeneous first primitives $b_a$ with $Q_0b_a=-e_a^{(1)}$. They can be chosen equivariant by restricting to the full invariant highest-weight slot and using the flavour lowering maps. The local solutions $T_a=r_a+u b_a$ satisfy $Q T_a=R_a^\star\bmod u^2$.

## 3. Cartan injectivity forces the only remaining repair

Only $16_{01}$ can have a second-order local residual. Denote its three lower-flavour components by $e_l=[u^2](QT_{16_{01},l}-R_{16_{01},l}^\star)$. They are $Q_0$-closed, and the highest, $l=2$ in zero-based indices, is in $(1,1,0;2,2)_4$, $C=+1$. Put
$$
 W=\sum_a c_a(-1)^{|f_a|}f_a^\star\star T_a-H .
$$
Then $QW=0\bmod u^2$, and its actual second residual is $\sum_l F_l e_l$ since the earlier local defects are zero. The two exact prefactors $s22,s23$ give
$$
 \sum_l F_l e_l =\sum_{a,b,c=0}^{2}\epsilon^{ab0}
   \left(\frac85 v^c_a v^{00}_c e_b+
         \frac{16}{15}v^0_a v^{0c}_b e_c\right).
$$
The [all-page Cartan-image criterion](../dot-cartan-image/PROOF.md) applies to the actual second residual of $W$, so its Cartan restriction is zero, including any change by previously exact classes. This is not a claim that a bare $Q_2 O_0$ has zero restriction.

Let $E=\rho(e_2)$, and $P$ be the even cyclic flavour permutation $2\mapsto0\mapsto1\mapsto2$, acting with the usual wedge signs on theta. The ambient Weyl-invariant flavour/dotted HW map
$$
 E\longmapsto \rho(F_2)E+\rho(F_0)P E+\rho(F_1)P^2 E
$$
is injective. In traceless diagonal coordinates the replay stacks two simple Weyl invariance tests, two flavour raising conditions, one dotted raising condition, and this multiplication in the same 300 raw monomial columns. Exact rational elimination selects 300 original rows whose determinants are $201312\bmod1000003$ and $1551631673\bmod2147483647$. The prefactors, with supports $14,36,36$, are independently rebuilt from the frozen $s22,s23$ trace expressions. Thus $E=0$ and covariance gives $\rho(e_l)=0$ for all three components. The last small-rank row now gives a second local primitive. There are no possible local residuals of order three and above, because the highest initial $\theta_3$ of $16_{01}$ is only two.

Inserting the corrected $T_{16_{01}}$ gives $QW=0$ identically. If needed, project $W$ to the combined loop/jet/flavour grading, which the differential preserves; the global highest initial minimum theta weight is four. The normalized fan coefficients of the resulting finite differential system are rational by LLC Eq. (2.5), so existence over an extension field (from the marked primitives) implies a rational solution by finite linear algebra. We do not need a separately normalized rational coefficient for every term of the bracket homotopy.

Finally the zero-jet projection composed with the [exact SU(2) block-separating dual](../dot-r007-exact-input) kills every formal boundary and has nonzero constant pairing with the repaired representative. The class is therefore not formal-coupling torsion.

## Replay and limits

The companion scripts verify the free polynomial cancellation, independent integer characters, and the joint small rank certificates, including the original 300-row Cartan minor. None requires a large $j=30$ matrix or explicit all-order coefficient polynomials. The conclusion applies only to the declared formal twisted operator complex. It leaves both the physical comparison and the instanton mixing coefficient open.

