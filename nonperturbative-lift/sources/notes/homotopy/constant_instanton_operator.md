# A limitation of constant ghost-contraction ansatzes

The source tree differential is $Q_0c^i=\tfrac12 f^i{}_{jk}c^jc^k$ and commutes with the two jet derivatives. The indices $i,j,k$ are Chevalley basis labels; $f^i{}_{jk}$ are the unchanged source structure constants. The coefficient field has characteristic zero.

**Proposition.** On the full jet algebra of $\mathfrak{sl}_N$, a finite constant-coefficient differential operator involving only derivatives with respect to $c$-jets, with positive differential order, cannot commute in the graded sense with $Q_0$ unless it is zero. The same conclusion holds for the reduced positive-ghost-jet algebra used by the source.

> **Proof.** Put the odd jet derivatives in a fixed order and discard products containing a repeated derivative, which vanish. Choose a nonzero term of highest differential order in this expression. Only finitely many ghost jets occur among its differentiation variables. Choose a positive integer $t$ larger than twice every first jet index in that finite set.
>
> For any jet $(r,s)$ appearing in those variables, the indicated two-factor part of the source tree transformation is exactly
>
> $$
> \left.Q_0(\partial_1^{r+t}\partial_2^s c^i)\right|_
> {(\partial_1^r\partial_2^s c)(\partial_1^tc)}
> =\binom{r+t}{r}f^i{}_{jk}
> (\partial_1^r\partial_2^s c^j)(\partial_1^tc^k).
> $$
>
> The two ghost jets are distinct. The two ordered Leibniz terms combine with the factor $1/2$ in $Q_0$ to give the displayed coefficient.
>
> In the commutator, contracting the first factor produces a coefficient $\partial_1^tc^k$ and replaces one differentiation variable at $(r,s)$ by a differentiation variable at $(r+t,s)$. No original differentiation variable has that first index. Different original jet pairs give different new jet pairs. Contributions with two contractions against the quadratic tree coefficient have one fewer differentiation variable and no remaining ghost; they cannot cancel this term.
>
> Fix the remaining differentiation variables and the new jet pair. Vanishing of the coefficient for every $i,k$ says that the vector of coefficients in the removed color slot commutes with every element of $\mathfrak{sl}_N$. Its center is zero, so every such coefficient vanishes. This contradicts the choice of a nonzero term of highest order. Repeating by descending order proves the claim. For the reduced algebra every original ghost jet is positive and both factors in the displayed equation remain positive; deleting bare ghosts therefore does not affect this argument. $\square$

This explains why the finite cubic constant-coefficient ansatz tested in `../instanton/test_cubic_ansatz.py` fails even after matching its value on $x_2$. It does not rule out infinite jet completions or cancellations with other loop coefficients in the full instanton equation. It is not a construction of $Q_{\rm inst}$.

**Proposition (finite polynomial coefficients on the full algebra).** Allow the ghost-only differential operator to have polynomial coefficients in finitely many ghost jets, while its differentiation variables still have finite jet support. If it graded-commutes with $Q_0$ on the full jet algebra, its positive-order part vanishes. Its remaining multiplication coefficient must be $Q_0$-closed. The same statement holds for the full reduced positive-jet algebra.

> **Proof.** Write each term with its polynomial coefficient on the left and its odd derivatives in a fixed order. Choose $t$ beyond both the coefficient-jet support and differentiation-jet support. Commuting $Q_0$ with a coefficient retains the original differentiation slots. It therefore cannot supply any of the shifted highest-order slots $(r+t,s)$ used above. The coefficient of such a slot still has a fresh factor $\partial_1^tc^k$, absent from every original coefficient. Although the coefficient algebra contains Grassmann zero divisors, multiplication by this fresh generator is injective on that subalgebra. Extract its coefficient and then every original ghost monomial: the zero-center argument applies to each ordinary color coefficient separately. The positive-order part vanishes by descending order. The remaining commutator is multiplication by the tree differential of the coefficient, proving the last assertion. $\square$

Finite jet support is essential to this statement. A locally finite sum acting on all jet orders, as in an ordinary field-theoretic transformation, is not covered. The polynomial-coefficient extension here is asserted on the full algebra only; it does not establish the corresponding restriction to gauge-invariant test polynomials.

**Proposition (restriction to gauge invariants).** The same conclusion holds if the graded commutator is only required to vanish on constant-gauge-invariant polynomials in the reduced positive ghost jets. Both the reduced-positive-jet hypothesis and finite jet support are required here.

> **Proof.** Use a highest-order coefficient and its shifted differentiation slot from the preceding proof. Pair each active ghost factor with a distinct fresh positive ghost jet by the source invariant trace, and multiply these two-traces. This is a constant-gauge-invariant test polynomial. Choose every spectator jet outside the finite original differentiation support, outside its selected shift, and different from $(t,0)$; there are arbitrarily many such jets.
>
> Extract the output coefficient containing every spectator ghost and exactly one ghost at $(t,0)$. This coefficient can arise only by differentiating all active factors. Differentiating a spectator instead leaves an active ghost in the output and removes that spectator. A commutator contribution with two contractions into the quadratic tree coefficient has no coefficient ghost and cannot supply $(t,0)$; a lower-order contribution leaves additional active factors. The extracted coefficient is consequently the same highest-order commutator tensor as above, with every remaining color slot paired against an independently variable spectator color through the nondegenerate trace form.
>
> Nondegeneracy forces that tensor to vanish. In particular the removed color slot commutes with every element of $\mathfrak{sl}_N$, and is zero. Descending in differential order proves the assertion. If several active factors have the same jet label, the fixed exterior ordering gives the corresponding antisymmetric tensor; the distinct spectator jets detect each of its coefficients. Positive active jets ensure that the shifted slot $(r+t,s)$ is different from the coefficient jet $(t,0)$. $\square$

The invariant-restriction statement does **not** extend as written to an algebra retaining bare ghosts: the graded commutator of $Q_0$ with $\partial/\partial c^i$ is the constant-gauge derivation, and therefore vanishes on gauge-invariant polynomials. This is an explicit counterexample to dropping the positive-jet hypothesis. The full-algebra proposition remains unchanged.
