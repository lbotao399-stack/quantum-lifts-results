# Characteristic residues: existing twisted periods and the cohomology bridge

Primary source: M. Feigin and A. Silantyev, *Singular polynomials from orbit spaces*, Compositio Mathematica 148 (2012), 1867–1879, [arXiv:1110.1946](https://arxiv.org/abs/1110.1946), DOI [10.1112/S0010437X1200036X](https://doi.org/10.1112/S0010437X1200036X). Checked live against the [author repository PDF](https://eprints.gla.ac.uk/68072/1/68072.pdf), Theorem 4.8, Corollary 4.9 and the first formula in Section 5 (printed page 9; PDF page 10).

The paper classifies invariant polynomial twisted periods of Coxeter orbit spaces. For type $A_{N-1}$ it gives fractional characteristic-polynomial residues; their nonconstant degrees are the invariant degrees $2,\ldots,N$ plus nonnegative multiples of $N$. Their gradients belong to the reflection-representation family of singular polynomials for the rational Cherednik algebra. These mathematical residue formulas are existing results.

**Proposition.** The already certified SU(3) and SU(4) coefficient polynomials are members of that same residue family, with the following exact exponent and degree identifications. This is an identification of the polynomial families, not a claim that the full quantum complex is a Cherednik module.

> **Proof.** On the Cartan, the source $b$ is traceless. Its characteristic polynomial is therefore precisely the monic polynomial with zero coefficient of $z^{N-1}$ used for type $A_{N-1}$. No field or trace is rescaled.
>
> For SU(3), the coefficient of the source $a_1$ in the generating representative is
>
> $$
> [z^{-1}]\det(z\mathbf1-b)^{(k-2)/3},
> \qquad \text{polynomial degree}=k-1.
> $$
>
> For SU(4), the coefficient of the quadratic-trace combination in `../su4_cartan/tree_module_and_radial_kernel.md` is
>
> $$
> [z^{-1}]\det(z\mathbf1-b)^{(k-3)/4},
> \qquad \text{polynomial degree}=k-2.
> $$
>
> Here $z$ is the auxiliary Laurent variable, and $[z^{-1}]$ extracts its indicated coefficient. With the usual differential-residue convention, each displayed expression equals minus the corresponding residue at infinity. This sign does not change any saved source cocycle.
>
> In both cases the exponent is exactly (polynomial degree minus one) divided by $N$. The allowed nonconstant degrees become
>
> $$
> \begin{array}{c|c|c}
> &\text{nonconstant polynomial degree}&\text{excluded leading sector}\\\hline
> SU(3)&k-1\in\{2+3m,3+3m\}&k\equiv2\pmod3\\
> SU(4)&k-2\in\{2+4m,3+4m,4+4m\}&k\equiv3\pmod4
> \end{array}
> \qquad m\in\mathbb Z_{\ge0}.
> $$
>
> The integer $m$ enumerates the nonnegative degree shifts. The initial constant coefficient occurs separately at SU(3) $k=1$ and SU(4) $k=2$; the still lower classes are the already treated $C_0$ and, for SU(4), $a_1$. The exact recurrence proofs establish that these coefficients reconstruct the entire respective Cartan-kernel class. Thus the equality of polynomial families is established in these two cases. $\square$

The research contribution must be assessed at the bridge from the original finite-$N$ one-loop differential to these coefficients, the full-$Q_0$ primitive constructions, higher-loop lifts and the physical instanton/BV problem. The residue formula alone is not new. The arbitrary-$N$ cohomology-to-period correspondence has not been proved. Two focused web searches found no direct semichiral/Cherednik bridge, which is not an exhaustive novelty review.

**Proposition: the entire leading representative from the ordinary Hessian.** Take the ordinary symmetric Hessian of the scalar polynomial

$$
[z^{-1}]\det(z\mathbf1-b)^{(k+1)/N}.
$$

Contract its first matrix slot with $\partial_1c$ and its second with $\partial_2c$, keeping that ghost order. This prescription uses the source trace pairing between the dual $b$ coordinates and the Chevalley ghost coordinates. The polynomial degree is $k+2$ and the resulting tensor degree is $k$. The Hessian operation here is the ordinary scalar Hessian followed by ordered ghost contraction; there is no factor $1/2$ and no replacement by two ordered odd directional derivatives.

> **Proof.** On the Cartan, differentiation of the determinant gives the ordered contraction
>
> $$
> \begin{aligned}
> [z^{-1}]\det(z\mathbf1-b)^{(k+1)/N}\Bigg[&
> \frac{(k+1)^2}{N^2}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\bigr)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_2c\bigr)\\
> &-\frac{k+1}{N}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-2}\partial_1c\partial_2c\bigr)\Bigg].
> \end{aligned}
> $$
>
> The residue of
>
> $$
> \frac{d}{dz}\left[
> \det(z\mathbf1-b)^{(k+1)/N}
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\partial_2c\bigr)\right]
> $$
>
> is zero. Therefore the same contraction is
>
> $$
> \begin{aligned}
> \frac{(k+1)^2}{N^2}[z^{-1}]\det(z\mathbf1-b)^{(k-N+1)/N}\Bigg[&
> \det(z\mathbf1-b)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\bigr)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_2c\bigr)\\
> &-\frac{d}{dz}\det(z\mathbf1-b)
> \operatorname{tr}\bigl((z\mathbf1-b)^{-1}\partial_1c\partial_2c\bigr)\Bigg].
> \end{aligned}
> $$
>
> The bracket is polynomial in $z$: at a simple eigenvalue its two diagonal simple poles cancel, while its off-diagonal entries are products of the other characteristic factors. This proves polynomiality on the dense regular Cartan set and hence as a polynomial identity. Its explicit values are
>
> $$
> N=3:\qquad 3(a_1-zC_0),
> $$
>
> $$
> \begin{aligned}
> N=4:\qquad -4\Bigg[&
> \left(z^2-\frac14\operatorname{tr}(b^2)\right)C_0-z a_1\\
> &+\operatorname{tr}(b^2\partial_1c\partial_2c)
> -\frac14\operatorname{tr}(b\partial_1c)\operatorname{tr}(b\partial_2c)\Bigg],
> \end{aligned}
> $$
>
> where the expressions are restricted to the Cartan, so the source scalar projection is already implicit. Inserting these values gives exactly $(k+1)^2/3$ times the saved SU(3) generating representative and $-(k+1)^2/4$ times the SU(4) generating representative. All original representatives remain unscaled.
>
> A homogeneous scalar polynomial of degree $k+2\ge2$ with zero ordinary Hessian must vanish. The Hessian map is thus injective on these homogeneous period spaces. The already proved complete Cartan-kernel classifications establish the correspondence for SU(3) and SU(4), including their lowest degrees. The general determinant identity is valid for arbitrary $N$, but its identification with the full Cartan $Q_1$ kernel for arbitrary $N$ remains open. $\square$

The primary semichiral reference, [Budzik et al., arXiv:2306.01039](https://arxiv.org/abs/2306.01039), Section 5.2, treats the full single-trace tower in the planar limit; its finite-$N$ calculations concern bounded SU(2)/SU(3) sectors. Do not conflate that planar result with the all-degree finite-$N$ scalar results recorded here, or with physical nonperturbative cohomology.

Independent symbolic verification: `verify_hessian_bridge.py/json` reconstructs the scalar residues in independent traceless eigenvalue coordinates and compares every ordinary Hessian entry with the ordered first-ghost/second-ghost coefficient of the saved generating formulas. All degrees $k=0,\ldots,16$ pass for SU(3) and SU(4), including each forbidden congruence degree. These finite exact checks support the all-degree determinant calculation above. The primary source Theorem4.8, Corollary4.9 and Section5 was read again directly from the author-repository PDF during this check.

**Degree consequence of the period construction.** For every $N\ge2$ and $k\ge0$, the displayed degree-$k+2$ scalar residue, and hence its ordinary Hessian, vanishes exactly when $k\equiv N-1\pmod N$.

> **Proof.** If $k=N-1+Nm$ with $m\ge0$, the exponent $(k+1)/N=m+1$ is an integer. The integrand is a polynomial in $z$, so its coefficient of $z^{-1}$ is zero.
>
> In every other degree the exponent is a noninteger. Expand the monic characteristic polynomial in its independent coefficients of weights $2,\ldots,N$. The coefficient sought has weight $k+2$. For $N\ge3$, every integer at least two is a nonnegative combination of $2$ and $3$, so a monomial of that weight exists. Its generalized multinomial coefficient is a falling product of noninteger factors divided by positive integer factorials; it is nonzero. Distinct monomials in these algebraically independent characteristic coefficients cannot cancel. For $N=2$, the remaining $k$ are even, and the power $(k+2)/2$ of the quadratic coefficient supplies the nonzero monomial. Finally, a homogeneous polynomial of degree $k+2\ge2$ with zero Hessian is zero. $\square$

This supplies the all-$N$ period-side degree pattern. The original finite-$N$ quantum-complex identification is proved for SU(3) and SU(4) at the Cartan-kernel stage described above; extending that identification to all $N$ remains a separate equation-level calculation.
