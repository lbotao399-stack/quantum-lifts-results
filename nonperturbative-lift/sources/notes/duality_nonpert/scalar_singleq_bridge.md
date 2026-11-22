# Scalar chiral-ring relations and one-supercharge cohomology

**Established:** in a unitary four-dimensional $\mathcal N=1$ SCFT, an actual full chiral-ring relation represented by a common-$\bar Q$-closed Lorentz scalar is trivial in the cohomology of a chosen $\bar Q$ component. This statement concerns renormalized physical local operators at the same conformal fixed point. It does not identify the full chiral ring with one-supercharge cohomology, prove a UV-to-IR operator map, or establish the proposed $SU(5)$ square-zero relation.

## 1. Exact radial algebra

[Dolan–Osborn, arXiv:0801.4947](https://arxiv.org/pdf/0801.4947), (1.1), (2.11), and Appendix A, use $H$ for the dilation generator, $R$ for the exact superconformal R-charge, and $\bar J_3$ for the Cartan generator of the dotted-spin rotation group. Their superscript $+$ is the adjoint for the positive radial-quantization inner product; it is not the Minkowski adjoint. Their component $1$ is the dotted weight $-1/2$, corresponding to the chosen $\dot-$ component. No rescaling of supercharges is made:

$$
Q=\bar Q_1,
\qquad Q^+=-\bar S^1,
\qquad
\{Q,Q^+\}=2\mathcal H,
\qquad
\mathcal H=H-2\bar J_3-\frac32R.
$$

Let $|O\rangle=O(0)|0\rangle$ be the state of a renormalized gauge-invariant Lorentz scalar, with scaling dimension $\Delta$ and exact R-charge $r$. The assumptions are

$$
H|O\rangle=\Delta|O\rangle,
\qquad R|O\rangle=r|O\rangle,
\qquad \bar J_3|O\rangle=0,
\qquad
\bar Q_1|O\rangle=\bar Q_2|O\rangle=0.
$$

The dotted component labels run over $1,2$. The second diagonal anticommutator has the opposite $\bar J_3$ sign, so on this scalar both give

$$
\{\bar Q_1,-\bar S^1\}|O\rangle
=\{\bar Q_2,-\bar S^2\}|O\rangle
=(2\Delta-3r)|O\rangle.
$$

## 2. Proof

Positivity and closedness give

$$
(2\Delta-3r)\langle O|O\rangle
=\|\bar Q_1|O\rangle\|^2
+\|(\bar Q_1)^+|O\rangle\|^2
=\|\bar S^1|O\rangle\|^2\geq0.
$$

If $2\Delta-3r>0$, there is the explicit radial-state primitive

$$
|O\rangle
=\frac{1}{2\Delta-3r}\bar Q_1(\bar Q_1)^+|O\rangle
=-\frac{1}{2\Delta-3r}\bar Q_1\bar S^1|O\rangle.
$$

This case does not require a prior chiral-ring relation.

If $2\Delta-3r=0$, applying positivity to both dotted components gives

$$
(\bar Q_1)^+|O\rangle=(\bar Q_2)^+|O\rangle=0.
$$

Now assume the *full operator chiral ring* makes $O$ zero. This means there are physical local states $|X^{\dot\alpha}\rangle$, in the dual dotted-spin representation, for which

$$
|O\rangle=\bar Q_{\dot\alpha}|X^{\dot\alpha}\rangle.
$$

The repeated upper/lower dotted index is the canonical spinor contraction. Componentwise radial adjunction then gives

$$
\langle O|O\rangle
=\langle O|\bar Q_1|X^{\dot1}\rangle
+\langle O|\bar Q_2|X^{\dot2}\rangle
=0.
$$

Unitarity identifies null states with zero. Thus $|O\rangle=0$ in the physical state space and its one-supercharge class is zero. A negative value of $2\Delta-3r$ is impossible for a nonzero common-closed scalar. This proves the claim. Finite decompositions into scaling eigenoperators are handled componentwise.

The primitive $-\bar S^1|O\rangle/(2\Delta-3r)$ is an IR state–operator construction. It is not an explicit microscopic UV composite or an instanton calculation.

## 3. Why the known counterexample does not contradict the proof

[Gukov–Nawata–Saberi–Stosic–Sulkowski, arXiv:1512.07883](https://arxiv.org/pdf/1512.07883), §1.4, use a three-generator graded complex. Their $\phi$ is bosonic and $\psi_1,\psi_2$ are fermionic, with

$$
Q_1\phi=\psi_1,
\qquad Q_2\phi=\psi_2,
\qquad Q_i\psi_j=0\quad(i,j=1,2).
$$

Consequently

$$
\ker Q_1=\operatorname{span}\{\psi_1,\psi_2\},
\qquad
\operatorname{im}Q_1=\operatorname{span}\{\psi_1\},
\qquad
H(Q_1)=\operatorname{span}\{[\psi_2]\},
$$

where the brackets denote a cohomology class. Both common-closed generators are killed after quotienting by both images:

$$
\frac{\ker Q_1\cap\ker Q_2}
{(\ker Q_1\cap\ker Q_2)\cap(\operatorname{im}Q_1+\operatorname{im}Q_2)}=0.
$$

Thus a general identification with one-$Q$ cohomology is false. The toy complex supplies no unitary superconformal radial algebra, and its surviving classes are fermions; if $\phi$ is a Lorentz scalar and $Q_i$ a spinor, the $\psi_i$ form a spinor. The proof above instead requires a scalar target and uses the superconformal positive anticommutator. Physical one-$Q$ cohomology also contains spinful semishort sectors and holomorphic translation descendants.

[Ye, arXiv:1706.02723](https://arxiv.org/html/1706.02723#S2), §2, cites this example through reference [34]. The reference is *Sequencing BPS Spectra*, not a Witten 2003 paper.

## 4. Which operator space is required

The argument uses the full physical SCFT local-operator space with a positive radial inner product. It does not apply directly to a vacuum-localized effective coordinate ring, a ring with fields inverted, a reduced ring of expectation values, or a quotient obtained by setting mesons to zero. Nilpotent operator classes can vanish in every vacuum without being chiral-ring zero.

A conformal primary is annihilated by the special conformal generators. A superconformal primary is additionally annihilated by the conformal supercharges. Only the latter assumption, together with chirality, directly gives $\Delta=3r/2$ from the algebra. A common-closed scalar which is chiral-ring exact can be a supersymmetry descendant with $\Delta>3r/2$.

Accordingly, $R((\operatorname{Tr}X^2)^2)=8/3$ does not by itself establish that a nonzero renormalized square has $\Delta=4$. Bare coincident polynomials, finite renormalized composites, and the product projected into the chiral-primary ring must be matched explicitly. The state proof requires that matching as an input; it supplies no contact-subtraction prescription.

## 5. Conditional application to the cubic-adjoint benchmark

For the $SU(5)$ theory with $N_f=4$ and cubic-adjoint superpotential, the anomaly-free charge assignment used before accidental-symmetry corrections is

$$
R(X)=\frac23,
\qquad
R(Q)=R(\widetilde Q)
=1-\frac54\frac23=\frac16,
\qquad
R(M_1)=\frac13,
\qquad M_1=\widetilde QQ.
$$

If this assignment were the exact superconformal assignment for a nonzero chiral-primary meson, it would imply

$$
\Delta(M_1)=\frac32\frac13=\frac12<1,
$$

violating the scalar unitarity bound. A surviving $M_1$ therefore cannot remain an interacting chiral primary with this R-charge. In the standard free-meson resolution it has

$$
R(M_1)=\frac23,
\qquad \Delta(M_1)=1,
$$

and an accidental symmetry changes the exact superconformal R generator. This exact generator is the one required in §§1–2. A product of a unitary interacting SCFT and free mesons still satisfies the proof.

There is a useful fixed-point consistency check. The magnetic gauge group is $SU(3)$. If its cubic $Y$ coupling and $M_2\widetilde qq$ coupling survive, with the charge-conjugation-symmetric assignment, then

$$
R(Y)=\frac23,
\qquad
R(q)=R(\widetilde q)=1-\frac34\frac23=\frac12,
\qquad R(M_2)=1.
$$

For free $M_1$, the other magnetic superpotential integrand has

$$
R(M_1\widetilde qYq)
=\frac23+\frac12+\frac23+\frac12=\frac73.
$$

If its chiral-primary class survives, its dimension is $7/2$ and the corresponding superpotential coupling has dimension $3-7/2=-1/2$. It cannot be a nonzero marginal coupling of that fixed point. This is a reason to distinguish the limiting SCFT superpotential from the full UV dual superpotential. It is not a justification for dropping that coupling in vacuum or finite-scale chiral-ring relations.

Similarly, the square of the trace has exact IR R-charge $8/3$, while a flavor-contracted monomial with two $M_1$ factors and two $M_2$ factors has R-charge $10/3$. An already-established SCFT operator identity with neutral constant coefficients must separate into its distinct R-eigenspaces. But unequal IR charges alone do not establish a finite multiplicative UV-to-IR map: matching coefficients can depend on irrelevant couplings, and the leading scaling projection need not preserve coincident products. A source relation involving a strong scale and mesons therefore requires an independent RG and contact-term analysis.

The benchmark implication is consequently conditional:

$$
\begin{gathered}
\text{full IR operator chiral-ring relation }
[(\operatorname{Tr}X^2)^2]=0,
\\
\text{correct identification of its renormalized scalar representative}
\end{gathered}
\quad\Longrightarrow\quad
[(\operatorname{Tr}X^2)^2]_{H(\bar Q_{\dot-})}=0.
$$

The antecedent has not been established here. In particular, KSS §6.1's meson-zero analysis and a classical magnetic characteristic/F-term relation cannot silently be promoted to a full quantum operator identity. Quantum characteristic relations, mesonic terms, and the precise SCFT limit must be checked separately. Even after this physical scalar implication is applied, comparison with the manuscript's renormalized holomorphic/BV complex remains an additional nonperturbative identification.

## Readback

The normalization was checked against Dolan–Osborn (1.1), (2.11), and (A.1)–(A.3); the counterexample was checked against arXiv:1512.07883 §1.4. Both PDFs and extracted text are saved under `sources/`. No original research manuscript was changed.
