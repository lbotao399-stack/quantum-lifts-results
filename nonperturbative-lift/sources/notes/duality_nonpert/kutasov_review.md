# Audit: k=2, SU(5), Nf=4 as a scalar-ring truncation benchmark

**Verdict, scope refined 2026-11-13.** This is a controlled Kutasov–Schwimmer–Seiberg candidate for nonperturbative lifting. Removal of one generic adjoint-eigenvalue branch has a direct SQCD explanation. The stronger nilpotent relation at $s_1=0$ is the source's duality-derived **trace-sector prediction**, with its mesonic-origin restriction. The chain-level argument in `zerojet_perturbative_audit.md` proves all-perturbative survival of the specified square within the stated holomorphic Ward complex. The full renormalized operator map and the use of characteristic and Konishi identities in one common scheme remain unestablished; a microscopic nonperturbative $Q$-lift has not been proved.

**Audit history.** The earlier wording called this an “operator-ring prediction” without distinguishing the source trace sector from the full ring. The source-critical check in `new_critical_source_audit.md` makes that distinction explicit. It preserves the all-perturbative result and the ADS branch-removal calculation. The conditional quantum characteristic expansion in that new note is a diagnostic, not a derived correction to the KSS theory.

Primary source: [Kutasov–Schwimmer–Seiberg, hep-th/9510222](https://arxiv.org/pdf/hep-th/9510222), §2.3 and §6.1, especially (6.1), (6.5), (6.10), (6.12), (6.14). Keep its $s_0,s_1$ and operator normalization.

## 1. Exact arithmetic and the source map

The source superpotential and magnetic rank are

$$
W=\frac{s_0}{3}\operatorname{Tr}X^3
+\frac{s_1}{2}\operatorname{Tr}X^2,
\qquad \bar N_c=2N_f-N_c=3.
$$

For the two critical-point multiplicities $r_1+r_2=5$, source (6.4) gives

$$
c_1=\frac{s_1r_2}{s_0(r_1-r_2)},
\qquad
c_2=-\frac{s_1r_1}{s_0(r_1-r_2)}.
$$

**Proof.** Substitution into $r_1c_1^2+r_2c_2^2$ gives

| $(r_1,r_2)$, up to exchange | $\langle\operatorname{Tr}X^2\rangle$ |
|---|---|
| $(5,0)$ | $0$ |
| $(4,1)$ | $\dfrac{20}{9}(s_1/s_0)^2$ |
| $(3,2)$ | $30(s_1/s_0)^2$ |

Thus (6.12) specializes to

$$
\operatorname{Tr}X^2
\left[\operatorname{Tr}X^2-\frac{20}{9}\left(\frac{s_1}{s_0}\right)^2\right]
\left[\operatorname{Tr}X^2-30\left(\frac{s_1}{s_0}\right)^2\right]=0.
$$

Crucially, (6.10) gives the **positive** operator map

$$
\operatorname{Tr}X^2
=\frac53\operatorname{Tr}Y^2
+\frac{20}{9}\left(\frac{s_1}{s_0}\right)^2.
$$

At $s_1=0$, this is $\operatorname{Tr}Y^2=\frac35\operatorname{Tr}X^2$, in agreement with the special $j=k$ rule (3.32). The negative-sign rule for $j<k$ does not apply.

The magnetic characteristic relation (6.14), with the source magnetic couplings, becomes

$$
\operatorname{Tr}Y^2
\left[\operatorname{Tr}Y^2-\frac{50}{3}\left(\frac{s_1}{s_0}\right)^2\right]=0.
$$

Applying the operator map gives the additional electric relation

$$
\left[\operatorname{Tr}X^2-\frac{20}{9}\left(\frac{s_1}{s_0}\right)^2\right]
\left[\operatorname{Tr}X^2-30\left(\frac{s_1}{s_0}\right)^2\right]=0.
$$

Its constant term is $\frac{200}{3}(s_1/s_0)^4$; hence for $s_1\ne0$ it excludes the zero-trace branch. At $s_1=0$, the classical and proposed exact trace-sector relations are respectively

$$
(\operatorname{Tr}X^2)^3=0,
\qquad
(\operatorname{Tr}X^2)^2=0.
$$

The classical nonzero class is therefore the square of the trace, rather than a classical expectation value distinguishable at this undeformed point.

## 2. What is independently nonperturbative

For generic nonzero $s_1$, the $(5,0)$ branch has a massive adjoint and an $SU(5)$ SQCD factor with four flavors. Its exact ADS superpotential is

$$
W_{\mathrm{ADS}}=\frac{\Lambda_L^{11}}{\det M}.
$$

It has no finite supersymmetric vacuum. This explains the branch exclusion without using magnetic duality. The $(4,1)$ and $(3,2)$ branches satisfy $r_i\le N_f$ and are not removed by this ADS criterion. [Seiberg, hep-th/9402044, (3.4)](https://arxiv.org/pdf/hep-th/9402044).

Two qualifications are essential:

- These are **adjoint-eigenvalue branches**, not three versus two isolated vacua of the entire theory. SQCD moduli remain on the surviving branches. In particular, the $SU(4)$ factor of the $(4,1)$ branch has the $N_f=N_c$ quantum constraint; all its mesons and baryons cannot simultaneously vanish.
- At $s_1=0$ the adjoint is no longer massive and the above SQCD reduction is invalid. The nilpotent relation cannot be obtained simply by continuing the list of vacuum expectation values. A reduced vacuum set cannot distinguish the ideals $((\operatorname{Tr}X^2)^2)$ and $((\operatorname{Tr}X^2)^3)$.

The stronger trace-sector statement uses the source's operator map and magnetic characteristic relation, continued in its specified chiral-ring family. The relevant source passage, printed p.36, explicitly says “setting $M_j=0$ as usual.” Assuming only holomorphy of vacuum values, or taking the closure of the surviving point set, establishes neither that family's nilpotent structure nor its promotion to the full renormalized ring.

## 3. What remains conditional and what is not yet a Q-lift calculation

Accepting the KSS trace-sector presentation, this predicts truncation of the classical adjoint-trace ring and retains information beyond the reduced moduli space. The cited passage does not independently establish the same relation in the full ring with mesons and gauge-field-strength operators. A derivation cannot silently combine classical characteristic identities with uncorrected generalized Konishi equations: the renormalized definitions must agree, as explained in `new_critical_source_audit.md`.

The perturbative and nonperturbative statements now have different proof status:

1. `zerojet_perturbative_audit.md` proves that the antichiral counterpart of $(\operatorname{Tr}X^2)^2$ is closed and nonexact at all perturbative orders, with no $\hbar$ torsion, in the full stated holomorphic Ward complex. This uses the known graph differential and the nonzero classical F-term-quotient class, not superpotential nonrenormalization alone.
2. Establishing the full quantum operator relation in a common characteristic/Konishi scheme, and then constructing the odd operator whose full $Q$-image gives that relation or deriving the corresponding operator Ward identity, remains open. The ordinary two-supercharge chiral ring and the source manuscript's one-supercharge holomorphic-jet cohomology must also be related explicitly; that nonperturbative identification is not supplied by the perturbative graph proof.

The absence of an explicit $\Lambda$ in the nilpotent equation does not determine an instanton order or a coefficient in a $\Lambda$-adic differential. No such coefficient is extracted here.

## 4. The chain-level detector closes the perturbative gap; nonrenormalization alone would not

Perturbative nonrenormalization of the Wilsonian superpotential controls the superpotential. It is not an all-order theorem that every classical scalar class remains nontrivial in the renormalized local-operator complex. Generalized Konishi identities already modify composite chiral relations through anomaly terms. See [CDSW, hep-th/0211170, §3](https://arxiv.org/pdf/hep-th/0211170), and [Seiberg, hep-th/0212225, (2.3)](https://arxiv.org/pdf/hep-th/0212225) for fundamental matter. Those papers use $U(N)$; applying their formulas to this $SU(5)$ benchmark requires the traceless projection and the actual matter content.

The sufficient chain map has now been constructed in `zerojet_perturbative_audit.md`. Every nonzero $n$-loop local kernel raises each holomorphic derivative total by $n$, including on derivative inputs. Setting positive jets to zero therefore intertwines $\mathbf Q$ with $Q_0$. Setting the remaining ghosts, $\mathsf b$ and $\beta$ to zero and reducing modulo the F-term ideal gives a detector $\pi$ with

$$
\pi\mathbf Q=0,
\qquad
\pi\big((\operatorname{Tr}\gamma^2)^2\big)\ne0.
$$

All outgoing loop contractions on pure scalar inputs vanish by the exact antifield-terminal count. These two facts prove the square's all-perturbative survival. The Konishi anomaly contains positive ghost jets in this field dictionary, so it is killed by the detector.

This does **not** prove a natural ring injection of the classical F-term quotient, nor that the entire perturbative trace ring retains the classical cubic relation. A class zero modulo the F-term ideal can still carry a hidden positive-jet quantum contribution. In particular this detector says nothing about perturbative exactness of $(\operatorname{Tr}\gamma^2)^3$.

**Safe benchmark claim:** generic branch removal is explained by nonperturbative ADS dynamics; KSS predicts square-zero in its undeformed trace-sector presentation. The specified square has a proved all-perturbative precursor in the stated holomorphic Ward complex. The full renormalized trace map, compatible characteristic/Konishi identities, and nonperturbative comparison with physical chiral cohomology remain open. This is a controlled candidate, not a proved microscopic nonperturbative lift.
