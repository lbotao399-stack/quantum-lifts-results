# One-loop common-root coherence audit

Source: the common-root operation and ordered bracket in `allloop_repo/main.tex`, equations `eq-ring-common-root-homotopy`, `eq-A-28`, and `eq-D-lambda`; the companion's displayed combination

$$
H_{A,B}(X)=\mathcal R(AB,X;1)+\mathcal R(X,A;B)+A\star\mathcal R(X,B;1).
$$

Original sources were not modified. All arithmetic below is exact over the rational numbers in the preserved component engine's conventions. The superscript `(1)` means the coefficient of its loop parameter.

## The attempted path-mark replacement fails

Putting the auxiliary marks along the path from A to B to X keeps X as the output root. Its naive collapse formula, however, omits derivatives of the surviving auxiliary field at the integrated B vertex. Consequently its ternary operation does not alone give the desired associator homotopy.

In SU(3), take A=C0 squared, B=C0, and X=tr(b squared). The full tree differential of X is zero. The incomplete candidate, even with its coefficient fitted on A=B=C0, leaves13894 nonzero rational monomials in the associator identity. The b_e12 test leaves1787. Both use the full differential, including undifferentiated ghosts.

The reachability and Laman-articulation argument still proves that a genuine root-X graph with pure-c other inputs contracts only b-jets at X. It does not establish that the incomplete path operation is a chain homotopy.

## Actual common-root operator

`common_triangle.py` evaluates the free common-root triangles and retains derivatives of all surviving letters. The ternary ordered-bracket sign cancels the explicit sign in the definition of R. A factor2 converts its raw current-kernel output to the stored engine convention, as checked by the complete associator equation. This is a kernel-to-engine conversion, not a change of the theory's normalization.

Checks include SU(2) and SU(3), A=B=C0, A=C0 squared and B=C0, b_e12 and b_h1, first b-jets, and the gauge-invariant tr(b squared) input.

For every single b-jet, the individual free common-root triangles vanish: its sole physical contraction consumes that insertion completely, so its source momentum is zero and the source wedge vanishes. Thus

$$
H_{A,B}^{(1)}(b_A)=0.
$$

Nevertheless, for SU(2), H_(C0,C0) applied to tr(b squared) has177 raw monomials, and its commutator with C0 on that gauge-invariant input has461. Hence the common-root homotopy is not a derivation in the b-jets with pure-c coefficients, and strict commutation is false already at this order.

## The C2 coherence obstruction is exact

Let a2 be the source's leading C2 representative. The raw commutator has306 monomials and vanishing full Q0 differential. There is an explicit2092-term primitive:

$$
[C_0,H_{C_0,C_0}^{(1)}]a_2=Q_0u.
$$

The raw-to-engine factor is retained in `common_commutator_C2_boundary_exact.pkl`; both sides acquire that factor. The full component differential gives exactly zero residual. This establishes trivial action on this one class at the first corrected order. It does not establish an endomorphism homotopy on the entire complex or any higher-order coherence.

The relevant gauge-and-rotation singlet sector has a stronger finite statement. At (Nb,Nc,n1,n2)=(1,7,6,6), exact character counts give dimension13. Twelve independent incoming Q0 images and a nonzero outgoing Q0 image give

$$
\dim H(Q_0)=13-12-1=0.
$$

The next sector, (1,8,6,6), has dimension1. All thirteen witnesses are exactly rotation invariant, and their full differential agrees with the implementation that discards bare ghosts. The evidence is `C2_coherence_target_acyclic.pkl`.

## The computed C2 partner is independent of these two choices at this order

The algebraically solved eta and the coefficient from the common-root graph satisfy, on a2,

$$
\eta_{\mathrm{algebraic}}(a_2)
-\left[\mu_1(x_2,a_2)-H_{C_0,C_0}^{(1)}(a_2)\right]
=Q_0u.
$$

The difference has704 monomials and the explicit primitive has1298. It is verified with the full differential in `eta_choice_boundary_exact.pkl`. Thus replacing the algebraic eta by the common-root choice preserves the first nonzero C2 partner class and its previously computed coefficient, with the original normalization labels retained.

On the source's leading C4 representative, the corresponding difference has2370 monomials and is also exactly Q0-exact. An explicit10598-term primitive is stored in `eta_choice_C4_boundary_exact.pkl` and checked with the full differential. This verifies the eta part of the change of homotopy choice for the leading C4 partner calculation; it does not by itself compare every commutativity homotopy entering that construction.

The scripts and pickle certificates in this directory are scratch evidence; no physical instanton measure or nonperturbative product has been inferred from these checks.
