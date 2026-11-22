# A discriminant test for the SU(3) instanton candidate

Let $F$ be the published 6,139-term BMN operator R-007 at level 30, $Q_0$ the classical supercharge and $A=\phi^1$ its first scalar. Write $\Delta$ for the discriminant of the characteristic polynomial of $A$ in the defining representation. There is a gauge-invariant polynomial cochain $H$ in the zero-derivative complex with

$$
Q_0H=\Delta^6F,\qquad
\Delta=\tfrac12(\operatorname{Tr}A^2)^3-3(\operatorname{Tr}A^3)^2 .
$$

Six is a sufficient exponent. The general localization theorem gives existence of a power for a Cartan-zero cocycle. The construction below fixes a power for this source and gives a finite Ward test. The full physical instanton matrix element remains open. [Classical localization, §6.2, Eq. (111)](https://arxiv.org/html/2610.22332v3#S6.SS2) · [Pinned classical closure and block witness](https://github.com/lbotao399-stack/quantum-lifts-results/tree/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-exact-input)

## A polynomial primitive

On the regular diagonal chart let $a_i$ be the three eigenvalues of $A$. Write $B=\phi^2$ and $C=\phi^3$ in this construction. The public fermion convention gives

$$
Q_0\psi_1=[B,C],\quad Q_0\psi_2=-[A,C],\quad
Q_0\psi_3=[A,B],\quad
Q_0 f=-[A,\psi_1]-[B,\psi_2]-[C,\psi_3].
$$

All three scalars are closed. Count every letter except $A$ with degree one, and split $Q_0$ into its degree-preserving and degree-raising parts. At fixed degree, the former sends the off-diagonal components of $(\psi_3,\psi_2,f)$ to the corresponding components of $(B,-C,-\psi_1)$ multiplied by $a_i-a_j$. On a monomial with $k>0$ off-diagonal entries, reverse these three maps by the graded derivation rule, dividing each replacement by its eigenvalue gap and the result by $k$. Set the result to zero for $k=0$. The anticommutator with the degree-preserving differential is the identity on the $k>0$ part. Each replaced slot contributes once and cross-slot replacements cancel in pairs.

Every source word has total second and third flavour degree eight. Every letter except $A$ carries at least one of these degrees, including $f$. The number of non-$A$ letters therefore stays at most eight throughout. The source is zero when every letter is diagonal. In fact, with just its scalars diagonal it has 3,986 monomials, each with exactly two off-diagonal fermions. At each degree take the lowest remaining term of $F-Q_0(\text{partial primitive})$. Its degree-preserving differential is zero by nilpotence. Its purely diagonal part vanishes because every tree boundary vanishes on commuting colours. Apply the preceding inverse to it and add the result to the primitive. The next residual has strictly greater degree. [Source expansion and independent matrix evaluation](cartan_contacts.py)

| Non-$A$ letters | Source trace words | Maximum inverse gap factors from the recursion | Maximum pole order per gap after restoring colour invariance |
|---:|---:|---:|---:|
| 5 | 49 | 1 | 6 |
| 6 | 866 | 2 | 8 |
| 7 | 3,096 | 3 | 10 |
| 8 | 2,128 | 4 | 12 |

To obtain the last column, restore each closed chain of colour indices as a trace with spectral projectors. For eigenvalue $a_i$ its projector is the product of $(A-a_j\mathbf1)/(a_i-a_j)$ over the two other eigenvalues. Each factor $a_i-a_j$ can occur at most once in the denominator of one projector. An index cycle with $n$ non-$A$ entries requires $n$ projectors. Insertions of $A$ already act by eigenvalues, and a trace containing only $A$ is a polynomial. Commutators split a matrix path into paths with one additional matrix entry, preserving closed contractions. Summing the paths is invariant under simultaneous relabelling of the eigenvalues. Multiplication by $\Delta^6$ clears all twelve possible poles per gap. The resulting expression is symmetric and polynomial in the eigenvalues with trace-polynomial coefficients. Replacing the elementary symmetric functions by the characteristic coefficients of $A$ gives a gauge-invariant polynomial on all matrices. Its tree derivative agrees with $\Delta^6F$ on the dense regular semisimple locus, hence coefficient by coefficient everywhere, including the exterior coefficients.

The replay reconstructs the full restricted polynomials on two embedded two-colour blocks. In each case it retains both roots, all diagonal and off-diagonal entries of the remaining traceless block fields, and the third diagonal entry of $A$. The two choices have distinct gaps. Each gives 97 nonzero rational source coefficients, 80 primitive coefficients, and zero coefficients in the residual. The unrestricted three-colour conclusion follows from the construction above. [Verifier](discriminant_ward.py) · [Exact output](discriminant_ward.json)

## A Ward identity at every finite input order

Let $D$ send polynomial cochains to polynomials in Cartan jets by applying a proposed leading contact and restricting the output to the Cartan. Let $\delta$ be the Cartan restriction of $\Delta$. Suppose $D Q_0=0$ on all gauge-invariant inputs and the iterated commutator with multiplication by $\Delta$ vanishes after $m+1$ steps, using multiplication by $\delta$ on the output. A polynomial contact that differentiates at most $m$ scalar inputs has this property, by the ordinary product rule. The bound must apply to the same operator on all scalar dressings of the source.

Inverting $\delta$, the sequence $\delta^{-t}D(\Delta^t F)$ is polynomial in the integer $t$ of degree at most $m$. Moving $D$ past $t$ copies of the multiplication operator gives the binomial coefficient $\binom{t}{s}$ on the $s$th commutator. For every $t\ge6$, $\Delta^t F=Q_0(\Delta^{t-6}H)$, so the polynomial vanishes. Evaluation at zero gives $D(F)=0$. The scalar $\delta$ is a non-zero-divisor on each exterior coefficient. For a contact with at most one scalar input, the two required tests take the form

$$
\delta^7 D(F)=7\delta D(\Delta^6F)-6D(\Delta^7F)=0.
$$

For Dot's 5,080-term outgoing class, keep $A$ and the same two scalar Taylor derivatives as in the three-jet test. Give them the diagonal entries $(1,2,-3),(1,0,-1),(0,1,-1)$ respectively. The discriminant is 400. Direct evaluation of all five surviving words gives $-14112/5$. Thus a globally defined finite-order contact satisfying the stated Ward identity cannot send the bare source to this Cartan polynomial with a nonzero coefficient. Establishing these hypotheses for the first coefficient of the matched instanton supercharge, including contributions from the perturbatively corrected source and physical operator mixing, is a separate calculation. A Cartan-invisible output also remains possible. [Five trace coordinates](README.md) · [Native target reconstruction](discriminant_ward.py)

At the original charge alone, resolve the mixed four-input colour tensor into its three independent flavour contractions. For the input flavours $i,j,k$, contracted scalar flavour $h$, and output $s$, the three weights are $\delta_i{}^h\epsilon_{jks}$, $\delta_j{}^h\epsilon_{kis}$ and $\delta_k{}^h\epsilon_{ijs}$ in the program's flavour basis. Together with the previous 86 pure-fermion columns these give 89 columns. All 120 specified source coordinates and the first four boundary projections have rank eight with or without the target. Two further boundaries, encoded by the trace layouts `0f0ad0bce` and `a0fb/c0d0e`, raise the coefficient rank to nine and the augmented rank to ten. The exact augmented minor and a primitive integer left relation give

$$
\det=\frac{83882060447599951121259}{410516225440},\qquad
\text{left relation on the required values}=\frac{7364131507059529316544}{5}.
$$

These six boundary tests already exclude this particular fit at fixed charge. The discriminant argument covers other colours and arbitrarily high finite input orders when the multiplication and Ward hypotheses hold. [Boundary conventions](local_contacts.md) · [Full rational test](mixed_flavours.py) · [Minor and relation](mixed_flavours.json)

From this directory run

~~~sh
python3 discriminant_ward.py | diff -u discriminant_ward.json -
python3 mixed_flavours.py | diff -u mixed_flavours.json -
~~~
