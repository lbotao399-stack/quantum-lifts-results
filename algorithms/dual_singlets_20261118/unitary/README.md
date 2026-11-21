# A protected family at every unitary rank

Let $Y$ be the number of letters in the full local-operator complex, and let $\Theta_1$ count differentiations with respect to the first odd coordinate of the superfield $C$. Every letter has at most one such differentiation. In $SU(N)$, at $Y=\Theta_1$ the complete classical cohomology consists of products of Taylor coefficients of

$$
\operatorname{Tr}\!\left[\left(\frac{\partial_L C}{\partial\theta_1}\right)^p\right],
\qquad 2\le p\le N.
$$

The coefficients of these traces are free generators of the resulting supercommutative algebra. Odd generators occur at most once in a monomial. [Superfield and classical differential](https://arxiv.org/abs/2209.06728)

## Calculation

At $Y=\Theta_1$ every letter is a Taylor coefficient of the even matrix $\partial_L C/\partial\theta_1$. The classical equation $Q_0 C=C^2$ gives

$$
Q_0\frac{\partial_L C}{\partial\theta_1}
=\left[C\big|_{\theta_1=0},\frac{\partial_L C}{\partial\theta_1}\right].
$$

The coefficients of $C|_{\theta_1=0}$ are independent ghosts for changes of basis that depend on the remaining two odd and two even coordinates and are the identity at the origin. The absent constant ghost is supplied by the original requirement of colour invariance. Thus $Q_0$-closed polynomials are precisely the conjugation invariants of the Taylor coefficients of this even traceless matrix. At each fixed weight all expressions use a finite jet, obtained by truncating the even-coordinate powers above that weight. The nilpotent ideal of positive coordinate degree makes exponentiation of a change of basis with identity constant term finite.

On the open set where the constant scalar matrix admits a chosen cyclic vector, that vector and its first $N-1$ iterated images form a basis over the even coefficient ring of the finite jet. The matrix in this basis is the companion matrix of its characteristic polynomial. Conjugation by an invertible matrix has the same action as conjugation by a determinant-one matrix after removing the scalar determinant root. This root exists locally over the complex numbers and extends through the nilpotent ideal. An invariant polynomial therefore equals its value on the companion matrix. The latter is polynomial in the Taylor coefficients of the characteristic polynomial, since the entries of the companion matrix are these coefficients, zeroes and ones.

The cyclic-vector determinant is a nonzero polynomial in the ordinary scalar entries. Multiplication by it is injective even after adjoining the remaining polynomial and exterior variables. Equality on the cyclic open set gives equality of polynomials everywhere. Conversely, substituting arbitrary even traceless companion matrices gives independent characteristic-coefficient jets with their prescribed parities. Newton identities over the rationals replace the nontrivial characteristic coefficients by traces of powers $2$ through $N$. The identities and the resulting dimensions descend to the rationals. There are no classical boundaries into this length, since $Y-1<\Theta_1$.

The physical charges obey $Y-\Theta_1=q_2+q_3-J_1-J_2=R_2+R_3$, where $q_i=R_i+J_L$ and $J_1+J_2=2J_L$. An order-$r$ supercharge term changes $(Y,\Theta_1)$ by $(1-r,-r)$. It raises $Y-\Theta_1$ by one for every $r$, including $r=0$. Each trace above is a component of a single graviton. The quantum product theorem supplies a perturbatively closed completion of every product. Consequently all outgoing perturbative differentials vanish on this family, at every charge and every $N$. An incoming class would require $Y-\Theta_1=-1$ at every derivative grade. The same charge constraint excludes an incoming instanton partner among perturbative local operators. Outgoing nonperturbative protection uses the exact single-graviton multiplets and their exact cohomological primary product. [Perturbative product theorem](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L4695) · [Cohomological primary product](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

## Integer checks

In the Taylor-degree order $(n_+,n_-,\Theta_1,\Theta_2,\Theta_3)$, a trace of power $p$ contributes one generator in each weight $(i,j,p,b,c)$ with $i,j\ge0$ and $b,c\in\{0,1\}$. Its parity as a fermionic variable is $b+c$ modulo two. Its colour charge conjugation is $(-1)^p$, so every shortest-length product has colour parity $(-1)^Y$. Form the supercommutative monomials and take the dotted-spin and flavour highest weights. In the flavour Weyl subtraction the terms increasing $\Theta_1$ are empty at this length.

| Group | Level | Taylor degrees | Colour parity | Cochains | Classical rank at least | Independent Casimir products |
|---|---:|---|---:|---:|---:|---:|
| $SU(4)$ | 28 | $(0,0,8,3,3)$ | odd | 110 | 110 | 0 |
| $SU(5)$ | 28 | $(0,0,10,2,2)$ | even | 290 | 260 | 30 |
| $SU(4)$ | 24 | $(1,1,5,2,2)$ | odd | 85 | 84 | 1 |
| $SU(5)$ | 24 | $(1,1,5,2,2)$ | odd | 89 | 88 | 1 |

Each lower bound is obtained from a nonzero minor modulo 1,000,003 of an integer matrix. Together with the independent product count it fixes each rational rank. The word calculation also gives zero on 112 direct images $Q_0,Q_1,Q_2,Q_3$ of 28 individual graviton components. [Minors and operator checks](saturated_sources.json) · [Product count](../saturated_flavour.py) · [Replay](../check_saturated_flavour.py)

## Colour characters and the next sectors

For $Y$ labelled matrix slots, let $\sigma$ be a slot permutation and $\chi_\lambda$ the symmetric-group character of the partition $\lambda$ of $Y$. Schur–Weyl duality decomposes the invariant endomorphism algebra into endomorphism spaces of the Specht representations with at most $N$ rows. Permuting slots acts by conjugation on each block. Reversing traces also transposes that block in an orthogonal Specht basis. The two traces on the invariant algebra are

$$
\sum_{\substack{\lambda\vdash Y,\ \ell(\lambda)\le N\\\mathrm{irreps}}}\chi_\lambda(\sigma)^2,
\qquad
\sum_{\substack{\lambda\vdash Y,\ \ell(\lambda)\le N\\\mathrm{irreps}}}\chi_\lambda(\sigma^2).
$$

The sums run over inequivalent Specht representations, and $\ell(\lambda)$ counts rows. To obtain traceless letters, remove the scalar summand by inclusion–exclusion on whole permutation cycles. Multiply the transpose trace by $(-1)^Y$ for physical charge conjugation. Symmetrizing bosonic species and antisymmetrizing fermionic species then gives both colour parities separately at arbitrary $N$. For all 434 species patterns through eight letters the results agree with the twisted $SU(3)$ Weyl character, the two orthogonal parities of $\mathfrak{so}(6)=\mathfrak{sl}(4)$ after complexification, and the vanishing odd parity of $SU(2)$. [Character algorithm](../unitary_characters.py) · [1,302 comparisons](saturated_sources.json)

The next searches use $Y-\Theta_1=1$. The previous result eliminates every positive-order incoming pairing to these sources. Their outgoing obstructions require the full filtered map. From the repository root, run

```sh
python3 algorithms/dual_singlets_20261118/check_saturated_flavour.py --out /tmp/saturated.json
cmp /tmp/saturated.json algorithms/dual_singlets_20261118/unitary/saturated_sources.json
python3 algorithms/dual_singlets_20261118/trace_descent.py --flavours 10 2 2 --deficit 2 --groups SU-4 --screen --dense-unitary --work /tmp/su4-level28
```
