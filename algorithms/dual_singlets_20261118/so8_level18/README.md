# The two parity components of the $SO(8)$ level-18 BMN generator

At doubled charges $(2J_1,2J_2,2q_1,2q_2,2q_3)=(n,n,6,6,6)$ the perturbative supercharge maps degree $n+1$ to $n$. The derivative grade $r$ is half the total dotted Taylor degree. There are $9-n-r$ letters. The table uses reflection-even spin and flavour singlets. Triality permutes the defining and two chiral spinor representations of the colour Lie algebra.

| $n$ | Cochains | Incoming classical ranks in increasing $r$ | Full incoming rank | Perturbative cohomology | Triality-invariant part |
|---:|---:|---|---:|---:|---:|
| 0 | 81 | 17, 38, 18, 4 | 77 | 4 | 3 |
| 1 | 173 | 16, 51, 26, 0 | 93 | 3 | 2 |
| 2 | 188 | 13, 51, 26, 0 | 90 | 5 | 4 |

The three full ranks equal the sums of their diagonal ranks. Elimination of the diagonal blocks therefore leaves rank zero at every positive filtration shift. In particular the leading term of a boundary in any grade belongs to the classical image there. Evaluation of the triality average on independent integer frames gives an idempotent commuting with the differential in all 78,309 entries. One reflection-even line in each cohomology is in the two-dimensional standard representation of the permutation group. [Rational complex](perturbative_complex.json.gz) · [Triality and ranks](triality.json)

## The odd representative

Write $O$ for the recorded 579-term BMN operator at $n=1,r=0$. The reflection-even classical quotient has dimension $34-17-16=1$, and the triality-invariant quotient in that grade has dimension zero. Let $\rho$ cyclically permute the three colour representations. Evaluate every trace of $O$ in the first chiral representation, then subtract the evaluation in the second. This gives $(\rho-\rho^{-1})O$. At the saved integer probe its value is $7{,}494{,}112-10{,}663{,}816=-3{,}169{,}704$. Reflection gives $3{,}169{,}704$. In classical cohomology,

$$
(1+\rho+\rho^{-1})[O]=0,\qquad
(\rho-\rho^{-1})^2[O]=-3[O]\ne 0.
$$

Thus the odd class and the even class are independent over the rationals. The earlier finite-field BMN census already detected the odd singlet. The spinor difference specifies its representative. [Trace evaluations and classical ranks](triality.json)

The universal trace variation $Q_0O$ has 207 terms and a coefficient $-60$. In the polarized defining indices, orthogonal metric pairings of $2k$ slots are independent for $k\le8$. Assign a different defining index to each pair of a chosen pairing. Its value is one and every other pairing vanishes. Applying the adjoint and identical-letter symmetry projectors preserves injectivity on their images. At eight letters an equality of even trace polynomials consequently lifts to the universal trace space, including any $Q_0$ boundary. Every proper even factor has at most seven letters, so its classical closure also lifts universally. A product of such factors cannot represent $O$. A BMN descendant with one matter fermion obtained from an even scalar by the flavour supersymmetry also has a universally closed representative.

Two odd factors must each use four letters to supply the orientation tensor. The scalar factor has flavour Dynkin labels $(4,0)$. The factor with a matter fermion lies in $(3,0)\otimes(0,1)=(3,1)\oplus(2,0)$. Neither summand is $(0,4)$, the dual of $(4,0)$, so the product contains zero flavour singlets. Triality preserves products, letter number and the BMN symmetry action. Applying $\rho-\rho^{-1}$ to a product or descendant expression for the odd partner would give such an expression for $-3[O]$. This adds one parity-resolved generator to the campaign. The two parity components belong to one colour-triality multiplet.

## Nonperturbative channels

The argument assumes the exact cohomological primary product of protected single gravitons, the superconformal Ward action, and flat-space triality equivariance on local adjoint operators. At $n=2,r=1$, five protected products restrict independently to the rank-four Cartan. Their determinant is $30{,}084{,}825{,}196{,}396{,}544$. Classical boundaries restrict to zero on these graded-commuting letters. The preceding rank test gives the same independence for the leading terms modulo all perturbative boundaries. These five classes span the reflection-even space at $n=2$. Thus their outgoing differential vanishes. Equivariance carries this result to the odd line in the same standard representation. At $n=1,r=2$, two products built from quadratic traces are triality-invariant and have determinant $-6{,}605{,}696$. They span the invariant middle space. [Products, colour assignments and determinants](protected_products.json) · [Protected product](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

The remaining possible outgoing partner is in the standard part at $n=0$. Its even line has coordinates $(-1/4,-1/6,1,0)$ in free columns $17,19,20,76$ of the full boundary quotient (numbered from zero). Acting with the two superconformal charges that replace the two gauginos by the first scalar leaves 23 scalar trace terms. On three commuting Cartan planes the defining and spinor traces before projection are $-8,-6,-6$. Subtracting the triality average gives

$$
-8-\frac{-8-6-6}{3}=-\frac43 .
$$

If this target were the leading instanton image, its twice-lowered image would be the leading term of an exact variation after the same Ward action. Its scalar expectation in this supersymmetric Coulomb vacuum would vanish. The coefficient $-4/3$ excludes that channel under the stated Ward identities. Since the standard representation occurs once in the target, equivariance excludes the pairing on both components. This is a protection argument. The instanton matrix element has not been evaluated. [Exact scalar evaluation](triality.json) · [Superconformal action, Appendix D, eq. (169)](https://arxiv.org/html/2610.22332v3#A4)

From the repository root, with NumPy, SciPy and python-flint installed, run

```sh
python3 algorithms/dual_singlets_20261118/so8_triality.py --out /tmp/so8-triality.json
cmp /tmp/so8-triality.json algorithms/dual_singlets_20261118/so8_level18/triality.json
python3 algorithms/dual_singlets_20261118/so8_instanton_sources.py --out /tmp/so8-products.json
cmp /tmp/so8-products.json algorithms/dual_singlets_20261118/so8_level18/protected_products.json
```
