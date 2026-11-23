# Classical products in the second obstruction to the sixth glueball of F₄

In the relative Taylor complex of undeformed four-dimensional $\mathcal N=4$ SYM, put $E_1=H(Q_0)$ for the cohomology of the classical supercharge. Write $S=\operatorname{Tr}_{26}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$. At level 36 the incoming second differential $d_2$ to $S^6$ vanishes on the entire span of products of positive-degree **classical cohomology classes** in its selected Taylor diagram $(4,4,3,3)$. This includes product classes whose survival under the first differential initially has to be determined. [Full source and first-page Ward identities](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential)

The new rank calculation is the surjection onto the eleven-letter terminal Taylor diagram $(5,4,2)$. The Weyl and symmetric-group characters give cochain dimensions $(834,311,55)$ at nine, ten and eleven letters. For the same 55 rows and 55 integer colour assignments, the determinants of the classical images are

$$
\det M=87\pmod {1009},\qquad \det M=474\pmod {1013}.
$$

Thus the rational rank from ten to eleven letters is 55. The terminal classical cohomology is zero.

## The 55 predecessors

Polarize in five ordinary scalar directions. Apply column alternation to the lists of labels $(1,2,3,4,5)$, $(1,2,3,4)$ and $(1,2)$. Every such flag is a highest vector of ordinary weight $(3,3,2,2,1)$, the transpose of the Taylor diagram. Replacing a bracket of first scalar jets by their odd second jet gives its classical predecessor. All scalar matrices are in the defining 26 or the adjoint 52. Each admits an invariant symmetric form.

* Each of the first 51 polynomials is a product of traces of total length eleven. There is an odd trace. Transposition equates it to minus its reverse. Successive exchanges along the reversal express twice the trace as a sum of words with a marked adjacent commutator.
* The next two polynomials pair two open words through the 26-dimensional Casimir summand of $\operatorname{End}(26)$. Its endomorphisms are symmetric in the defining metric. Both recipes have word lengths five and six. Transposition of the length-five word inside this pairing gives minus its reversal, with the same commutator primitive.
* The last two polynomials use four defining-representation tensor slots between two invariant contractions of Albert cubics. Five consecutive inputs are alternated across four slots. Two act successively on the same slot. Pairing permutations by exchange of these two labels gives an explicit commutator in every summand. The two lower and upper tensors satisfy all eight simple-generator Ward identities in the replay. [Colour frame and Albert tensors](https://github.com/lbotao399-stack/quantum-lifts-results/tree/2a19364c8b5c32da2f713268b19a3dcc2580c6ed/reviews/independent-f4/fifth-glueball-partners)

The matrices have rational entries regular at the two primes. The Casimir projector on the 26 is unique. Contracting a basis of its image with the inverse trace form defines the same rational tensor independently of modular pivot choices. The replay rebuilds both representations, the projector, all ordered entries and both determinants.

## Products on the classical page

Write $D$ for total Taylor degree and $Y$ for the number of letters. The selected source has $(D,Y)=(14,13)$. In a homogeneous product of two positive-degree classes one factor has $D-Y=0$ and the other has $D-Y=1$. Every representative at $D=Y$ uses only first jets. Every bare $Q_r$ with $r\ge1$ would need $3r$ flavour slots on $r+1$ first jets and so annihilates this representative. [Fan supercharge and marked product](https://arxiv.org/abs/2611.24160)

For a factor with $D-Y=1$, the induced first differential changes its Taylor highest weight by $(1,1\mid-1,-1,-1)$. This follows on classical cohomology from the first-order Ward identity. The third column must therefore contain at least three boxes. A factor of degree at most eight has zero first differential. If the first-jet factor has degree at least six, the remaining factor satisfies this bound. [Filtered Taylor action](https://github.com/lbotao399-stack/quantum-lifts-results/blob/3843a42d05f064b9fb855571c03917707dee9a20/reviews/independent-g2-and-r007/excess-two-ward/README.md)

At degrees one through five the rational dimensions of the invariant tensors on labelled adjoint inputs are $0,1,1,5,16$. At degree four three metric pairings and two bracket pairings span. At degree five the ten products of a metric with the invariant trilinear bracket and six nested brackets paired with the last input span. The verifier uses both the colour characters and square evaluation minors. In degree three the sole invariant is a bracket paired with the third input. All tensors containing such a bracket become classical boundaries after polarization. Therefore the positive-degree terminal cohomology through degree five is generated by the invariant bilinear. Specialization to repeated labels with their graded signs preserves each identity.

The bilinear in two first jets has Taylor diagram $(1,1)$. By the vertical-strip Pieri rule, removing it from $(4,4,3,3)$ leaves just three complementary diagrams. The other factor has eleven letters.

| Complement | First target | Reason its first differential vanishes |
|---|---|---|
| $(4,4,2,2)$ | none | Third column of height two |
| $(4,3,3,2)$ | $(5,4,2)$ | Terminal rational dimension $55-55=0$ |
| $(3,3,3,3)$ | $(4,4,3)$ | Classical source dimension $141-80-61=0$ |

The last equality is from the earlier level-24 computation. Classical Taylor isotypic projections are exact on these covariant tensor powers and commute with $Q_0$. Decompose the two factors before multiplying and project their product to the selected diagram. Only the three listed complements contribute. This also treats a degree-four metric product by associating one bilinear with the other factor on the classical page. Every such projected class is consequently a sum of projected products of $d_1$-closed factors. A factor that is already a $d_1$ boundary contributes zero on $E_2$. [Earlier $(3,3,3,3)$ ranks](https://github.com/lbotao399-stack/quantum-lifts-results/blob/fa4192141124a4721f80a487f11205856784f525/algorithms/dual_singlets_20261118/bmn_primary/f4_scalar/README.md#application-at-level-24)

The marked ring gives the Leibniz rule for $d_2$ on $E_2=H(E_1,d_1)$. Its component along $S^6$ vanishes on every product of two positive-degree $E_2$ factors. Indeed, any unacted-on factor contributing to a pure first-gaugino target must be $S^k$ for $k\ge1$. A differential into $S^{6-k}$ would, after multiplication by $S^{k-1}$, make $S^5$ a second boundary. Its survival to perturbative infinity was established previously. The excess-one Ward identity permits projection of this argument to the selected Taylor diagram. Together these facts prove the claim on **classical** products. [Earlier argument on second-page products](https://github.com/lbotao399-stack/quantum-lifts-results/blob/03fc9a5d8a68a1b70072242a0770335df1d862d3/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential)

If a closed raw representative is divisible by $S$, the quotient is also classically closed. Multiplication by the nondegenerate exterior two-form $S$ is injective on the exterior degree of these two gaugino species occurring in the tree equation. Thus the conclusion also applies to raw cycles carrying a common $S$. The remaining source is the quotient of the 13-letter classical cohomology by the span of its positive-degree products. Its corrected quantum pairing remains open. [Exterior injectivity and earlier trace-family tests](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/scalar_colour_tests.md)

## Replay

Use Python 3.13, NumPy and SciPy. Each terminal command rebuilds the native rows independently. The first prime also checks every entry against the stored matrix.

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify_products.py | diff -u products.json -
python3 -B verify_terminal.py --prime 1009 | diff -u terminal-1009.json -
python3 -B verify_terminal.py --prime 1013 | diff -u terminal-1013.json -
~~~
