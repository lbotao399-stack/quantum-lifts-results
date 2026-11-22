# Twenty first-differential multiplets in SO(19)

In the relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, two $SO(19)$ determinant sources generate independent first-differential images of dimension 64 each. Each contains ten highest-weight multiplets under dotted spin and matter flavour. The two terminal lengths are 30 and 32. The result uses supercharge order. The physical comparison is the one in Appendix B of the perturbative manuscript. [Supercharge conventions](https://arxiv.org/abs/2611.24160) · [Coefficients and highest weights](six_substitutions.json)

For $k=15,16$, take nine scalar single paths and $k-4$ double gaugino paths in the polarized Pfaffian. Each double path is the auxiliary identity times the ordered difference of the two dotted words in the existing determinant convention. Project the nine scalar markers with the cube of their three-by-three determinant. The coefficient has Pfaffian degree $k+5>19$, so it vanishes identically in nineteen colours. Its trace length $2k+1$ is odd. Applying graded reversal to an odd trace gives a primitive $O_0$ whose tree variation is this relation. In finite colours $Q_0O_0=0$. Both singlet projections are available over the rationals. Project the primitive to dotted-spin and flavour singlets if its word ordering breaks either symmetry. [Pfaffian and classical primitive](README.md) · [Graded reversal](../orthogonal_odd_trace_primitive.py)

Use the six odd substitutions that replace each of the three scalar flavours by either dotted gaugino. On higher jets their coefficients are those in the divided Taylor basis. They anticommute in pairs, square to zero and anticommute with $Q_0$. In a source or its substituted sources there is exactly one second jet. The only active pair of $Q_1$ consists of a matter second jet and a complementary first scalar. Acting on the pair with any substitution destroys its flavour content, and both gaugino outputs are annihilated. The actions on spectators cancel with the graded signs. Thus each substitution takes a first image to a first image on these grades. The checker tests all 36 ordered substitution pairs on 17 jet letters, as well as the first kernel on all 145 ordered input pairs for each substitution. [Letter action](../so8_level24_descendant.py) · [Replay](six_substitutions.py)

Apply all six substitutions to $d_1[O_0]$, in the replay order, and evaluate in the defining-trace glueball coordinate. Only products of bilinears can contribute. In every longer graded-commuting trace at least one of the two gauginos repeats. This gives the following full-colour nonzero multiples of $S^k$, where $S$ is the defining trace of the positive and negative dotted gauginos in that order. The top substitution is a spin and flavour singlet, so singlet projection of the primitive preserves its coefficient.

| $k$ | Level at the source before substitutions | Native columns tested | Ordered second-jet primitive transpose terms | Coefficient in the ordered terminal convention |
|---:|---:|---:|---:|---:|
| 15 | 84 | 2661 | 849 | 48660480 |
| 16 | 90 | 2661 | 849 | -123863040 |

The ordered terminal sign for the top product is fixed directly by the replay. The classical glueball theorem gives $[S^k]\ne0$ for $k<19-2$. Repeated double paths are evaluated with one commuting marker and its factorial multiplicity. The original labelled-path implementation gives identical coefficients **on every ordered word** in the three controls $(N,k)=(7,4),(9,5),(10,6)$. The top coefficients also agree with $(-1)^k1440(k-4)!$ times the previously derived pure first coefficient in all five cases. No fit in the colour dimension is used. [Classical glueball theorem](https://arxiv.org/abs/math/0312452) · [Pure coefficient](first_family.md)

At a fixed substitution degree, suppose a linear combination of the images vanishes. Select one subset with a nonzero coefficient and apply its complement. Every other subset at that degree repeats a substitution and vanishes. The chosen coefficient now multiplies the nonzero top image, a contradiction. Different degrees have different total dotted Taylor degree. The 64 images therefore carry the exterior powers of the dotted doublet tensored with the dual flavour triplet. Their irreducible dimensions in order of degree are

$$
1,\quad 6,\quad 9+6,\quad 4+16,\quad 9+6,\quad 6,\quad 1.
$$

The determinant source for the second $k$ has a different letter length, preserved by the first differential and all six substitutions. The two images are independent even where their physical levels coincide. This adds twenty $d_1$ highest-weight multiplets to the finite campaign count. The smaller-colour controls are excluded from this increment. The calculation supplies zero new BMN generator families. Any $d_2$ search in these sectors can quotient these first images before constructing a second-order matrix.

Run from this directory in the repository's NumPy, SciPy and python-flint environment.

~~~sh
python3 six_substitutions.py > /tmp/so19-six.json
cmp /tmp/so19-six.json six_substitutions.json
~~~
