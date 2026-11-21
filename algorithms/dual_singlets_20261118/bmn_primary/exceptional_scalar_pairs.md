# Four two-scalar BMN generators at exceptional gauge groups

The four flavour highest weights below give indecomposable rational scalar classes in classical BMN cohomology. Each has a nonzero expectation value on the reduced Higgs branch. They add four explicit representatives to the campaign. Their multiplicities and quantum numbers appear in Bonetti, Meneghelli and Rastelli, equations (B.4e)–(B.4h). Hunziker also gives a bidegree-$(4,4)$ polynomial for the Weyl group of $F_4$. The campaign counts a generator beyond its audited Casimir families. These calculations establish neither a new literature priority nor rank-changing fortuity. [Bonetti, Meneghelli and Rastelli, Appendix B.1](https://arxiv.org/pdf/1810.03612#page=65) · [Hunziker, §5](https://citeseerx.ist.psu.edu/document?doi=0d69156f2d25ae96f993ef62b21a3b810d8af22e&repid=rep1&type=pdf)

| Group | Level $j=2(E+J_L)$ | Flavour partition | Cochain dimensions in increasing nonempty letter length | Ordinary index | Free even-length Casimir multiplicity |
|---|---:|---|---|---:|---:|
| $F_4$ | 16 | $(4,4,0)$ | $1,1,4,6$ at $5,6,7,8$ | $-1+1-4+6=2$ | 1 |
| $E_6$, outer even | 16 | $(4,4,0)$ | $1,1,4,6$ at $5,6,7,8$ | $-1+1-4+6=2$ | 1 |
| $E_7$ | 20 | $(6,4,0)$ | $2,8,11,17,15$ at $6,7,8,9,10$ | $2-8+11-17+15=3$ | 2 |
| $E_8$ | 24 | $(6,6,0)$ | $3,11,18,16,20,16$ at $7,8,9,10,11,12$ | $-3+11-18+16-20+16=2$ | 1 |

At weight $(\lambda_1,\lambda_2,0)$ only the two scalars $X,Y$ and their matter-fermion letter of weight $(1,1,0)$ enter. With $k$ fermions, the diagonal weight space is the colour singlet in $\operatorname{Sym}^{\lambda_1-k}\mathfrak g\otimes\operatorname{Sym}^{\lambda_2-k}\mathfrak g\otimes\Lambda^k\mathfrak g$, where $\mathfrak g$ is the complexified gauge algebra. Subtract the result at $(\lambda_1+1,\lambda_2-1,0)$ to obtain the flavour highest weight. The length sign is $(-1)^{\lambda_1+\lambda_2-k}$. The $E_6$ outer-odd cochain dimension vanishes at every listed length. For a primitive Casimir of degree $d$, its two-flavour scalar polarizations have weights $(d-i,i)$, $0\leq i\leq d$. The one-fermion descendants have weights $(d-i,i+1)$, $0\leq i<d$. The two-raising representation $(d,1,1)$ has zero two-variable character. Taking symmetric products of the scalar weights and exterior products of the fermionic weights, retaining length $d$ for either kind of factor, gives the last column. [Exact lengths and multiplicities](exceptional_pair_indices.json) · [Two-flavour cycle-index calculation](exceptional_pair_index.py)

## An adjoint polynomial and its rank

Write $[\ ,\ ]:\mathfrak g\otimes\mathfrak g\longrightarrow\mathfrak g$ for the Lie bracket and $[\ ,\ ]^{\mathsf t}$ for its transpose using the Killing form and its tensor product. At flavour partition $(\lambda_1,\lambda_2,0)$ take the scalar polynomial

$$
-\operatorname{Tr}_{\mathfrak g\otimes\mathfrak g}\!\left(
[\ ,\ ]^{\mathsf t}[\ ,\ ]
(\operatorname{ad}_X\otimes 1)^{\lambda_1-\lambda_2}
(\operatorname{ad}_X\otimes\operatorname{ad}_Y-
 \operatorname{ad}_Y\otimes\operatorname{ad}_X)^{\lambda_2}\right).
$$

It is gauge invariant, has the displayed bidegree, and is invariant under $Y\mapsto Y+tX$ for every scalar $t$. Independence of the third flavour gives the other highest-weight condition. Since $Q_0 X=Q_0Y=0$ it is classically closed. At commuting Cartan values of $X,Y$ every classical boundary vanishes. The middle tensor operator in parentheses with exponent $\lambda_2$ acts on root spaces $\alpha,\beta$ by the determinant $\alpha(X)\beta(Y)-\alpha(Y)\beta(X)$. Pairs with a Cartan entry or opposite roots give zero. For $\alpha+\beta$ a root, take $p$ to be the largest nonnegative integer with $\beta-p\alpha$ a root. The diagonal coefficient of minus the transposed bracket times the bracket is

$$
(p+1)^2\frac{(\alpha,\alpha)(\beta,\beta)}{2(\alpha+\beta,\alpha+\beta)} .
$$

Here $(\ ,\ )$ on roots is the inverse Killing form. The Chevalley root-string coefficient has absolute value $p+1$. Pairing opposite root vectors gives $2/(\alpha,\alpha)$. Applying these identities to the three root spaces gives the displayed coefficient. For each fixed output root $\gamma$, the sum of these coefficients over ordered root inputs with $\alpha+\beta=\gamma$, plus the two Cartan inputs $2(\gamma,\gamma)$, equals one in the replay.

For each proper bidegree contributing to the target, the two-vector Weyl invariant dimension is the conjugacy-class average of the product of the two symmetric-power traces. We recover those traces from $1/\det(1-tw)$ for a Weyl group element $w$. Products of polarized ordinary traces have an evaluation minor attaining the whole invariant dimension at every such bidegree. In type $E_6$ the degree-five trace uses the 27-dimensional representation. In each other required degree the adjoint trace suffices. Thus the restrictions of **all** products of lower-degree classes are contained in the evaluated product space. Positive BMN raising descendants preserve letter length and increase Taylor degree, so cannot generate a pure-scalar polynomial at this length. Adding the displayed adjoint polynomial gives the following full ranks in the actual bidegree, before the highest-weight subtraction. [Weyl classes and generating functions](exceptional_weyl_pairs.py) · [Pinned conjugacy data](https://github.com/jay-taylor/charliepy/blob/0a70d17a7f101ad5c9e8d07bfdf907d7b06083b3/charliepy/data/rawdata.py)

| Group | Polarized products | Rank with the adjoint polynomial | Augmented determinant modulo 1000003 | Separate determinant modulo 1000033 |
|---|---:|---:|---:|---:|
| $F_4$ | 7 | 8 | 584731 | 376321 |
| $E_6$ | 7 | 8 | 302384 | 846483 |
| $E_7$ | 13 | 14 | 542951 | 312542 |
| $E_8$ | 11 | 12 | 507327 | 773082 |

The two runs use different integer Cartan points. Each also evaluates the selected determinant over the rationals. The complete point sets, polarizations, pivots and rational determinants are archived. On the exact $\mathcal N=4$ moduli space the commuting hypermultiplet scalars give the reduced Higgs branch $(\mathfrak h\oplus\mathfrak h)/W$, with $\mathfrak h$ the Cartan subalgebra and $W$ its Weyl group. A supersymmetric vacuum expectation annihilates exact operators. The four nonzero independent Higgs-branch functions therefore have protected completions. This use of the reduced branch requires no assertion about nilpotents of the commuting scheme. [First rational evaluation](exceptional_scalar_pairs_1000003.json) · [Second rational evaluation](exceptional_scalar_pairs_1000033.json) · [Higgs-branch ring and Molien series, §2 and Appendix B](https://arxiv.org/abs/1810.03612)

From this directory, using Python 3.9 or later, run

```sh
python3 exceptional_pair_index.py --out /tmp/exceptional-pair-indices.json
python3 exceptional_scalar_pairs.py --prime 1000003 --seed 70319 --out /tmp/exceptional-pairs-1.json
python3 exceptional_scalar_pairs.py --prime 1000033 --seed 96131 --out /tmp/exceptional-pairs-2.json
cmp exceptional_pair_indices.json /tmp/exceptional-pair-indices.json
cmp exceptional_scalar_pairs_1000003.json /tmp/exceptional-pairs-1.json
cmp exceptional_scalar_pairs_1000033.json /tmp/exceptional-pairs-2.json
```
