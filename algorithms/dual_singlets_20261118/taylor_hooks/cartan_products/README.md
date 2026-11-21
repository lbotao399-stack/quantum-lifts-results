# Cartan coordinates for three second-order targets

Let $Y$ be the number of letters. In the Taylor diagram $(5,5,1^c)$, the notation $1^c$ denotes $c$ rows of length one. At each of the following physical levels, the classical target at the indicated $Y$ is spanned over the rationals by products of single gravitons. Its Cartan restriction is injective.

| Gauge group and level | $c$ | $Y$ | Five spin-highest product multiplicities | Rational target dimension |
|---|---:|---:|---|---:|
| $SU(3)$, 30 | 3 | 9 | $110,198,163,72,15$ | $110-198+163-72+15=18$ |
| $SU(3)$, 32 | 4 | 10 | $155,270,221,97,19$ | $155-270+221-97+19=28$ |
| $USp(4)$, 36 | 6 | 12 | $317,526,410,175,35$ | $317-526+410-175+35=61$ |

The products use the defining quadratic and cubic traces for $SU(3)$ and the quadratic and quartic traces for $USp(4)$. Take all components of these traces of powers of $dC$, with all outer Taylor derivatives of the specified weight. This span and its restriction to Cartan-valued fields commute with linear changes of the Taylor directions. On the Cartan every classical letter variation, a sum of graded commutators, vanishes. Independent restricted products therefore give independent classical cohomology classes. [Trace factors and Cartan evaluation](../../rank_two_taylor_products.py)

All rank differences are taken in the single field of $2{,}147{,}483{,}647$ elements. In each of the five Taylor weights, subtract the rank in the next higher spin weight from the rank in the weight itself. This gives the multiplicity of the corresponding spin representation in the product image. For example, the first $SU(3)$ level-32 entry is $789-634=155$. The five alternating entries then give the indicated Taylor representation by the hook character formula. The complete product enumerations determine the modular ranks, including the subtracted ranks. [Twenty-seven product enumerations and minors](jobs.json) · [All weight ranks](resolved.json) · [Hook character calculation](../README.md)

The degrees in Taylor slots are $13,14,16$. Their factorials are units at the calculation prime. On the multilinearized slots the Young idempotent thus commutes with reduction at this prime. A nonzero minor of the **projected** integral Cartan map after reduction is a rational lower bound on that same Schur component. The projected classical differentials give rational upper bounds $3683-2174-1491=18$, $6828-4185-2615=28$ and $46860-28039-18760=61$. These bounds coincide with the product multiplicities. The unitary rows have odd and even colour conjugation respectively. [Classical ranks and comparison](../../rank_two_taylor_product_bounds.py) · [Unitary colour sectors](../colour_parity/2147483647/) · [Symplectic weight sectors](symplectic_classical/)

The perturbative product theorem gives a closed completion for each of these target products. An incoming second differential is still possible. For a surviving source with first correction $Q_0\mathcal O_1=-Q_1\mathcal O_0$, Cartan coordinates now detect the classical class of $Q_2\mathcal O_0+Q_1\mathcal O_1$ in the indicated Taylor component. Its second-page value additionally takes the quotient by the image of $d_1$. Assignment of a physical source to this component uses filtered Taylor equivariance. [Quantum product theorem, parts (a) and (c)](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L384)

From the repository root, rebuild all product images and their saved minors with Python, NumPy, SciPy and python-flint.

```sh
d=algorithms/dual_singlets_20261118
python3 "$d/rank_two_taylor_product_bounds.py" --rebuild-ranks --replay-minors --out /tmp/n4-cartan-products.json
cmp /tmp/n4-cartan-products.json "$d/taylor_hooks/cartan_products/resolved.json"
```
