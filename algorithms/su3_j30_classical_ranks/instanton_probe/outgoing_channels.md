# One classical ring generator in a full-jet SU(3) channel at level 30

In the relative holomorphic jet complex of $SU(3)$, take flavour degrees $(3,0,0)$, dotted degrees $(4,4)$, dotted spin zero, charge conjugation odd and length seven. The classical cohomology over the rationals has dimension **34**. Products of positive-degree classical classes span **33** dimensions. The remaining one-dimensional quotient proves the existence of one fortuitous ring generator in this grading. Its superconformal parent and an explicit representative of the quotient remain open.

## Classical cohomology

Here $Q_0$ is the classical supercharge and increases length by one. The table gives cochain dimensions in the indicated charge-conjugation sector and certified lower bounds on the ranks of the adjacent maps.

| Dotted degrees | Flavour degrees | Length | Cochains | Incoming rank | Outgoing rank | Cohomology dimension |
|---|---|---:|---:|---:|---:|---:|
| $(1,1)$ | $(6,3,3)$ | 10 | 34430 | 18003 | 16427 | 0 |
| $(3,3)$ | $(4,1,1)$ | 8 | 26400 | 13063 | 13337 | 0 |
| $(4,4)$ | $(3,0,0)$ | 7 | 1896 | 825 | 1037 | 34 |

The first two rows saturate the cochain dimension because $Q_0^2=0$. For the last row the adjacent length-six and length-eight dimensions are $1169=344+825$ and $1588=1037+551$. These equalities also give the matching upper bounds on both ranks at length seven. Each finite-field determinant is a nonzero reduction of a characteristic-zero minor. The cochain dimensions are calculated separately from characters. A8 jobs 5708505, 5708370 and 5710346 produced the large minors for the first two rows. [Native odd-sector differential](../../dual_singlets_20261118/taylor_hooks/odd_colour_minor.cpp) · [Minors and seeds](outgoing_minors.json.gz) · [Replay](outgoing_channels.py)

Restricting every colour to a diagonal matrix annihilates classical boundaries. Include both charge parities at length seven. There are 3758 cochains, an incoming minor of size 1698 and a minor of size 2060 for the outgoing differential stacked with Cartan restriction. Thus the classical kernel of restriction has dimension at most $3758-1698-2060=0$. All 34 classes are detected on the Cartan. [Stacked matrix and second-prime replay](outgoing_channels.py)

## Products and the one-dimensional quotient

On the traceless Cartan, evaluate each graviton trace on the three diagonal entries $x,y,-x-y$ with the graded signs. At length seven all possible trace lengths are covered by lengths two through seven. In dotted weights $(4,4)$ and $(5,3)$ the ranks of the product evaluations are 176 and 143. Complete reducibility of the dotted $SL(2)$ gives singlet multiplicity $176-143=33$. The 33 archived singlets are constructed from traces of lengths two and three. Each is a proper product or a linear combination of proper products. The integer minor evaluating them is

$$
26631369071794128280428233087422233415424680686256128 .
$$

The calculation includes all 447 and 361 graviton recipes at the two weights. [Exact product expansion](outgoing_products.py) · [Recipes and integer minor](outgoing_products.json)

It remains to allow products with arbitrary finite-colour factors. A length-one trace vanishes in $SU(3)$, so the two factor lengths are $2+5$ or $3+4$. Their nonnegative flavour degrees split $(3,0,0)$. Their dotted spins must agree to contain a singlet. Charge conjugation acts on a diagonal word of length $Y$ by $(-1)^Y$. Factors of the opposite parity have zero Cartan image, and restriction is injective in the target.

These conditions give 56 pairs and 112 factor sectors with nonempty cochain spaces. For each factor the classical cochain dimension minus the two adjacent rank bounds equals the number of independent graviton products measured on the Cartan in its highest weight. All 112 bounds saturate. Both factors have nonzero cohomology for 26 of the pairs. Consequently every proper product lies in the 33-dimensional span above, while the 33 explicit singlets attain it. This proves $34-33=1$ for the quotient by all proper products. [All 114 differential matrices](outgoing_factor_minors.json.gz) · [Enumeration and native replay](outgoing_factors.py)

## Instanton interpretation

The vanishing rows exclude two charge-conjugation-odd outgoing channels of the level-30 operator R-007. The length-seven block is available classically. A leading instanton map satisfying the previously derived Cartan-vanishing Ward constraint on R-007 has zero projection to this block, since restriction here is injective. The renormalized physical instanton matrix element at the conformal origin remains open. [Ward hypotheses](completed_discriminant.md) · [Classical restriction to the two-colour block](levi_block.md)

## Reproduction

Run these commands in this directory with the repository's Python dependencies (NumPy, SciPy, SymPy and python-flint), a C++17 compiler and assertions enabled.

~~~sh
python3 -B outgoing_channels.py
python3 -B outgoing_channels.py --rebuild small --threads 4
python3 -B outgoing_factors.py --rebuild
python3 -B outgoing_products.py | diff -u outgoing_products.json -
~~~

The small native replay regenerates the four-pair endpoint. Use `--rebuild full` in its place to regenerate the large A8 minors as well.
