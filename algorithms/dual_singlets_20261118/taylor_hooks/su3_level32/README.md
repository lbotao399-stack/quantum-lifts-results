# The level-32 SU(3) BMN source

The published operator at flavour weight $(8,4,4)$ has eleven letters. Expanding the 47 summands in the printed v2 formula and taking its charge-conjugation-even component gives a highest-weight representative $\mathcal O_0$ with 8,640 rational trace terms. Its classical variation $Q_0\mathcal O_0$ has 26,267 formal trace terms. At three colours, all 761,859 matrix monomials in the 14 letter contents have zero coefficient after adding 1,010,262,640 integer contributions. [Published operator](https://arxiv.org/html/2412.08695v2#S2.SS1) · [Upstream trace conventions](https://github.com/davidsd/fortuity-ai/tree/4256f207986029540060fb44dd494746828f65ce/code/fortuity/bmn) · [Source expansion](../../su3_level32_source.py) · [Rational representative](O32_C_even.json.gz) · [Classical closure](classical.json)

On restriction to $SU(2)\times U(1)$ the coefficient of four abelian scalars equals $-1020600$ times the published $SU(2)$ threshold operator in every one of its 117 nonzero matrix coordinates. Applying the two raising charges and projecting to the level-36 highest weight $(8,5,5)$ gives an 82,088-term descendant. Applying the two previously archived block functionals to this descendant and the level-34 descendant gives

$$
\det\begin{pmatrix}-163296000&0\\816480000&400982400\end{pmatrix}
=-65478821990400000 .
$$

[Two-colour boundary test](../su3_bmn_sources/quartic.json) · [Direct descendant calculation](descendant.json)

The coefficients $Q_1$ and $Q_2$ of the supercharge give the following direct images at doubled charges $(4,4,16,8,8)$. Generic traceless diagonal substitution gives zero in both cases.

| Cochain | Length | Formal trace terms | Common denominator | Data |
|---|---:|---:|---:|---|
| $Q_1\mathcal O_0$ | 11 | 331,201 | 3 | [Image](Q1.jsonl.gz) · [Digest and control](order1.json) |
| $Q_2\mathcal O_0$ | 10 | 31,347 | 1 | [Image](Q2.jsonl.gz) · [Digest and control](order2.json) |

In the ten-letter classical Taylor representation with rows $(5,5,1,1,1,1)$, the cochain dimension minus the two adjacent image ranks is

$$
6828-4185-2615=28\quad(\text{even}),\qquad
6639-4019-2620=0\quad(\text{odd})
$$

in both fields of 2,147,483,647 and 1,000,033 elements. A8 jobs 5684039 and 5684216 used different seeds. All five native weights give identical adjacent ranks at the two primes. The source has even charge conjugation. The first equation over the full three-colour space is $Q_0\mathcal O_1=-Q_1\mathcal O_0$, where $\mathcal O_1$ is the first correction. On a solution the second obstruction is the class of $Q_2\mathcal O_0+Q_1\mathcal O_1$ after quotienting classical boundaries and first-differential images. [First prime](second_target/) · [Independent prime](second_target_prime1000033/) · [Taylor decomposition](../README.md)

## Three scalar restrictions of the first equation

Set one scalar to $x\operatorname{diag}(1,1,-2)$ and put every other Taylor letter in the upper traceless two-by-two block. Here $x$ is its commuting scalar coordinate. The classical commutators with the chosen scalar vanish. The three choices of scalar give simultaneous necessary equations for the same $SU(3)$ correction. The sign of the archived $Q_1\mathcal O_0$ is reversed before solving.

At length ten there are 14,088 distinct products of five quadratic traces at the required weight. An exact linear combination of 180 matrix coefficients pairs to $-1/3$ with $-Q_1\mathcal O_0$. It annihilates the classical variation of every one of these products. Exactly 8,086 of them have a restriction on the scalar loci and powers entering the functional. This excludes any first correction made entirely of products of quadratic traces. [Functional](first_correction/joint.json) · [Rational verification](first_correction/verified_pairs.json)

For a scalar $S$ and an ordered active pair of letters $U,V$ containing a derivative, the needed additional polynomial can be chosen to contain

$$
-\frac12\operatorname{tr}(S^2UV)+\frac14\operatorname{tr}(S^2)\operatorname{tr}(UV).
$$

On the locus for $S$ it equals $x^2\operatorname{tr}(UV)$. On either other locus, $S$ is traceless in the upper block, where $S^2=\tfrac12\operatorname{tr}(S^2)\mathbf1_2$. At least one of $U,V$ is a derivative letter with zero lower entry. The polynomial therefore vanishes there. Multiplication by powers of $\operatorname{tr}(S^2)/6$ supplies higher even powers of $x$.

Adding 2,399 such combinations gives a 1,602-term rational correction with 1,444 products of five quadratic traces and 158 products of one quartic and three quadratic traces. A fresh expansion on generic traceless two-by-two matrices gives zero residual in every coefficient on all three loci. On the archived noncommuting three-colour assignment the residual $Q_0\mathcal O_1+Q_1\mathcal O_0$ is $-5497700040000$. This particular correction requires additional terms in the full space. [Shared correction](first_correction/quartic.json) · [Independent matrix replay](first_correction/verified_quartic.json) · [Three-colour test](first_correction/three_colour.json) · [A8 component solves](first_correction/receipt.json)

Exchanging the two derivative directions and exchanging the two equally weighted flavours both act by $-1$ on all 331,201 first-image terms. Each exchange commutes with all 180 classical letter splits. Taking the negative eigenspace of both exchanges gives a 2,612-term correction satisfying the same three scalar restrictions. Its variation cancels the previous residual content, which has only one derivative-carrying letter of bidegree $(1,1)$. A different full three-colour coefficient is $1077753600000$. [Projected coordinates and symmetry checks](first_correction/projected.json) · [Scalar restrictions](first_correction/verified_projected.json) · [Full three-colour coefficient](first_correction/projected_three_colour.json)

At second order a separate three-colour evaluation on five Grassmann generators gives

$$
(Q_0Q_2,\ Q_1Q_1,\ Q_2Q_0)\mathcal O_0
\longmapsto (-428558902632000,\ 428558902632000,\ 0).
$$

The bare $Q_2\mathcal O_0$ fails classical closure. The last entry is zero also by letter content. Three output letters in this evaluation carry derivatives, while $Q_2$ on underived inputs can place derivatives on only two output letters. [Integer evaluations and matrices](second_consistency.json)

## A three-colour functional and a simultaneous correction

Choose the three traceless colour matrices in the saved functional. The trace of every pair and all 27 traces of a commutator times a third matrix vanish. On the selected four-fermion coefficient the 46 surviving trace words in $Q_1\mathcal O_0$ give $-9676800000$. Direct evaluation of all 5,589 incident charge-conjugation-even parent trace orbits gives zero on the 2,499 orbits containing a quadratic trace. The 1,774 nonzero variations have trace lengths $(3,3,4),(3,7),(4,6),(5,5)$ or $(10)$. Thus a correction written as a sum of trace products requires at least one summand without any quadratic trace factor. [Matrices and exterior coefficient](contents_pair/trace_functional.json) · [Parent variations and examples](contents_pair/trace_parents.json) · [Complete enumeration](../../su3_level32_trace_parents.py)

For this letter content and its derivative-exchanged partner, eight parent contents give 1,011 even finite-colour coordinates against 266 target coordinates. The classical image and the augmented matrix both have rank 264 modulo 1,000,003. A8 job 5684433 reconstructed a shared 493-term rational correction. Expanding its classical variation on arbitrary traceless three-colour matrices cancels $Q_1\mathcal O_0$ in all 2,030 monomials in each content. There are 293 allowed target contents in the full first equation. [Correction and exact polynomial checks](contents_pair/content_certificate.json.gz) · [Frame and receipt](contents_pair/receipt.json) · [Letter contents and incidence](contents_pair/graph.json)

## Replay

Run from the repository root with the dependencies of the archived SU(3) runtime. The two JSONL digests in the records refer to the uncompressed streams.

```sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_level32"
c++ -std=c++17 -O3 algorithms/su3_j34_20261118/source_closure.cpp -o /tmp/su3-level32-closure
python3 "$d/su3_level32_replay.py" classical --engine /tmp/su3-level32-closure --out /tmp/su3-level32-classical.json
cmp /tmp/su3-level32-classical.json "$p/classical.json"
python3 "$d/su3_level32_replay.py" descendant --out /tmp/su3-level32-descendant.json
cmp /tmp/su3-level32-descendant.json "$p/descendant.json"
python3 algorithms/su3_independent_prime_20261118/build_native.py
python3 "$d/su3_level32_actions.py" --order 1 --out /tmp/su3-level32-actions
python3 "$d/su3_level32_actions.py" --order 2 --out /tmp/su3-level32-actions
c++ -std=c++17 -O3 -pthread -I "$d/taylor_hooks/native" "$d/taylor_hooks/unitary_parity.cpp" -o /tmp/su3-level32-parity
python3 "$d/unitary_hook_parity.py" --engine /tmp/su3-level32-parity --spin 5 --tail 4 --length 10 --prime 2147483647 --seed 783311 --out /tmp/su3-level32-target
cmp /tmp/su3-level32-target/resolved.json "$p/second_target/resolved.json"
```

To reconstruct the source from the published TeX, pass its v2 main file to `su3_level32_source.py --paper`, the pinned Fortuity repository's `code` directory to `--fortuity-code`, and an output gzip path to `--out`. The program checks the hashes of the TeX and the two imported upstream modules.

The three scalar tests can be regenerated with the same dependencies. The first two verification commands below use the checked-in correction and functional, respectively.

```sh
python3 "$d/su3_level32_first.py" extract --work /tmp/su3-level32-walls
python3 "$d/su3_level32_walls.py" verify --work /tmp/su3-level32-walls --certificate "$p/first_correction/quartic.json" --out /tmp/su3-level32-quartic.json
python3 "$d/su3_level32_walls.py" verify-pairs --work /tmp/su3-level32-walls --certificate "$p/first_correction/joint.json" --out /tmp/su3-level32-pairs.json
python3 "$d/su3_level32_walls.py" pairs --work /tmp/su3-level32-walls --certificate /tmp/su3-level32-functional.json
python3 "$d/su3_level32_walls.py" quartic --work /tmp/su3-level32-walls --certificate /tmp/su3-level32-correction.json
python3 "$d/su3_level32_residual.py" --certificate /tmp/su3-level32-correction.json --out /tmp/su3-level32-residual.json
python3 "$d/su3_level32_symmetry.py" --certificate /tmp/su3-level32-correction.json --out /tmp/su3-level32-projected.json
r=algorithms/su3_independent_prime_20261118/runtime/work
PYTHONPATH="$r/sp3-native.nosync:$r/cluster-next-v2.nosync/work/so8-bmn/vendor" python3 algorithms/su3_j34_20261118/second_consistency.py --level32 --out /tmp/su3-level32-nilpotency.json
```

The two-content calculation can be checked in a fresh directory. The first command regenerates the right-hand sides from the archived full first variation. The colour exporter uses generic traceless matrices and integer Grassmann coefficients.

```sh
c++ -O3 -std=c++17 algorithms/su3_j34_20261118/colour_columns.cpp -o /tmp/su3-colour-columns
python3 "$d/su3_level32_trace_parents.py" --certificate "$p/contents_pair/trace_functional.json" --out /tmp/su3-trace-parents.json
cmp /tmp/su3-trace-parents.json "$p/contents_pair/trace_parents.json"
python3 "$d/su3_level32_content.py" prepare --work /tmp/su3-content-pair
cp "$p/contents_pair/content_certificate.json.gz" /tmp/su3-content-pair/
python3 "$d/su3_level32_content.py" verify --work /tmp/su3-content-pair --engine /tmp/su3-colour-columns
```
