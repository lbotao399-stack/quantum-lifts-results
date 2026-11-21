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

## Perturbative survival of the two fermionic tower seeds

Section 4.2 of arXiv:2610.20340v2 reports the unrestricted calculation for the two $SU(3)$ operators $\mathcal O_{30}^{(3)}$ and $\mathcal O_{32}^{(3)}$. For each, the image under the first correction $Q_1$ is a boundary for the classical supercharge $Q_0$. Thus their first supercharge obstruction $d_1$ vanishes. The same paper proves first-order protection for their scalar dressings. The level-32 operator here is the even projection of the published second seed. Charge conjugation commutes with both supercharges, so projecting a primitive gives a primitive for this representative. Our unrestricted polynomial reconstruction remains partial. [Published first-order computation, §4.2](https://arxiv.org/html/2610.20340v2#S4.SS2) · [Seed construction, §2](https://arxiv.org/html/2412.08695v2#S2)

Under filtered Taylor equivariance, $d_2$ on the level-32 seed takes values in the even length-ten representation $(5,5,1,1,1,1)$. Its rational classical dimension and the rank of its restriction to commuting Cartan letters are

$$
6828-4185-2615=28,\qquad 155-270+221-97+19=28.
$$

Every perturbative differential has zero Cartan restriction by the discriminant argument. Injectivity therefore gives $d_2=0$. At level 30 the corresponding two counts are $3683-2174-1491=18$ and $110-198+163-72+15=18$. From order three onward a polynomial Taylor partner would require a negative flavour degree. At third order the three shifted degrees are $(3,-1,-1)$ at level 32 and $(2,-1,-1)$ at level 30. Each further order decreases all three entries by one. An incoming order at either seed would require a negative spacetime derivative degree. This proves perturbative survival of both seeds under the stated selection. [Cartan minors](../cartan_products/README.md) · [Discriminant argument and higher-order degrees](../../../perturbative_cartan_targets.md#three-second-order-target-components)

The published $F_1$ tower consists of nonzero symmetric flavour projections of products of these seeds with $u_2$ and $u_3$, the quadratic and cubic scalar gravitons. Take closed formal completions of the seed and graviton factors, form the quantum product and project to the same flavour representation. The product theorem gives a closed completion with the published classical leading term. Every tower member has zero spacetime derivative degrees, so the same incoming-degree exclusion applies. This removes the constructed $F_1$ tower from the perturbative search under filtered Taylor equivariance. The separate level-34 source of flavour $(8,5,4)$ and the nonperturbative action are independent questions. [Tower construction](https://arxiv.org/html/2412.08695v2#S2) · [Quantum product](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L3495)

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

## The six- and seven-scalar equations and the five-scalar restriction

Count occurrences of the undifferentiated first scalar $S=\phi^1$. A single rational first correction with 12,528 trace terms solves all 78 six-scalar and all 13 seven-scalar contents on arbitrary traceless three-colour matrices. Of these 91 residuals, 73 vanish as formal trace expressions. The remaining 18 vanish in all 153,356 collected polynomial coefficients after imposing the three-colour identities. There are 293 target contents in total. [Input corrections](scalar_layers/) · [Complete three-colour checks](scalar_layers/six_and_seven.json) · [Construction](../../su3_level32_layers.py)

A further 880-term change has no formal six-scalar variation. The combined correction has 13,400 terms. It solves 55 five-scalar evaluations on $S=\operatorname{diag}(1,1,-2)$ using 51 independent variation directions. Keeping the original right-hand side, direct reconstruction of the integer matrix gives the saved nonzero determinant and zero in all 55 equations. [Trace recipes and rational coefficients](scalar_layers/levi_solved.json) · [Independent reconstruction](../../su3_level32_levi.py) · [Determinant](scalar_layers/levi_replay.json)

This correction solves the entire five-scalar equation on the centralizer of $S$. In each of the other letters put an arbitrary $2\times2$ matrix $A$ in the upper block and $-\operatorname{tr}A$ in the lower block, with the parity of that letter. In a trace containing $m$ copies of $S$ and the ordered remaining matrices $A_1,\ldots,A_k$, the upper block gives $\operatorname{tr}(A_1\cdots A_k)$. The lower block gives $(-2)^m(-1)^k\operatorname{tr}(A_1)\cdots\operatorname{tr}(A_k)$. At $k=0$ their sum is $2+(-2)^m$. Keeping the order of odd entries, the 143,633 five-scalar residual terms reduce to 11,570 formal two-colour trace terms across the contents. Every coefficient vanishes on generic matrix entries in all 119 contents. Homogeneity restores every multiple of $S$. For unrestricted three-colour matrices one recorded coefficient has residue $884535$ modulo $1000003$. [Symbolic expansion and unrestricted witness](scalar_layers/centralizer_five.json) · [Four variables per matrix, with exterior signs](../../su3_level32_centralizer.py)

### A polynomial primitive on the diagonal scalar slice

Keep $S=\operatorname{diag}(x,y,-x-y)$ and restore **arbitrary** three-colour matrices for the other letters. Consider the part of $Q_0$ that splits off one copy of $S$. The 30 other letter species form 15 pairs. On one member of a pair this part is zero, and on the other it is the signed commutator of $S$ with the first. On the matrix entry $(i,j)$ the coefficient is the difference of the $i$th and $j$th diagonal entries of $S$. [Signed split table and pair construction](../../su3_level32_levi.py)

The five-scalar residual is closed under this part of $Q_0$. Indeed, the six-scalar component of $Q_0(Q_1\mathcal O_0+Q_0\mathcal O_1)=0$ has only two terms, and the term obtained from the six-scalar residual vanishes by the complete calculation above. Now fix the total occupation of each pair of matrix entries, and fix the diagonal entries' individual occupations. Each occupied off-diagonal pair has one even and one odd variable. At fixed positive total occupation it has two possible monomials. The differential between them is the corresponding eigenvalue difference times a nonzero rational number. The full block is the tensor product of these two-term complexes over $\mathbb Q[x,y]$.

An invertible rational change of basis in these Koszul factors retains just the independent linear forms, with zero differentials in the other factors. Two independent eigenvalue differences generate $(x,y)$, so the cohomology in positive scalar degree vanishes. For one independent difference the cohomology is detected by reducing modulo that difference. All occupied off-diagonal entries then belong to one $2\times2$ block. Its reduction vanishes by the full centralizer calculation above and permutation of the three colours. With no off-diagonal entries the residual vanishes on the full Cartan since $Q_1\mathcal O_0$ has zero Cartan restriction and $Q_0$ consists of commutators. Thus a four-scalar polynomial primitive exists on this diagonal slice, including the coincident-eigenvalue lines. Extending a primitive to a conjugation-invariant polynomial for arbitrary $S$ is a separate equation. [First-image Cartan calculation](order1.json) · [Five-scalar identity](scalar_layers/centralizer_five.json)

## Five-scalar equations on unrestricted matrices

A rational change with 238 trace terms and four copies of $S$ extends the 13,400-term first correction. The combined vector has 13,638 terms. Its maximum scalar count after applying $Q_0$ to the change is five, preserving all 91 higher-scalar equations.

Group the 119 five-scalar contents by removing $S$ and identifying the two members of each scalar-split pair. Three families with $5+8+5=18$ contents already vanish formally. In four further families the residual has two spectator traces whose scalar splits vanish. Factoring them off reduces the complete three-colour coefficient equations to systems with 452 or 669 parent words. Their ranks modulo 1,000,003 are 172 and 259. The returned rational corrections have $112+110+8+8=238$ terms after restoring the spectators and their Grassmann signs. Each cofactor identity was verified over the rationals; after restoring the spectators all 30 contents vanish, ten by collecting 84,696 integer-polynomial coordinates and twenty formally. For these invariant polynomials it suffices to diagonalize $S$ with generic eigenvalues and leave every other matrix arbitrary, by density and conjugation invariance.

The result on unrestricted matrices is $48/119$ at five scalars. The other 71 contents have nonzero evaluations for this particular vector, since the four families of added terms do not enter them. At six and seven scalars the result remains $91/91$. A8 job 5685041 independently tested all fourteen input families and reduced their formal traces. Its separate formal correction leaves 75 nonzero five-scalar residuals. A nonzero residual of a particular vector leaves the full $d_1$ equation open. [Rational change](scalar_layers/five_scalar_descent/delta.json.gz) · [Exact cofactor and polynomial checks](scalar_layers/five_scalar_descent/factored.json) · [A8 evaluations and payload hashes](scalar_layers/five_scalar_descent/a8_descent.json)

Rebuilding $Q_1\mathcal O_0+Q_0\mathcal O_1$ from the full 331,201-term first image gives the following results for the same rational vector. Every nonzero entry has a nonzero rational matrix evaluation.

| Copies of $S$ | Formal zeros | Further three-colour identities | Nonzero residuals |
|---:|---:|---:|---:|
| 3 | 3 | 0 | 11 |
| 4 | 21 | 0 | 48 |
| 5 | 38 | 10 | 71 |
| 6 | 60 | 18 | 0 |
| 7 | 13 | 0 | 0 |
| **Total** | **135** | **28** | **130** |

The vector satisfies $135+28=163$ of the 293 equations. [All 293 coefficient tests](scalar_layers/five_scalar_descent/all_contents.json) · [Full residual replay](../../su3_level32_full_residual.py)

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

Replay the successive scalar calculations with the three-colour engine built above.

```sh
l="$p/scalar_layers"
python3 "$d/su3_level32_layers.py" --pair "$l/component_pair/content_certificate.json.gz" --four "$l/component_four/content_certificate.json.gz" --engine /tmp/su3-colour-columns --work /tmp/su3-high-layers
cmp /tmp/su3-high-layers/six_and_seven.json "$l/six_and_seven.json"
python3 "$d/su3_level32_levi.py" --out /tmp/su3-five-scalar
cmp /tmp/su3-five-scalar/levi_correction.json.gz "$l/levi_correction.json.gz"
python3 "$d/su3_level32_centralizer.py" --correction /tmp/su3-five-scalar/levi_correction.json.gz --out /tmp/su3-five-scalar/centralizer.json
cmp /tmp/su3-five-scalar/centralizer.json "$l/centralizer_five.json"
```

The unrestricted five-scalar calculation uses the same three-colour exporter. From the repository root, set the algorithm module paths as follows (NumPy, SciPy and python-flint must be installed).

```sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_level32"
r=algorithms/su3_independent_prime_20261118/runtime/work
export PYTHONPATH="${d}:algorithms/su3_j34_20261118:algorithms/so7_cartan_target_20261118/native:$r/sp3-native.nosync:$r/cluster-next-v2.nosync/work/so8-bmn/vendor${PYTHONPATH:+:$PYTHONPATH}"
c++ -O3 -std=c++17 algorithms/su3_j34_20261118/colour_columns.cpp -o /tmp/su3-colour-columns
python3 "$d/su3_level32_scalar_descent.py" prepare --correction "$p/scalar_layers/levi_correction.json.gz" --folder /tmp/su3-five-families
cmp /tmp/su3-five-families/groups.json "$p/scalar_layers/five_scalar_descent/groups.json"
python3 "$d/su3_level32_factored_descent.py" --folder /tmp/su3-five-families --engine /tmp/su3-colour-columns --work /tmp/su3-five-factored
cmp /tmp/su3-five-factored/factored_delta.json.gz "$p/scalar_layers/five_scalar_descent/delta.json.gz"
cmp /tmp/su3-five-factored/factored_descent.json "$p/scalar_layers/five_scalar_descent/factored.json"
python3 "$d/su3_level32_full_residual.py" --delta "$p/scalar_layers/five_scalar_descent/delta.json.gz" --engine /tmp/su3-colour-columns --work /tmp/su3-all-contents
cmp /tmp/su3-all-contents/all_contents.json "$p/scalar_layers/five_scalar_descent/all_contents.json"
```
