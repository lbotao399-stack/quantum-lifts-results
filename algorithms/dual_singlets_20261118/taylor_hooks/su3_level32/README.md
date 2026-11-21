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

over the field of 2,147,483,647 elements. A8 job 5684039 evaluated five native weights. The source has even charge conjugation. The first equation over the full three-colour space is $Q_0\mathcal O_1=-Q_1\mathcal O_0$, where $\mathcal O_1$ is the first correction. On a solution the second obstruction is the class of $Q_2\mathcal O_0+Q_1\mathcal O_1$ after quotienting classical boundaries and first-differential images. [Five weights and Slurm receipt](second_target/) · [Taylor decomposition](../README.md)

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
