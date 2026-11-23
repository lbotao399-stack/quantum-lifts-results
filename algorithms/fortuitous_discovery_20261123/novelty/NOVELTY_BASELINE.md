# Fortuitous-class novelty baseline — 2026-11-23

**Critical later source-lineage correction:** the upstream SO13 open-problem statement below is stale relative to `research/codex-independent-mixed-lifts-20261122` at `78d170a26afc2a16cae43f123df86d5880be6196`. Its `reviews/orthogonal-mixed-lifts/bmn13/README.md` already proves two independent scalar-product classes, forcing the SO13 j30 parents `(5,5,5)` and `(7,5,3)` to be nonzero. A fresh exact comparison with that complete quotient also finds the degree17 all-singleton `(7,5,5)` relation nonzero at SO13, so its persistence to SO15 excludes that backup under the stated novelty criterion. See `SOURCE_LINEAGE.md` and `source_lineage/` for the actual circuits and rank receipts. These lanes must not be counted as new fortuitous classes.

## Conventions / 记号沿用 source

`j=2(E+J_L)`; `Y` is letter length. Paper tuples are `(theta1,theta2,theta3;n+,n-)_Y`; native `trace_descent` tuples put `(n+,n-,theta1,theta2,theta3)` first. BMN means `n+=n-=0`; its existing labels are `(2J;q1,q2,q3)` with `Y=q1+q2+q3-2J`. No renormalization is used here.

## Read-only source pins

| Source | Exact pin | What was inspected |
|---|---|---|
| quantum-lifts-results latest main | `dd99912081c05cf9ee6507df4e3244e32e7ec821` | `RESEARCH_DISCUSSION.md`, complete tree, comparison to provided source pin |
| Same repository runtime / historical tables | `082be7070d207e2cdecfd0a91b03d5e014fd303d` | `supplement/supplement.tex`, `data/N4_SO_Sp_levels.csv`, runtime code |
| Same repository detailed branch | `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9` | 26 selected research Markdown files, 42 campaign ledgers, complete tree, native Taylor engine |
| Latest USp(4) existence bound linked by main | `d4608bf0270b7692dff90f154241ac44db9538f0` | `algorithms/dual_singlets_20261118/bmn_primary/second_odd_generator.md` |
| Primary Fortuity repository HEAD | `4256f207986029540060fb44dd494746828f65ce` | `summary/fortuity-summary.tex`, complete tree, SO11 representative metadata, SO13 index metadata |
| Current main paper source | `bba7a4857710c1baef433341f1428b63ba7642f0` | Entire `All-loop-holomorphic-kernels/main.tex`; blob `028e76b765d84f73985b46d5af64c8b84415ace4` |

The latest main tree does not contain the detailed `dual_singlets_20261118` branch. Its research log links immutable branches. Absence from the main tree alone is therefore not novelty evidence. All tool actions were read-only remotely.

## Immediate candidate priorities

| Priority | Exact proposed source | Novelty baseline | Available route and remaining requirement |
|---|---|---|---|
| 1 | `Sp(4)=USp(8)`, native `(3,3,7,0,0)`, `Y=11`, `j=32` | Not present in the inspected catalogs or current paper. The paper states no fortuitous class through `j=24`; historical data have only 12/116 sectors completed at `j=28`. No all-rank Sp threshold/omega theorem or conjecture was found. | `trace_descent.dimensions/highest_weights/action`; native engine command uses paper order `chain Sp4 7 0 0 3 3 11 11`. Need full adjacent ranks, actual finite closure and quotient by stable cycles. Distinguish from the known Sp3 omega and its Taylor multiplet. |
| 2 | `Sp(5)=USp(10)`, native `(3,3,9,0,0)`, `Y=13`, `j=36` | Same apparent catalog gap; no general Sp-family claim in the current paper. This is a separate rank/charge class, but two successful Sp cases might constitute one general family rather than two unrelated families. | Same route; native order `chain Sp5 9 0 0 3 3 13 13`. Larger exact cochain spaces require sparse/native evaluation. |
| Excluded after lineage check | `SO(13)`, `(5,5,5;0,0)_14`, `j=30`; `(7,5,3;0,0)_14` | The upstream index prediction and open-problem entry are superseded by an existing two-product cohomology certificate on another research branch. Its independent nonzero products force these parents to be nonzero. | New direct certificates improve explicitness but do not supply new classes. Same-group BMN primitiveness remains true, yet does not establish novelty. |

The scalar-agent dimensions and positive-agent new rank computations are not evidence used in this catalog. This note establishes prior-art exclusions and available code, rather than reviewing either exploratory implementation.

## Known classical-group classes that must be excluded

| Group | Existing class / sector | Existing status and counting restriction | Exact source |
|---|---|---|---|
| SU2 = SO3 = Sp1 | `O0: (4,4,4;0,0)_7`, `j=24`; `O1` and `O1′: (6,6,6;0,0)_9`, `j=36` | Published cores. The upstream exact census has two non-graviton singlets at `(2J;q)=(9;6,6,6)`; `O1` lifts at one loop and `O1′` survives. Their descendants and graviton products are not additional generators. | `supplement/supplement.tex` at `082be707`, section `sec-E-SU-classes`; upstream summary section SU2 and `data/su2/bmn_hw_cohomology_j44.json` at `4256f207` |
| SU2 | `P39: (5,5,5;3,0)_9`; `X40: (7,5,5;1,1)_10` | Already recorded full-letter primaries. P39 quotient is one new direction beyond 12 products/descendants. | Same supplement |
| SU3 | `X7: (2,2,2;2,2)_7`, `j=24`; member `(0,0,0;5,5)_7`, `j=30` | One known multiplet, 35 recorded highest-weight blocks; do not recount its members. | Same supplement |
| SU3 | `F1: (n+4,4,4;0,0)_(n+7)`, `j=24+2n`, `n>=3` | Published brane-fused tower; first source is R-007 at `(7,4,4)_10`, `j=30`. | Same supplement; `algorithms/su3_j30_classical_ranks`; arXiv 2412.08695 and 2610.20340 |
| SU3 | `B1: (n+5,5,4;0,0)_(n+8)`, `j=28+2n`, `n>=3` | Published bosonic tower; first level34 sector has an independently identified non-scalar-fused core direction. | Same supplement; upstream summary `sec:su3core`; `algorithms/su3_j34_20261118/README.md` |
| SU3 | `(5,4,4;1,1)_10` j32; `(7,4,4;1,0)_11` j33; `(6,4,4;1,1)_11` j34; `(7,5,4;1,0)_11` j35 | Already catalogued non-BMN primaries/sectors. Equal charges alone do not establish multiple primitive families. | Same supplement, `sec-E-SU-classes` |
| SU4 | `(4,2,2;2,2)_9`, j28 | Known full-letter fortuitous class; its absence in the BMN census is irrelevant because it is outside BMN. | Same supplement and current main paper |
| SU4 | BMN `(10,4,4)_13`, j36, `[6,0]` | Already constructed by Kim–Lee–Lee–Song; generalized all-N and multiple-core constructions are also published. | arXiv 2610.20340, 2610.21592, 2611.04722 |
| SO5 = Sp2 | `(3,2,2;1,1)_7`, j20 | Known full-letter threshold. | Same supplement/current paper |
| Sp3 | `(5,2,2;1,1)_9`, j24; omega `(5,0,0;3,3)_9`, j28 | Same known Taylor-related family. Latest main gives full rational omega fortuity quotient `1198-116-1081=1`; stable H is zero. Same grading vanishes at all Sp(N), N>=4. | `sp3_omega/README.md` at `9a554a6d`; `DUAL-CHATGPT-20261123` at `dd999120` |
| SO7 | `(3,3,3;0,0)_8`, j18 | Known O_f; one-loop-lifted. j20/22/24 blocks are descendants; j26 `(5,4,4)_10` is known graviton dressing; j28 two classes are its descendants. | Upstream summary `sec:so`, `data/so_sp/so7_direct_cohomology_1333.json`; GLRT 2506.13887 |
| SO9 | `(5,3,3;0,0)_10`, j22 | Explicit known fused threshold, nonzero SO7 block restriction. | Upstream `data/so_sp/so9_nongraviton_rep_1533.json` |
| SO11 | `(7,3,3;0,0)_12` AND `(5,5,3;0,0)_12`, j26 | Both explicit known classes, `[4,0]` and `[0,2]`; the second is omitted by the companion's abbreviated table. | Upstream `data/so_sp/so11_nongraviton_rep_1733.json`, `so11_nongraviton_rep_1553.json`; summary `sec:so` |
| SO8 | j18 `(3,3,3)_8`, reflection-even/odd triality doublet | Published prior art. The even member is trace-covering fortuitous; the odd member belongs to the generalized monotone Pfaffian family. Count interpretation depends on covering and triality/parity convention. | `so8_level18/README.md` at `9a554a6d`; Miao–Tang 2611.23173 sections6–8 |
| SO8 | j24 `(4,4,4)_10`, reflection-even/odd pair | Already rationally constructed and classified beyond lower audited BMN products and descendants. One triality standard representation, two parity components. | `so8_level24/README.md`, `generator.json`, `chain.json` at `9a554a6d` |
| Sp2 = USp4 | j30 `(5,5,5)_10` | Already named as the ten-letter singlet core. Initial campaign's scoped 3-family count includes the early BMN cores; it is not a global minimal-ring census. | `b2_level36/README.md` at `9a554a6d`, section Selection of higher-order targets |
| Sp2 = USp4 | j36 `(10,4,4)_13`: at least two independent BMN generators | Existing rational rank/index/product/parent proof, even without choosing representatives. Extracting two representatives would not discover two new classes. | `b2_level36/README.md`, `parents/result.json`, ledger `...20261121_primitive.json` at `9a554a6d` |

Also exclude all already constructed orthogonal glueball sources: SO9 j30 source into S5; its j36 product with S; SO10 j36; SO11 j36 d3 source; SO12 j42 d3 source; SO13 j42 d3 source; the SO13 j48 source is the j42 source multiplied by S; SO14 j48 is already recorded in D014. The SO17/SO18 determinant sources at j54/j60 already die by d1. These are known sources or product consequences; a new nonzero pairing is not automatically a new primitive generator.

The aligned SO11 closed combination at j36 has zero classical fortuitous component in latest main: `9*(1/6)+(25600/25200)*(-189/128)=0`. This excludes that proposed novelty lane, without claiming its entire finite leading class is zero.

## Existing existence bounds and zero sectors

- USp4: at least four odd BMN generators by j48 is already proved at `d4608bf0`. The six full odd columns `(6,5,4)`, `(7,4,4)`, `(7,5,3)`, `(8,5,2)` at j30 and `(6,5,5)`, `(7,5,4)` at j32 are zero over Q. Merely rerunning these columns cannot add a class.
- G2: an odd generator at `(6,6,6)`, j36, at Y=11,13 or15, and at least two odd generators by j44 are already proved as existence bounds. The `(6,5,4)` and `(6,6,4)` odd columns are zero. The j30 `(5,5,5)` cohomology is exhausted by four known products, with `3030-1960-1066=4` and `1070-1066=4`.
- G2 `det Tr(gamma^I gamma^J)` and S3 are non-Coulomb monotones. Cartan invisibility does not establish fortuity.
- Exceptional scalar additions in the campaign ledger are not a global classical-family fortuity count. The ledger itself explicitly disclaims exhaustive full-letter primitiveness and absolute literature novelty.

## SO13: exactly what was predicted, and what remains new

The following describes only the older upstream summary. The later branch and exact product-ancestry argument in the opening correction supersede its novelty conclusion.

Primary upstream `summary/fortuity-summary.tex`, pin `4256f207`, gives first SO13/USp12 index difference `-56 t^30+168 t^32` and candidate SU3_R content `[6,0]+[2,2]+[0,0]`. Its `data/so_sp/so13_sp6_bmn_index.meta.json` says the refined content was tested at three sampled fugacities and calls this “strong evidence, not a proof”. The open-problem list says “direct SO(13) cohomology at j=30 (three fermionic classes predicted)”. The complete tree contains only SO13 index files, not SO13 class representatives or a completed cohomology scan.

The three predicted weights are `(9,3,3)`, `(7,5,3)`, `(5,5,5)`, all Y14 and 2J1. Upstream interprets the whole rank pattern through `Sym^(N-3)[2,0]`: the SO7 core with one scalar bilinear per separated SO2 block. Thus a certified SO13 representative is new direct mathematics in an anticipated brane-fused family. Calling it an entirely new general mechanism would overstate novelty.

### Same-group BMN primitive test that does not need a full lower census

For a proposed fortuitous SO13 class at j30,Y14:

1. A nonconstant trace invariant has at least two letters. Every factor in a positive-level product of total length14 has length at most12.
2. Metric pairings on at most26 vector slots are linearly independent in defining dimension13: assign distinct defining indices to the13 pairs of a selected pairing. The selected pairing evaluates to1 and every other pairing to0. Antisymmetrizing the two indices in each adjoint slot and projecting identical letters preserve injectivity. For SO13, one epsilon cannot contract an even total number of vector indices; two epsilons reduce to metrics.
3. Consequently the finite specialization is injective on both the Y-letter factor and its Y+1-letter classical image for every factor with Y<=12. A finite closed factor therefore has a universally closed lift. Its products are monotone.
4. BMN raising operations preserve Y and raise j by2. A lower ancestor of Y14,j30 has j<=28. At j28,Y14 it is purely scalar and has a universally closed lift; at j<=26 no such cochain exists. Its raising descendants remain monotone.

Therefore any nonzero fortuitous SO13 class actually constructed in these Y14 sectors is independent of all same-group lower positive products and BMN raising descendants. This conclusion does not use index positivity or assume the existence of a complete SO13 j28 scan. It does not exclude a lower-rank brane-fused interpretation, or establish a full all-letter superconformal-generator classification.

The companion supplement separately states the SO(2n+1) fortuity-free stable range `Y<2n+2`. For n6, together with the pure-scalar Y14 endpoint, that statement would imply no lower BMN fortuity through j28. The product/descendant argument above uses only the strictly injective Y<=12 range and avoids relying on the borderline Y13 interpretation of that statement. A blanket full-letter no-fortuity claim through j29 would additionally require a proved injectivity statement for the first relation module under restriction to a larger Cartan; this catalog does not establish that endpoint claim. No complete SO13 cohomology enumeration through j28 was found in the upstream source.

## Sp4/Sp5: check against the current paper

The entire current `main.tex` at `bba7a485` was inspected. The census table records Sp4/Sp5 cutoff j24 and “no fortuitous class in scan range”. The threshold table specifies only Sp2 and Sp3, including Sp3 omega. The only adjacent general threshold formula is explicitly for SO(2N+1), N=3,4,5. There is no stated universal Sp(N) omega formula or conjecture in this paper, the detailed omega note, or the primary Fortuity summary.

Thus the Sp4/Sp5 sectors are absent from the inspected baseline. This is a bounded literature statement, not a proof that no unpublished calculation exists. A rank extension of omega still requires proving a new finite class and distinguishing whether the two ranks are being counted as two classes or as one generator-family pattern.

## Runnable source locations

`trace_descent.py` at `9a554a6d` accepts SO defining dimension and Sp defining dimension: `SO13`, `Sp8`, `Sp10`. Its `--flavours` argument parametrizes a full charge orbit, not a single native weight. For isolated sources use `dimensions`, `highest_weights`, `basis`, `action` directly, or the native engine.

The native engine uses rank for symplectic names (`Sp4`, `Sp5`) and paper tuple order `(theta1,theta2,theta3,n+,n-)`. Its exact fetched files are in:

`catalog/pinned_9a554a6d/algorithms/dual_singlets_20261118/taylor_hooks/native/`

The engine source blob is `a1e4b90432dfeaa02ec059f1905037730cbee172`; `core.h` is `3c2b4ebdffbb193c1800650bc9c8b23cd7201044`. These nine files were fetched read-only in one batch. The engine exposes `dims`, `chain`, `d1`, `selftest`; a full `chain` returns adjacent Q0 ranks and a Cartan rank. Cartan rank alone is not a fortuity test.

## Primary literature anchors

- [Miao–Tang, 2611.23173v1](https://arxiv.org/html/2611.23173v1), especially sections6–8: SO8 j18 triality pair and all-rank generalized Pfaffian monotones.
- [Kibe–van Leuven–Mathieson, 2610.22332](https://arxiv.org/abs/2610.22332): SO7 S4 and a known lifting partner; non-Coulomb does not imply fortuitous.
- [Kim–Lee–Lee–Song, 2610.20340](https://arxiv.org/abs/2610.20340): SU3/SU4 BMN families, including SU4 j36.
- [de Mello Koch–Kim–Liu, 2610.21592](https://arxiv.org/abs/2610.21592): existence of fortuitous BMN operators at every SU(N), N>=3.
- [de Mello Koch–Duary–Kim–Liu, 2611.04722](https://arxiv.org/abs/2611.04722): multiple lower-rank core constructions, including SU4.
- [Gadde–Lee–Raj–Tomar, 2506.13887](https://doi.org/10.1007/JHEP02(2026)026): known SO7 core and SO/Sp index comparison.

The arXiv search records and the fully opened 2611.23173v1 text were checked. Several other arXiv full-text URLs were unavailable through the web reader, so their detailed sector claims here are cross-supported by the fetched primary repository and current manuscript source, not by inferred unseen passages.
