# Novelty evidence: exact primary locators

**Correction after a broader branch search:** the SO13 open-problem entry is stale. Existing `reviews/orthogonal-mixed-lifts/bmn13` at `78d170a26afc2a16cae43f123df86d5880be6196` proves nonzero products of the proposed scalar parents. See `SOURCE_LINEAGE.md`; the snippets below must not be used alone as evidence that those classes are new.

## SO13 is predicted, with direct cohomology unresolved

Repository: `davidsd/fortuity-ai`.
Pin: `4256f207986029540060fb44dd494746828f65ce`.
Path: `summary/fortuity-summary.tex`.
Blob: `7a6f2e1d5c21831db92f089dae2870a200f1852c`.

- Lines 526–531 give the first index difference `-56 t^30+168 t^32`, the candidate content `[6,0]+[2,2]+[0,0]`, and the rank-pattern expression `-Sym^(N-3)[2,0]` at `j=4N+6`. The text qualifies the three-point character check as **“strong evidence, not a proof”**.
- Line 1009 lists **“direct SO(13) cohomology at j=30 (three fermionic classes predicted)”** among the open tasks; the displayed quotation converts only the source's TeX group notation to plain text.
- The complete-scan table at lines 576–582 has SO7 through j18, SO9 through j22 and SO11 through j24. There is no SO13 row.

Independent data locator at the same pin: `data/so_sp/so13_sp6_bmn_index.meta.json`, blob `8f45bc7e817f9d41bd7c7512d777c956d223e083`. Its only SO13 refined-content test is evaluation at three fugacity points; a full refined character expansion was not computed. This is index evidence, not an actual cohomology representative or a minimal-generator proof.

## The SO11 `(5,5,3)` sector is already constructed

Same upstream pin. Path: `data/so_sp/so11_nongraviton_rep_1553.meta.json`.
Blob: `61a95e639edeec089ea9b8860530647e78401294`.

The recorded highest-weight sector is `(2J;q)=(1;5,5,3)` with SU3_R label `[0,2]`. Its dimensions are `2354` total, `1820` closed, `1808` exact, `12` cohomology, `11` graviton, `1` non-graviton. The supplied representative passed fresh SO11 closure tests and has nonzero two-block restriction with coefficients in the ratio `1:-2:1`. This excludes it as a new class, despite its omission from the companion's short census table.

## Sp4/Sp5 are absent from the current manuscript census

Repository: `lbotao399-stack/All-loop-holomorphic-kernels`.
Pin: `bba7a4857710c1baef433341f1428b63ba7642f0`.
Path: `main.tex`.
Blob: `028e76b765d84f73985b46d5af64c8b84415ace4`.

- Line 474 has the census row `Sp(4), Sp(5)`, complete cutoff `24`, with **“no fortuitous class in scan range”**.
- Lines 579–585 contain the Sp2 and Sp3 thresholds and the known Sp3 sector `(5,0,0;3,3)_9`.
- Line 593's rank formula concerns `SO(2N+1)` and explicitly restricts its established cases to `N=3,4,5`.
- The full fetched manuscript contains no universal Sp(N) threshold or omega conjecture. This is a bounded source inspection, not a global absence theorem.

The candidate Sp4 `(7,0,0;3,3)_11`, j32, and Sp5 `(9,0,0;3,3)_13`, j36, do not appear in these inspected sources. Their existence remains a separate computation; two successful ranks may belong to one common family.

## SU2 has two known singlets at j36

Upstream pin `4256f207986029540060fb44dd494746828f65ce`.
Path: `data/su2/bmn_hw_cohomology_j44.json`.
Blob: `332f24c8c86a750005c60121592df804855ffb25`.

At `lam=[6,6,6]`, `twoJ=9`, the record gives `m=25,h=2,grav=0,BH=2`. The summary's SU2 one-loop discussion identifies the two known singlets as O1 and O1′, with different lifting behavior. Re-extracting the second singlet would not be a new generator.

## Scope of the new local length argument

`NOVELTY_BASELINE.md` proves same-group BMN product/raising exclusion for a fortuitous SO13 class of length14 and j30 using only specialization injectivity through length13. It neither assumes nor claims a completed SO13 scan through j28. A full-letter minimal-level statement covering the first relation module at j29 is an additional theorem and has not been established by this catalog.
