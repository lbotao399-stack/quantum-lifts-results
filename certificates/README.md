# Exact rational (`[Q, cert]`) lifting certificates

Each `.json.gz` file stores an exact rational lifting chain `x = (x_0, ..., x_{r-1})` and a dual rational witness `phi = (phi_1, ..., phi_r)` (with rational entries `[num, den]`) certifying `d_r[x_0] != 0` over $\mathbb{Q}$ via three polynomial identities (Supplement Section S2.4):

- **(V1)** Every monomial of `x_k` lies in the stage-`k` block `B(k)`, and $\sum_{a+b=k} Q_a x_b = 0$ over $\mathbb{Q}$ for $k = 0, \dots, r - 1$.
- **(V2)** For every monomial $y \in B(m)$ ($m = 1, \dots, r$), $\sum_{k=m}^r \phi_k(Q_{k-m} y) = 0$ over $\mathbb{Q}$.
- **(V3)** The target $t = \sum_{k=1}^r Q_k x_{r-k}$ satisfies $Q_0 t = 0$ and $\phi_r(t) = 1 \ne 0$ over $\mathbb{Q}$.

## Included certificates and verifiers

| File | Description |
|---|---|
| `singlet_cubic_self_16_1_1_m8_d4_cert.json.gz` | **Four-loop lift $d_4 \ne 0$** (`r = 4`) in the abelian cubic-singlet theory (`singlet_cubic_self`, $W = s x y + \tfrac13 s^3$) at block `[16, 1, 1, -8, [0], [0]]` ($R = 8/3, J_L = 5, J_R = 0$), with two-term source $x_0 = (\mathsf b_0^5 s - 5\,\beta^x\beta^y\mathsf b_0^4)\,\partial_2\mathsf c_0\,\partial_1\mathsf c_0$; reconstructed across 11 primes (`max_height_bits = 251`) and verified over $\mathbb{Q}$ across all `106,587` stage monomials |
| `singlet_cubic_self_13_1_1_m6_d3_cert.json.gz` | **Three-loop lift $d_3 \ne 0$** (`r = 3`) in `singlet_cubic_self` at block `[13, 1, 1, -6, [0], [0]]` ($R = 8/3, J_L = 4, J_R = 0$); verified over $\mathbb{Q}$ across all `22,418` stage monomials |
| `n2su2_10_0_0_m5_m2_d2_cert.json.gz` | **Non-abelian two-loop lift $d_2 \ne 0$** (`r = 2`) in pure $\mathcal N=2$ $SU(2)$ Yang–Mills (`n2su2`) at block `[10, 0, 0, -5, -2]` ($R = 1, J_L = 5/2, J_R = 0$); verified over $\mathbb{Q}$ across all `193` $H=0$ stage monomials |
| `verify_abelian_cert.py`, `qengine.py` | Standalone Python 3 verifier (standard library only, `fractions.Fraction` for `PRIME=0`) for abelian $d_r$ certificates |
| `spanp_dual.cpp` | Sparse $\mathbb{F}_p$ Gaussian elimination solver with reverse pivot back-substitution (`mode=2`) for extracting sparse dual witnesses $\phi$ |
| `cert_dr.py` | Reference $SU(2)$-invariant subcomplex certificate driver and `PRIME=0` verifier |

## Standalone exact-$\mathbb{Q}$ verification (`PRIME=0`)

Requires only Python 3 (standard library):

```sh
# Verify the three-loop (r = 3) certificate at weight 13 over Q (~6 s):
PRIME=0 python3 verify_abelian_cert.py singlet_cubic_self '[13,1,1,-6,[0],[0]]' 3 singlet_cubic_self_13_1_1_m6_d3_cert.json.gz

# Verify the four-loop (r = 4) certificate at weight 16 over Q (~55 s, 106,587 stage monomials):
PRIME=0 python3 verify_abelian_cert.py singlet_cubic_self '[16,1,1,-8,[0],[0]]' 4 singlet_cubic_self_16_1_1_m8_d4_cert.json.gz
```

## Exact rank runs of 12 November 2026 (`exact_runs_2026-11-12/`)

Records behind the entries marked $[\mathbb Q]$ in Section 2.3 of the article (Tables of the Leigh–Strassler lifts and of the lifts at two, three and four loops). Every run is an exact elimination over $\mathbb Q$ (`PRIME=0`) of the complete block; the reported ranks are therefore rational ranks.

| File | Content |
|---|---|
| `LS_P14_exact.log` | $SU(3)$ Leigh–Strassler theory, source $(7,0;0,0)_5$: $\dim E_1=5$, $\operatorname{rank}d_1=3$, $\operatorname{rank}d_2=1$, with the rational certificate of the two-loop lift (field $\mathbb Q$, all checks passed) |
| `LS_dP14_exact.log` | the same theory, source $(7,0;1,0)_5$: rational certificate of the two-loop lift (field $\mathbb Q$, all checks passed) |
| `exact_<theory>_<block>.json.gz`, `abelian_exact_block_runs.log` | abelian theories (massive SQED, quartic SQED, cubic singlet, $\gamma_s^2\gamma_+\gamma_-$, $U(1)^2$ triangle and its mass deformation, $\mathcal N=2$ SQED with $N_f=1$): dimensions, $E_1$, and the ranks of $d_2$, $d_3$ of the block named in the file |
| `su2_exact_rank_runs.log` | pure $\mathcal N=2$ $SU(2)$ Yang–Mills theory and $SU(2)$ with a flipping singlet: exact $d_2$ ranks of nine blocks |

