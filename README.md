# Fortuitous Q-cohomology and quantum lifting: computational results

Botao Li, October 2026. Companion to [All-loop-holomorphic-kernels](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels) (the operators Q_n used below).

**Master report (updated in place): [Fortuity_Lift_Master.pdf](Fortuity_Lift_Master.pdf).** Per theory: twisted action and components, Lie-algebra realization, highest-weight algorithm (which weights), frontier in conserved charges, fortuitous generators and classes, quantum-lift status, per-level data.

Red = class not found in publications.

![Summary table](summary_table.png)

## Theories

| theory | complete frontier (all sectors) | fortuitous generators and classes | quantum lift | PDF |
|---|---|---|---|---|
| N=4, SU(2) | E+J_L ≤ 33/2 | O_0 and its multiplet, O_1, products; $\color{red}{\text{P39}}$, $\color{red}{\text{X40}}$ | O_0, P39, X40: all orders; O_1: one loop | §1 |
| N=4, SU(3) | E+J_L ≤ 15 | $\color{red}{\text{X7 multiplet}}$, F1 and B1 towers; $\color{red}{\text{new primaries}}$ (5,4,4;1,1) Y10, (7,4,4;1,0) Y11, (6,4,4;1,1) Y11, (7,5,4;1,0) Y11 | X7, F1 (n=3): all orders; (7,4,4;1,0): all orders [cond]; (5,4,4;1,1): d_2 running | §1 |
| N=4, SU(4) | E+J_L ≤ 27/2 | $\color{red}{\text{first SU(4) class}}$ (4,2,2;2,2) Y9, K=1 over Q | all orders | §1 |
| N=4, SU(5) … SU(15) | E+J_L ≤ 14, 15, 15, 31/2, 33/2, 17, 18 (N=5…10, 11–15) | none | – | §1 |
| N=4, SO(2N+1), Sp(N) | E+J_L ≤ 23/2 (SO(5), SO(7), SO(9)); ≤ 12 (Sp(3), Sp(4), Sp(5), SO(11)) | O_0 (SO(3)=Sp(1)); $\color{red}{\text{SO(5)=Sp(2)}}$ (3,2,2;1,1) Y7; O_f (SO(7)); $\color{red}{\text{Sp(3)}}$ (5,2,2;1,1) Y9 and $\color{red}{\omega}$; $\color{red}{\text{SO(9)}}$ (5,3,3;0,0) Y10; $\color{red}{\text{SO(11)}}$ (7,3,3;0,0) Y12 | SO(5) class: all orders; O_f and the SO(9) class: one loop; Sp(3), SO(11): running | §2 |
| N=2 Z2 quiver, SU(N)×SU(N) | E+J_L ≤ 10 (N=2), ≤ 9 (N=3) | fortuity split open | d_2 = 0 on every candidate pair; E_r(κ) ≅ E_r(1) | §3 |
| N=2 SCQCD, SU(N), N_f=2N | E+J_L ≤ 12 (N=2), ≤ 9 (N=3) | fortuity split open | SU(2): d_1 = 0; d_2 = 0 on every candidate pair | §4 |
| N=4 SU(3), β and Leigh–Strassler deformations | E+J_L ≤ 9 | $\color{red}{\text{fortuitous}}$ (Θ,c;n+,n−)_Y = (7,0;0,0)_5, (7,0;1,0)_5 (LS) | $\color{red}{\text{LS: rank }d_2=1}$ at (7,0;0,0)_5, (8,0;0,0)_6, (7,0;1,0)_5; β: d_2 = 0 on 11 pairs | §6 |
| massive-adjoint flows | w ≤ 16 (N=4 SU(2), one mass), w ≤ 12 (SU(3), one mass), w ≤ 13 (SCQCD SU(2) and Z2 quiver), w ≤ 11 (SCQCD SU(3)) | O_7, O_16 (Choi–Kim) | O_7, O_16: all orders; d_2 = d_3 = 0 on every candidate pair | §5 |

## Conventions

```
λ = (θ1,θ2,θ3; n+,n−),   Θ = θ1+θ2+θ3,   Y = number of letters,   d = (1,1,1; −1,−1)

Q_n : C^Y(λ) → C^{Y+1−n}(λ − n d),        d_r : E_r^Y(λ) → E_r^{Y+1−r}(λ − r d)

E   = (Θ+Y)/2 + n+ + n−          R_I = θ_I − (Θ−Y)/2
J_L = (Θ − Y + n+ + n−)/2        J_R = (n+ − n−)/2
E   = R1 + R2 + R3 + 2 J_L       (λ + k d, Y + k) has the same (E, R_I, J_L, J_R)

massive flows, (A,B; n±; gh):  A = θ1+θ3,  B = θ2+θ3,  gh = Y − Θ
R1+R3 = A+gh,   R2+R3 = B+gh,   R_IR = (R1+R2)/2 + R3,   E_IR = (3/2) R_IR + 2 J_L
w = A + B + 2(n+ + n−) = (4/3)(E_IR + J_L)
```

Certification: [Q] exact over Q; [p] ranks modulo p = 2^31 − 1; [p,p'] also p' = 2147483629; [cond] res(R_r) = 0 at targets with h = ρ.

## Data

| file | content |
|---|---|
| [data/N4_SU2_levels.csv](data/N4_SU2_levels.csv) | SU(2) per level: sectors, highest-weight classes, bosonic/fermionic states, non-graviton blocks |
| [data/N4_SU2_nongraviton_blocks.csv](data/N4_SU2_nongraviton_blocks.csv) | SU(2) odd-Y classes with charges and classification |
| [data/N4_SO_Sp_levels.csv](data/N4_SO_Sp_levels.csv) | SO(3…11), Sp(1…5) per level: sectors, h, f, fortuitous sectors |
| [data/N2_Z2_quiver_orbifold_point_levels.csv](data/N2_Z2_quiver_orbifold_point_levels.csv) | Z2 quiver per level: h, d_1, E_2, d_2 pairs and ranks |
| [data/N2_SCQCD_levels.csv](data/N2_SCQCD_levels.csv) | SCQCD per level |
| [data/massive_flows_levels.csv](data/massive_flows_levels.csv) | massive flows per w |
