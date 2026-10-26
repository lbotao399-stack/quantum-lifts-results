# Fortuitous $Q$-cohomology and quantum lifting: computational results

Botao Li, October 2026. Companion to [All-loop-holomorphic-kernels](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels) (the operators $Q_n$ used below).

**Master report (updated in place):** [Fortuity_Lift_Master.pdf](Fortuity_Lift_Master.pdf). Each theory has seven items there: twisted action and components, Lie-algebra realization, highest-weight algorithm (which weights), frontier in conserved charges, fortuitous generators and classes, quantum-lift status, and per-level data.

$\color{red}{\text{Red}}$ = class not found in publications.

![Summary table](summary_table.png)

## Theories

| theory | complete frontier (all sectors) | fortuitous generators and classes | quantum lift | PDF |
|---|---|---|---|---|
| $\mathcal N=4$, $SU(2)$ | $E+J_L\le\frac{33}{2}$ | $\mathcal O_0$ and its multiplet, $\mathcal O_1$, products; $\color{red}{\text{P39}}$, $\color{red}{\text{X40}}$ | $\mathcal O_0$, P39, X40: all orders; $\mathcal O_1$: one loop | §1 |
| $\mathcal N=4$, $SU(3)$ | $E+J_L\le15$ | $\color{red}{X_7}$ multiplet, $F_1$, $B_1$ towers; $\color{red}{\text{new primaries}}$ $(5,4,4;1,1)Y10$, $(7,4,4;1,0)Y11$, $(6,4,4;1,1)Y11$, $(7,5,4;1,0)Y11$ | $X_7$, $F_1$ ($n=3$): all orders; $(7,4,4;1,0)$: all orders [cond]; $(5,4,4;1,1)$: $d_2$ running | §1 |
| $\mathcal N=4$, $SU(4)$ | $E+J_L\le\frac{27}{2}$ | $\color{red}{\text{first }SU(4)\text{ class}}$ $(4,2,2;2,2)Y9$ ($K=1$ over $\mathbb Q$) | all orders | §1 |
| $\mathcal N=4$, $SU(5..15)$ | $E+J_L\le14,15,15,\frac{31}{2},\frac{33}{2},17,18$ | none | – | §1 |
| $\mathcal N=4$, $SO(2N+1)$, $Sp(N)$ | $E+J_L\le\frac{23}{2}$ ($SO(5,7,9)$), $\le12$ ($Sp(3,4,5)$, $SO(11)$) | $\mathcal O_0$ ($SO(3)=Sp(1)$); $\color{red}{SO(5)=Sp(2)}$ $(3,2,2;1,1)Y7$; $O_f$ ($SO(7)$); $\color{red}{Sp(3)}$ $(5,2,2;1,1)Y9$, $\color{red}{\omega}$; $\color{red}{SO(9)}$ $(5,3,3;0,0)Y10$; $\color{red}{SO(11)}$ $(7,3,3;0,0)Y12$ | $SO(5)$ class: all orders; $O_f$, $SO(9)$ class: one loop; $Sp(3)$, $SO(11)$: running | §2 |
| $\mathcal N=2$ $\mathbb Z_2$ quiver, $SU(N)^2$ | $E+J_L\le10$ ($N=2$), $\le\frac{17}{2}$ ($N=3$) | fortuity split open | $d_2=0$ on all 44 candidate pairs; $E_r(\kappa)\cong E_r(1)$ | §3 |
| $\mathcal N=2$ SCQCD, $SU(N)$, $N_f=2N$ | $E+J_L\le\frac{23}{2}$ ($N=2$), $\le9$ ($N=3$) | fortuity split open | $SU(2)$: $d_1=0$; $d_2=0$ on all 41 candidate pairs | §4 |
| massive-adjoint flows | $w\le16$ ($\mathcal N=4$ $SU(2)$, one mass), $w\le12$ ($SU(3)$, one mass), $w\le13$ (SCQCD $SU(2)$), $w\le11$ (SCQCD $SU(3)$), $w\le13$ ($\mathbb Z_2$ quiver) | $\mathcal O_7$, $\mathcal O_{16}$ (Choi–Kim) | $\mathcal O_7$, $\mathcal O_{16}$: all orders; $d_{2,3}=0$ on every candidate pair | §5 |

## Conventions

$$\lambda=(\theta_1,\theta_2,\theta_3;n_+,n_-),\qquad \Theta=\theta_1+\theta_2+\theta_3,\qquad Y=\#\text{letters},\qquad \mathsf d=(1,1,1;-1,-1)$$

$$Q_n:\ C^{Y}(\lambda)\to C^{Y+1-n}(\lambda-n\mathsf d),\qquad d_r:\ E_r^{Y}(\lambda)\to E_r^{Y+1-r}(\lambda-r\mathsf d)$$

$$E=\tfrac12(\Theta+Y)+n_++n_-,\qquad R_I=\theta_I-\tfrac12(\Theta-Y),\qquad J_L=\tfrac12(\Theta-Y+n_++n_-),\qquad J_R=\tfrac12(n_+-n_-)$$

$$E=R_1+R_2+R_3+2J_L,\qquad (E,R_I,J_L,J_R)(\lambda+k\mathsf d,\,Y+k)=(E,R_I,J_L,J_R)(\lambda,Y)$$

Massive flows, $(A,B;n_\pm;\mathrm{gh})$ with $A=\theta_1+\theta_3$, $B=\theta_2+\theta_3$, $\mathrm{gh}=Y-\Theta$:

$$R_1+R_3=A+\mathrm{gh},\qquad R_2+R_3=B+\mathrm{gh},\qquad w=A+B+2(n_++n_-)=\tfrac43\left(E_{\mathrm{IR}}+J_L\right),\qquad E_{\mathrm{IR}}=\tfrac32R_{\mathrm{IR}}+2J_L,\qquad R_{\mathrm{IR}}=\tfrac12(R_1+R_2)+R_3$$

Certification: [Q] exact over $\mathbb Q$; [p] ranks modulo $p=2^{31}-1$; [p,p'] also $p'=2147483629$; [cond] $\mathrm{res}(R_r)=0$ at targets with $h=\rho$.

## Data

| file | content |
|---|---|
| [data/N4_SU2_levels.csv](data/N4_SU2_levels.csv) | $SU(2)$ per level: sectors, highest-weight classes, bosonic/fermionic states, non-graviton blocks |
| [data/N4_SU2_nongraviton_blocks.csv](data/N4_SU2_nongraviton_blocks.csv) | $SU(2)$ odd-$Y$ classes with charges and classification |
| [data/N4_SO_Sp_levels.csv](data/N4_SO_Sp_levels.csv) | $SO(3..11)$, $Sp(1..5)$ per level: sectors, $h$, $f$, fortuitous sectors |
| [data/N2_Z2_quiver_orbifold_point_levels.csv](data/N2_Z2_quiver_orbifold_point_levels.csv) | $\mathbb Z_2$ quiver per level: $h$, $d_1$, $E_2$, $d_2$ pairs and ranks |
| [data/N2_SCQCD_levels.csv](data/N2_SCQCD_levels.csv) | SCQCD per level |
| [data/massive_flows_levels.csv](data/massive_flows_levels.csv) | massive flows per $w$ |
