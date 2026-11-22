# Nonperturbative lifting: monopoles and worldsheet instantons

Primary-source check, 2026-11-13. 本文仅讨论 local operators；没有把 monopole disorder operators 当作 line operators。

## 1. Minimal (0,2) CP1: the unit becomes exact

无 left-moving fermions。Yagi 的 $Q=\bar Q_+$ 是 physical right-moving supercharge；$[Q,O\}$ 是对 bosonic $O$ 取 commutator、对 fermionic $O$ 取 anticommutator。Chiral algebra 是 local-operator $Q$-cohomology，带 holomorphic OPE；它不是 ordinary topological field theory。

在 target coordinate $\phi$ 中，$g_{\phi\bar\phi}$ 和 $R_{\phi\bar\phi}$ 分别是 metric、Ricci tensor，$\bar\psi_+$ 是 right-moving fermion。原文 eq. (4.5) 的代表元为

$$
\theta=-\frac12R_{\phi\bar\phi}\partial_z\phi\,\bar\psi_+.
$$

Degree-one instanton $\phi_0:\mathbb P^1\to\mathbb P^1$ 的 fermion bundle 是 $K_\Sigma^{1/2}\otimes\phi_0^*TX$；$K_\Sigma^{1/2}$ 是 spin bundle，$\mathcal O(n)$ 是 degree-$n$ line bundle，$h^0,h^1$ 是其 holomorphic-section 和 first-cohomology dimensions：

$$
K_\Sigma^{1/2}\otimes\phi_0^*TX
=\mathcal O(-1)\otimes\mathcal O(2)=\mathcal O(1),
\qquad h^0(\mathcal O(1))=2,\quad h^1(\mathcal O(1))=0.
$$

两个 fermion zero modes 由 supercurrent 和 $\theta$ 饱和。计算要保留 zero-mode wavefunctions 对 fluctuating boson 的依赖，才能得到 supercurrent OPE 的 antiholomorphic pole；随后 instanton-moduli integral 降为 constant-map integral。

$\Lambda$ 是 dynamical scale，$t_0$ 是在该 scale 评价的 degree-one topological action。原文 eq. (4.2) 的 $O_k$ 是 degree-$k$ correction 的 local operator：

$$
\{Q,\theta\}
=\sum_{\substack{k\ge1\\\mathrm{instanton\ degree}}}
\Lambda^{2k}e^{-kt_0}O_k,
\qquad
\dim O_k=(1-k,0),\qquad R(O_k)=2-2k.
$$

不存在 negative-dimension local operators，故

$$
O_k=0\ (k>1),\qquad
\{Q,\theta\}=\Lambda^2e^{-t_0}O_1.
$$

$O_1$ 是非零 constant operator。原文 eqs. (4.15)–(4.16) 没有计算 overall numerical coefficient；不能在保留原 $\theta$ normalization 时把右边直接改成 $1$。[Yagi, *Chiral Algebras of (0,2) Models*, §§2.2, 4.1, eqs. (4.2), (4.5), (4.15)–(4.16)](https://arxiv.org/pdf/1001.0118).

以下是直接 algebraic consequence。固定非零 $\Lambda^2e^{-t_0}$，对 $[Q,O\}=0$，由 graded derivation 和 $Q^2=0$：

$$
[Q,\{Q,\theta\}^{-1}]=0,
$$

$$
[Q,\theta\{Q,\theta\}^{-1}O\}
=\{Q,\theta\}\{Q,\theta\}^{-1}O
-\theta[Q,\{Q,\theta\}^{-1}O\}
=O.
$$

这里 inverse 只作用于 nonzero constant operator，multiplication 取原 chiral algebra 的 renormalized local product；没有引入新的 $\theta$ normalization。

若保留 nonnegative instanton powers，即 coefficient ring 是 $\mathbb C[[e^{-t_0}]]$，则上述 primitive 使用 $e^{t_0}$，不属于该 ring。此时严格得到

$$
[Q,\theta O_1^{-1}O\}=\Lambda^2e^{-t_0}O,
\qquad e^{-t_0}[O]=0.
$$

若 formal $[1]=0$，在 $e^{-t_0}=0$ specialization 后会得到 perturbative $[1]=0$，与已知 perturbative unit class 矛盾。因此 formal torsion 与 fixed-coupling collapse 必须分开。

对于 Hilbert-space states，unit exactness 通过 operator action 使 $Q$-closed states exact；在适当定义的 normalizable state space 上意味着没有 supersymmetric states。Original treatment 在 cylinder/plane 比较 physical 与 twisted theory，再 untwist 后 compactify；不能忽略 curved-worldsheet twist anomaly。[Tan–Yagi, *Chiral Algebras of (0,2) Sigma Models: Beyond Perturbation Theory*, §§4.2.1–4.2.2, eq. (4.29)](https://arxiv.org/pdf/0801.4782).

## 2. Aharony duality: selected local monopoles are removed by F-terms

Electric $U(N_c)$ 有 $N_f>N_c$ fundamental/antifundamental pairs；magnetic theory 是 $U(N_f-N_c)$。$V_\pm$ 是 electric local monopoles，对应 magnetic elementary singlets；$\widetilde V_\pm$ 是 magnetic local monopoles。$M_i{}^{\tilde i}$ 是 meson singlet，$q^i,\widetilde q_{\tilde i}$ 是 magnetic quarks；$i,\tilde i$ 是两个 flavor spaces 的 component labels。原文 eq. (3.3)：

$$
W=M_i{}^{\tilde i}q^i\widetilde q_{\tilde i}
+V_+\widetilde V_-+V_-\widetilde V_+.
$$

仅对 elementary singlets 求导：

$$
\frac{\partial W}{\partial V_+}=\widetilde V_-,\qquad
\frac{\partial W}{\partial V_-}=\widetilde V_+,
\qquad
\frac{\partial W}{\partial M_i{}^{\tilde i}}=q^i\widetilde q_{\tilde i}.
$$

因此 quantum chiral ring 中

$$
[\widetilde V_-]=[\widetilde V_+]=[q^i\widetilde q_{\tilde i}]=0.
$$

这里是 ordinary local chiral ring；未计算 single-supercharge cohomology 的所有 spin/derivative sectors。这是 F-term ideal 去掉特定 classes；没有得到 $[1]=0$。不能把 composite $\widetilde V_\pm$ 当作 elementary fields 随意求导，进而误删 electric $V_\pm$。[Aharony, *IR Duality in d=3 N=2 Supersymmetric USp(2Nc) and U(Nc) Gauge Theories*, §3, eq. (3.3)](https://arxiv.org/pdf/hep-th/9703215).

Microscopic monopole-instanton generation 的直接例子是 pure $SU(2)$：两个 gluino zero modes 产生原文 eq. (6.2) 的 $W=1/Y$。作者明确 suppressed dimensional scale；所以该公式不能用于报告 absolute instanton coefficient。在该 source convention 内，$Y$ 是 Coulomb-branch chiral coordinate：

$$
\partial_Y W=-Y^{-2}\ne0
\qquad(0<|Y|<\infty).
$$

这是 finite-field supersymmetric vacuum 消失/runaway，不能直接替换成一个全局 UV local-operator cohomology vanishing theorem。

同一文的 $U(1)$ one-flavor mirror 使用 source normalization $W=-MV_+V_-$；三者在 dual Wess–Zumino description 中都是 elementary fields：

$$
\partial_MW=-V_+V_-,\qquad
\partial_{V_+}W=-MV_-,\qquad
\partial_{V_-}W=-MV_+.
$$

故三个 pairwise products 为零，三个 coordinate branches 仍在。Abelian SQED mirror 的这个 interaction 不应误称为原 $U(1)$ theory 的 smooth microscopic monopole-instanton amplitude。[Aharony–Hanany–Intriligator–Seiberg–Strassler, *Aspects of N=2 Supersymmetric Gauge Theories in Three Dimensions*, eqs. (3.2), (6.2); nonabelian generalization eqs. (8.3)–(8.4)](https://arxiv.org/pdf/hep-th/9703110).

## 3. (2,2) CP(N−1): nonperturbative deformation without collapse

Hori–Vafa 的 $\Sigma$ 是 $U(1)$ field-strength twisted chiral superfield；$Y_i$ 是 dual periodic twisted chiral fields，$i=1,\ldots,N$ 标记 charge-one matter species；$\mu$ 是 scale，$t$ 是 FI–Theta parameter。原文 eqs. (3.38), (5.1)：

$$
\widetilde W=
\Sigma\left(\sum_{\substack{i=1,\ldots,N\\\mathrm{matter\ species}}}Y_i-t\right)
+\mu\sum_{\substack{i=1,\ldots,N\\\mathrm{matter\ species}}}e^{-Y_i}.
$$

Exponentials 是 vortex-generated twisted F-terms；phase dualization、two-zero-mode correlation function、holomorphy 是其计算 ingredients。以下直接消去 $Y_i$，没有另定义 quantum parameter：

$$
\partial_{Y_i}\widetilde W=\Sigma-\mu e^{-Y_i}=0,
\qquad
\partial_\Sigma\widetilde W=
\sum_{\substack{i=1,\ldots,N\\\mathrm{matter\ species}}}Y_i-t=0,
$$

$$
\Sigma^N
=\mu^N\prod_{i=1}^{N}e^{-Y_i}
=\mu^Ne^{-t}=\Lambda^N,
\qquad\Lambda=\mu e^{-t/N}.
$$

$\sigma$ 是 $\Sigma$ 的 lowest component；$\ell$ 标记 distinct vacua：

$$
\sigma=\Lambda e^{2\pi i\ell/N},
\qquad \ell=0,\ldots,N-1.
$$

对 $N=2$，原文 mirror coordinate $\Theta$ 给出

$$
\widetilde W=\Lambda(e^\Theta+e^{-\Theta}),\qquad
\partial_\Theta\widetilde W=\Lambda(e^\Theta-e^{-\Theta})=0,
$$

$$
e^\Theta=1\ \text{or}\ -1,\qquad
\partial_\Theta^2\widetilde W=2\Lambda\ \text{or}\ -2\Lambda\ne0.
$$

两 isolated vacua；local A-model quantum cohomology/twisted-chiral ring 被 deformation，没有整个 unit 消失。[Hori–Vafa, *Mirror Symmetry*, §§3.2.2, 5.3, eqs. (3.38), (5.8)–(5.9)](https://arxiv.org/pdf/hep-th/0002222).

## Scope comparison for the parent synthesis

| Model | Verified nonperturbative statement | What it does not establish |
|---|---|---|
| Minimal (0,2) CP1 | At fixed nonzero instanton factor, every local $Q$-class is exact | Formal power-series cohomology is zero; (2,2) CP1 has the same behavior |
| Aharony duality | Magnetic monopole classes and magnetic mesons vanish by elementary-singlet F-terms | All local cohomology vanishes; local monopoles are lines |
| Pure 3d N=2 SU2 | Monopole instanton lifts the finite Coulomb-branch vacua | Global UV unit-exactness theorem with a computed absolute coefficient |
| Hori–Vafa CP(N−1) | Vortex terms produce $\sigma^N=\Lambda^N$ and $N$ massive vacua | Instantons must lift every BPS class |

CP1 minimal (0,2) 是 lifting mechanism benchmark，本身不是这里给出的 target-to-target duality theorem。以上任何例子都不能单独证明当前 4d project 的 instanton differential；需要该 project 自己的 zero-mode saturation、operator insertion、nonzero matrix element 和 coefficient-ring specification。
