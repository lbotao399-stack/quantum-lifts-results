# Exact zero-jet nonexactness criterion in the holomorphic Ward complex

**Audit result.** The proposed argument is valid for the stated perturbative holomorphic BF plus beta-gamma model with derivative-free cubic superpotential. It proves that a specified gauge-invariant scalar polynomial nonzero modulo the classical F-term ideal is closed and nonexact under the full perturbative differential. It does **not** define a natural ring injection from the classical F-term quotient into quantum cohomology.

The external input is the full local differential constructed by the source Ward/holomorphic-kernel formalism. In `/Users/libotao/np_lifts_engine/work/allloop_repo/main.tex`, the anchors are `eq-1-1`, `eq-C-local-Ward-identity`, `eq-C-current-delta-integral`, `eq-C-loop-count`, `eq-antifield-loss`, and the zero-form images at lines 1248–1254. The argument below independently verifies the grading and field-count consequences. It does not re-prove completeness of that perturbative construction or its equivalence to the full untwisted physical operator complex.

## 1. Every nonzero local n-loop kernel has derivative bidegree (n,n)

The graph $\Gamma$ has one local insertion, $n+1$ interaction vertices, and $2n+1$ Wick edges. Equivalently,

$$
|\Gamma_0|=n+2,
\qquad |\Gamma_1|=2n+1,
\qquad |\Gamma_1|-|\Gamma_0|+1=n.
$$

These equalities use the local heat-edge form-degree condition, not a classification into fans. For the source current representation, $\lambda_{v,+\dot a}$ are output momenta, $z_e^{+\dot a}$ generate derivatives on contracted input letters, and $y_{e,+\dot a}$ are edge currents.

**Proof by independent scaling of the two current components.** Let $r,t>0$ be numerical rescaling parameters. Apply

$$
\begin{array}{lll}
y_{e,+\dot+}\mapsto r y_{e,+\dot+},&
\lambda_{v,+\dot+}\mapsto r\lambda_{v,+\dot+},&
z_e^{+\dot+}\mapsto r^{-1}z_e^{+\dot+},\\
y_{e,+\dot-}\mapsto t y_{e,+\dot-},&
\lambda_{v,+\dot-}\mapsto t\lambda_{v,+\dot-},&
z_e^{+\dot-}\mapsto t^{-1}z_e^{+\dot-}.
\end{array}
$$

The loop-current equations and the positive-Schwinger-time domain are preserved. Each edge measure contributes $rt$; each independent vertex current constraint contributes $(rt)^{-1}$. There are $|\Gamma_1|$ edge measures and $|\Gamma_0|-1$ independent vertex constraints. Hence

$$
\mathcal I_\Gamma(r\lambda_{\dot+},t\lambda_{\dot-};
r^{-1}z_{\dot+},t^{-1}z_{\dot-})
=r^{|\Gamma_1|-|\Gamma_0|+1}
t^{|\Gamma_1|-|\Gamma_0|+1}
\mathcal I_\Gamma(\lambda;z)
=r^nt^n\mathcal I_\Gamma(\lambda;z).
$$

Holomorphic locality makes every displacement Taylor coefficient a polynomial in the output momenta. A coefficient with $a_+$ derivatives in $z_{\dot+}$ and $a_-$ in $z_{\dot-}$ therefore has momentum bidegree

$$
(n+a_+,n+a_-).
$$

The integers $a_+,a_-$ count derivatives on contracted input letters. Uncontracted spectator letters retain their original derivatives. Consequently, if $n_+,n_-$ are the total input holomorphic derivative counts, every nonzero output term of $Q_n$ has totals

$$
(n_++n,n_-+n).
$$

This includes derivative inputs: differentiating their displacement factors adds momentum powers; it does not remove them. Possible cancellations can make an entire homogeneous term zero but cannot produce a lower-degree term. The proof applies to every contributing Laman topology, including any nonfan or multivertex graph admitted by the local Ward formula.

## 2. The zero-jet projection is a chain map

Set every positive holomorphic jet to zero and keep the zero-form letters $\gamma_i,\beta^i,\mathsf b_A,\mathsf c^A$. This is an explicit evaluation map, distinct from the manuscript's homotopy-transfer projections. Derivative-free interactions imply that $Q_0$ preserves the two nonnegative derivative totals, so for every local cochain $B$,

$$
\left.Q_0B\right|_{\text{positive jets}=0}
=Q_0\left(\left.B\right|_{\text{positive jets}=0}\right).
$$

Section 1 gives, for every $n\ge1$,

$$
\left.Q_nB\right|_{\text{positive jets}=0}=0.
$$

Thus for the source $\mathbf Q=Q_0+\hbar Q_1+\hbar^2Q_2+\cdots$,

$$
\left.\mathbf QB\right|_{\text{positive jets}=0}
=Q_0\left(\left.B\right|_{\text{positive jets}=0}\right).
$$

The identity holds coefficientwise for formal perturbative series. It does not require a spectral-sequence convergence assumption.

## 3. Projecting the zero-jet complex to the F-term quotient

The source zero-jet images have the following exact structural content:

$$
\begin{aligned}
Q_0\mathsf c^A&=-\frac g2C^A{}_{BC}\mathsf c^B\mathsf c^C,\\
Q_0\gamma_i&=-ig\mathsf c^A(t_A)^j{}_i\gamma_j,\\
Q_0\beta^i&=ig\mathsf c^A(t_A)^i{}_j\beta^j
-\frac{\partial\bar W}{\partial\gamma_i},\\
Q_0\mathsf b_A&=-gC^B{}_{AC}\mathsf c^C\mathsf b_B
-ig\gamma_j(t_A)^j{}_i\beta^i.
\end{aligned}
$$

The matter slot $i$ and its dual use the source variance. Set $\mathsf c=\mathsf b=\beta=0$ after applying $Q_0$ to a zero-jet cochain $B$. Only the source F term can survive:

$$
\left.Q_0B\right|_{\mathsf c=\mathsf b=\beta=0}
=-\frac{\partial\bar W}{\partial\gamma_i}
\left.\frac{\partial_LB}{\partial\beta^i}
\right|_{\mathsf c=\mathsf b=\beta=0}.
$$

The derivative is left differentiation in the odd variable. All other source images still contain a discarded letter. In particular, the matter moment map in $Q_0\mathsf b$ retains one $\beta$.

Consequently, compose zero-jet evaluation with this evaluation and reduction modulo the ideal generated by $\partial\bar W/\partial\gamma_i$. Calling the resulting nonexactness detector $\pi$, its target is

$$
\mathbb C[\gamma]\big/
\left(\frac{\partial\bar W}{\partial\gamma_i}\right),
$$

with zero differential, and

$$
\pi\mathbf Q=0.
$$

This is a nonexactness detector. No assertion that this quotient embeds as a quantum ring has been used.

## 4. No outgoing loop action on a pure scalar polynomial

Every Wick edge, including the cut heat edge, consumes one $\mathsf b$ or $\beta$ terminal. Each gauge or matter interaction contributes one such terminal; a superpotential vertex contributes none. For an insertion containing only $\gamma$ fields, there are at most $n+1$ available terminals but $2n+1$ edges.

More explicitly, the source terminal count is

$$
\#(\mathsf b,\beta\text{ in the output})
=0+(n+1-\#(\bar W\text{ vertices}))-(2n+1)
=-n-\#(\bar W\text{ vertices})<0
$$

for $n\ge1$. Therefore no such Wick contraction exists, and

$$
Q_n f(\gamma)=0\qquad(n\ge1).
$$

If $f(\gamma)$ is gauge invariant, then also $Q_0f=0$. Thus the uncorrected polynomial itself obeys $\mathbf Qf=0$. This proof counts every interaction vertex and requires neither a fan assumption nor a single-superpotential-vertex assumption.

## 5. Exact survival statement, including the absence of hbar torsion

Suppose the specified gauge-invariant scalar polynomial satisfies

$$
\pi(f)\ne0.
$$

Then $\mathbf Qf=0$ by §4. If there were a perturbative cochain $B$ and a nonnegative integer $m$ such that

$$
\mathbf QB=\hbar^m f,
$$

applying the detector gives

$$
0=\pi(\mathbf QB)=\hbar^m\pi(f)\ne0.
$$

The last inequality holds in the formal-series target over $\mathbb C[[\hbar]]$. Hence $f$ is nonexact and not killed by a power of $\hbar$.

For the massless cubic $SU(5)$ benchmark, the source classical F-term/finite-rank calculation gives

$$
(\operatorname{Tr}\gamma^2)^2\ne0
\quad\text{modulo the classical F-term ideal}.
$$

The antichiral scalar representative $(\operatorname{Tr}\gamma^2)^2$ therefore survives to every perturbative order in this holomorphic complex. Conjugating back gives the corresponding statement for the source's holomorphic scalar class, provided the stipulated physical/twisted correspondence is applied.

## 6. Scope and the important non-conclusion

- The proof is for the polynomial local-operator complex with nonnegative holomorphic jet degrees and the source's complete perturbative local Ward differential. It does not assume an instanton coefficient, a nonperturbative counterterm prescription, or a nonperturbative duality map.
- Gauge-anomaly cancellation is needed for $\mathbf Q^2=0$. It is satisfied by the adjoint plus vectorlike fundamental matter in the stated benchmark.
- Disconnected vacuum or bulk-anomaly components do not evade the count: the stated QME/Ward construction removes components not connected to the local insertion. Every term of the local differential is a connected rooted collision graph; uncollapsed vertices only dress its effective representative.
- Cubic massless $s_1=0$ is sufficient. Formal derivative-free mass insertions do not alter the scaling count, but no claim about a mass-resummed propagator or an integrated-out nonpolynomial effective theory is needed here.
- Physical gaugino terms in the Konishi anomaly become positive ghost jets under the source field dictionary. The detector kills them. The projection is therefore a concrete reason the usual anomaly concern does not defeat this particular nonexactness criterion.
- Completeness of the holomorphic graph differential and its comparison with renormalized physical composite operators remain external construction/correspondence inputs. A finite redefinition need not preserve the literal coordinate projection; when a valid chain equivalence is supplied, the detector must be transported with it.

**There is no natural ring-injection claim.** If a scalar representative is zero modulo the classical F-term ideal, this detector says nothing about its quantum class. In particular,

$$
\mathbf Q(\beta^i f_i)
=Q_0(\beta^i f_i)+\hbar Q_1(\beta^i f_i)+\cdots
$$

can contain positive-jet Konishi terms. Classical F-term-zero representatives and products can therefore carry hidden quantum contributions. The argument proves survival of the specified square; it does not prove that the entire perturbative trace ring remains $\mathbb C[\operatorname{Tr}\gamma^2]/((\operatorname{Tr}\gamma^2)^3)$, nor that the cube is perturbatively exact.

Within these boundaries, the prior perturbative-survival gap for the square is closed by a chain-level argument, not by superpotential nonrenormalization alone. The nonperturbative Kutasov relation still requires the separate operator-duality/physical matching input, and no microscopic nonperturbative odd primitive has been constructed here.
