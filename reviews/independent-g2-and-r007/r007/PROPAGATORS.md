# The first direct propagator channels for the SU(3) level-30 Ward matrix

Take the bare 6,139-term R-007 source, the radial conjugate of Dot's 5,080-term target, and the improved chiral supercurrent of undeformed four-dimensional $\mathcal N=4$ Yang–Mills theory. At the conformal origin, for charge-one instantons of either orientation, allow Gaussian matter propagators between the endpoints and evaluate the separated current on the classical background. In this sector the self-dual orientation requires at least four direct propagators and the anti-self-dual orientation at least six. The calculation below sets the three- and five-propagator contributions respectively to zero before bosonic moduli integration. [Dot's target](https://github.com/lbotao399-stack/quantum-lifts-results/blob/04fe44e3a06f3fc589f72fe7bfd8568da5afc9a3/nonperturbative-lift/contributions/dot-r007-t2-primary/PROOF.md) · [Earlier bounds of three and five and field conventions](https://github.com/lbotao399-stack/quantum-lifts-results/blob/6a6f5a708ab1f73f00f8baafd9a133e824fe8238/algorithms/su3_j30_classical_ranks/instanton_probe/flavour_modes.md)

## Self-dual orientation

In flavour order $(4,1,2,3)$ the minimum occupations of the source, conjugate target and current are respectively $(0,3,6,6)$, $(9,5,2,2)$ and $(0,1,1,1)$. One same-species scalar or fermion propagator removes one complete quartet from their combined occupations. At exactly three propagators the lower bound is

$$
(0,3,6,6)+(9,5,2,2)+(0,1,1,1)-3(1,1,1,1)=(6,6,6,6).
$$

The full one-instanton exterior algebra has precisely six generators of each flavour. Every surviving elementary field and the current must therefore be taken at its minimum, and the quartic lifted-mode action contributes its constant term. All six remaining factors of flavour four must come from the six uncontracted letters at the radial target. This conclusion is unchanged when a derivative acts on a bosonic-background propagator. [Charge-one modes, orientations and ADHM projections, §§3 and C](https://arxiv.org/abs/hep-th/0212332)

All undotted oscillator indices of the ket singleton have the same polarization. In a jet with divided derivative labels $(p,q,I)$, where $I$ is its set of flavour derivatives, there are $p+q+|I|-1$ such indices. In the conjugate letter all dotted indices therefore lie on the same spinor line. In particular its covariant derivatives are the two $D_{\alpha\dot+}$ in a frame adapted to that line. The independent undotted indices on the conjugate play no role in this restriction. The receipt checks the occupation assignment for the 45,720 occurrences of 33 distinct jets. [Beisert, equations (4.9)–(4.14)](https://arxiv.org/abs/hep-th/0307015)

At the insertion choose one broken special supersymmetry of flavour four, combined with its Poincaré translation so its conformal Killing spinor vanishes there. Choose its dotted parameter on the same line as the upper dotted polarization. The antisymmetric two-spinor pairing on this line is zero. On the complex null plane through the insertion tangent to both polarized derivatives the conformal Killing spinor vanishes identically: it is the displacement bispinor contracted with that dotted parameter. Thus the variation of the connection and of each scalar vanishes on the plane. In the fermion variation the curvature term is multiplied by the vanishing Killing spinor. Its remaining scalar term has the chirality and dotted index of the special parameter. It vanishes on an undotted fermion and its contraction with the common dotted polarization vanishes on a dotted fermion. The polarized dotted curvature is the commutator of the two covariant derivatives tangent to the plane. Every iterated tangential covariant derivative of these fields obeys the same variation. This argument applies to open colour indices at the lowest field degree. A compensating gauge rotation of positive fermionic degree cannot alter that degree. [Dymarsky–Pestun, equations (2.4)–(2.6)](https://arxiv.org/abs/0911.1841) · [Broken superconformal orbit and elementary degrees, §2](https://arxiv.org/abs/1605.06346)

Equivalently, the two components of the local exact Poincaré combination of flavour four can occur, together with at most one independent extra component of its exact special parameter. Translation to the insertion is an invertible change of the four exact coordinates. The broken chiral Poincaré and antichiral special generators mutually anticommute on gauge invariants, and the leading covariant field variations just given fix the same statement for the open letters. The two additional ADHM modes of flavour four are independent of these exact coordinates. Consequently the uncontracted right factors use at most $2+1+2=5$ independent odd coordinates of flavour four. Their degree-six part vanishes. The exterior dimensions in degrees zero through six are $(1,5,10,10,5,1,0)$. Restoring the missing dotted component gives a nonzero exterior control. [Exact and lifted collective coordinates, Appendix C](https://arxiv.org/abs/hep-th/0212332)

## Anti-self-dual orientation

The gamma-traceless conserved current for the undotted supercharge has Lorentz spins $(1,\tfrac12)$. A putative contribution of fermionic degree one is made from a degree-zero anti-self-dual curvature and one dotted gaugino. Its Lorentz content is

$$
(0,1)\otimes(0,\tfrac12)=(0,\tfrac32)\oplus(0,\tfrac12).
$$

There is no component of undotted spin one. A single twice-differentiated adjoint gaugino has zero defining trace. The terms with a scalar, a derivative and an undotted matter fermion require at least $2+3=5$ modes. Two scalars and a dotted gaugino require at least $2+2+1=5$. Covariant derivatives with a bosonic background preserve these lower bounds. This accounts for the dimension-$7/2$ gauge-invariant current, including improvements with the stated spins. The current has the same R weight at every allowed degree, so the next allowed occupation is at least $(2,1,1,1)$. The three endpoints then require

$$
(10,7,4,4)+(0,4,7,7)+(2,1,1,1)=(12,12,12,12).
$$

Five propagators leave at least seven coordinates of each flavour. This exceeds six. The argument concerns a separated improved current. In a Ward identity its contact terms and fluctuations contracted with the current require independent operator matching. Neither of the new lower bounds determines a projected instanton differential.

## Replay

The input hashes are checked before reconstruction. Species multiplicities count individual pairings before signs and colour contractions. There are 7,929,439,422 such three-propagator pairings, of which 7,097,115,018 pass the separate leading source capacity bound. The polarization argument applies to them all.

~~~sh
python3 r007_polarized_propagators.py --out standalone.json
cmp standalone.json r007-polarized-propagators.json
~~~
