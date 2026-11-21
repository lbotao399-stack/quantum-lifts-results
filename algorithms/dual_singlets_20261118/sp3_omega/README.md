# The level-28 incoming differential into the USp(6) class ω

At energy 13 and left spin 1 the level is $j=2(E+J_L)=28$. The class $\omega$ has nine letters and Taylor degrees $(3,3,5,0,0)$, in the order of the two spacetime derivatives followed by the three odd-coordinate derivatives. The [earlier census](../../../Fortuity_Lift_Master.pdf) identifies the candidate in §2.6.

## Rational classical cohomology

There are 2,624 formal highest-weight predecessors of length eight. Exact trace-word elimination gives 1,672 independent images under the classical supercharge $Q_0$ and 952 independent integer relations. The corresponding 952 closed polynomials have a minor of determinant 348046 modulo 1000003 on the 2,559-dimensional finite-colour source. Every relation vanishes before imposing the finite-colour identities. The target boundary minor has order 1,607 and determinant 624328 at the same prime. Consequently

$$
1607\leq \operatorname{rank}_{\mathbb Q}Q_0\leq2559-952=1607.
$$

The nine-letter target has dimension 2,689. Its outgoing minor has order 1,081 and determinant 110010 modulo 1000003. The saved 819 rational source coefficients, of common denominator 13,824, give a polynomial with 4,264 terms. Its classical variation has 4,502 formal terms and vanishes on all 1,281 independent integer coordinates in each of two frames. Its value on the modular functional annihilating the incoming images is one. The rank equality above therefore makes the augmented boundary rank 1,608 over the rationals. This proves

$$
\dim_{\mathbb Q}H(Q_0)=2689-1607-1081=1 .
$$

The corresponding $SO(7)$ dimension is $2692-1608-1084=0$. A8 job 5681842 reproduced the complete integer cycle list and the finite-colour minor. [Classical matrices and frame definitions](target.json.gz) · [Predecessor frame](previous.json.gz) · [Exact relations and minor](boundary_bounds.json.gz) · [Rational cycle](cycle.json.gz) · [Independent integer replay](cycle_replay.json.gz) · [Orthogonal comparison](so7_target.json.gz)

## Lowering to the known nine-letter class

Replace one derivative of each spin by the second flavour, then one derivative of each spin by the third flavour. On all 2,870 formal highest-weight traces the four reverse substitutions in reverse order return exactly 72 times the input. On the representative of $\omega$, the lowered polynomial has 23,208 terms. Its first nonzero flavour-raising variation has 12,569 terms. The extremal projector onto the highest weight $(1,1,5,2,2)$ gives 63,508 terms. Its pairing with the lower target functional is 6 modulo 1000003. The reverse substitutions followed by the original highest-weight projector pair to 72 in the original block. The substitutions and weight projectors respect $Q_0$. A rational boundary at the lower weight would return to a rational boundary at the original weight, which the nonzero pairing excludes. Thus the lower cohomology also has rational dimension one, attaining the bound $7294-4336-2957=1$. [Formal substitutions](substitutions.json) · [Lower frame](lowered.json.gz) · [Projected representative checks](cycle_replay.json.gz)

The lowered first-order source has rational dimension $1861-712-1149=0$. An order $r$ incoming differential requires derivative bidegree $(1-r,1-r)$ there, which is empty for $r\geq2$. Transferring this exclusion back to $\omega$ on the second page requires a filtered quantum Ward action for the four odd generators. [Lowered source](lowered_source.json.gz) · [Selection rule and its scope](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex)

## Direct quantum actions

At prime 1000003 the chosen target functional pairs nontrivially with 4,488 of 25,147 formal first-action inputs. Its pairing with the second action vanishes on all 38,592 formal inputs. A8 job 5681845 found a linear combination of unrestricted first-action inputs equal to the target modulo classical boundaries. A cochain-level first-action preimage must also pass the classical closure equation to give a $d_1$ pairing. The earlier independent calculation at prime 2147483647 gives the first-order source dimension $23520-9937-13583=0$.

For a classically closed second-order source $\mathcal O_0$, the contribution on $\omega$ includes $Q_1\mathcal O_1$ with $Q_0\mathcal O_1=-Q_1\mathcal O_0$. The following calculation accounts for this correction by constructing protected representatives of the entire source. [First action](direct1.json.gz) · [Second action](direct2.json.gz) · [Unrestricted preimage coefficients](target_only.json.gz) · [A8 receipts](jobs.json)

## All incoming perturbative orders

The first column gives either one of the two equal spacetime Taylor degrees of the predecessor. The second gives its letter length. The two classical ranks come from finite-field images of integral maps, so they give lower bounds over the rationals. The cochain dimensions come from the colour characters. At derivative degree zero, A8 job 5681942 reproduced both ranks with independent rows and a different prime. [Archived source ranks](prior_source_ranks.json) · [Archived third-order source](prior_third_source.json) · [Independent third-order minors](third_source.json.gz)

| Derivative degree | Length | Cochains | Incoming classical rank | Outgoing classical rank | Rational cohomology upper bound | Independent graviton products |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 11 | 7,569 | 3,227 | 4,342 | 0 | 0 |
| 1 | 10 | 33,630 | 14,108 | 19,068 | 454 | 454 |
| 2 | 9 | 23,520 | 9,937 | 13,583 | 0 | 0 |
| 3 | 8 | 2,559 | 891 | 1,607 | 61 | 61 |

At length ten the 454 highest-weight combinations use 17,161 products of 612 single-trace factors. Each factor is a Taylor derivative of a graded-symmetrized trace of first superfield derivatives. Rebuilding every trace gives $Q_0=0$ on all 612 factors and verifies all 1,836 simple-root actions. Restricting to the three Cartan planes kills every classical boundary. A8 job 5681944 reconstructed the product coefficients by partitioning their factors among the planes and evaluated the determinant over the integers. Its residues are 901968 modulo 1000003 and 76599 modulo 1000033. At length eight the same computation gives residues 121552 and 859753. Thus all four rational upper bounds are attained. [454 product recipes](products.json.gz) · [Integer determinant](product_replay.json) · [61 product recipes](previous_products.json.gz) · [Second integer determinant](previous_product_replay.json)

The perturbative cohomological ring supplies a closed completion for each of these products. Every incoming second-order source therefore has zero obstruction. The first- and third-order sources have dimension zero. Higher incoming orders require negative derivative degree. Every positive outgoing order on $\omega$ requires negative flavour degree. Hence $\omega$ survives the full perturbative supercharge in this complex. [Perturbative ring theorem](https://github.com/lbotao399-stack/All-loop-holomorphic-kernels/blob/dcab65c7d758aef018b61e0b71ce398dca20c4d4/main.tex#L4695)

The four rows also exhaust the physical predecessor with the same dotted spin and flavour highest weight when instanton effects can mix letter lengths. If the protected single-trace factors have an exact cohomological product with these classical leading terms, the 515 incoming directions all admit exact closed representatives. This excludes an incoming instanton differential under that hypothesis. An outgoing instanton differential on $\omega$ remains open. [Exact primary product](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

To replay the classical result from the repository root with Python, NumPy, SciPy, python-flint and a C++ compiler, run

```sh
d=algorithms/dual_singlets_20261118
python3 "$d/sp3_omega_boundaries.py" --previous "$d/sp3_omega/previous.json.gz" --target "$d/sp3_omega/target.json.gz" --out /tmp/sp3-omega-boundary
python3 "$d/sp3_omega_replay.py" --page "$d/sp3_omega/target.json.gz" --cycle "$d/sp3_omega/cycle.json.gz" --lower "$d/sp3_omega/lowered.json.gz" --out /tmp/sp3-omega-replay
python3 "$d/sp3_omega_product_replay.py" --products "$d/sp3_omega/products.json.gz" --out /tmp/sp3-products.json
python3 "$d/sp3_omega_product_replay.py" --products "$d/sp3_omega/previous_products.json.gz" --out /tmp/sp3-previous-products.json
```

Use `sp3_omega_actions.py --task first` or `--task second` with the same target page and an output directory to regenerate the projected actions. [Paired quantum equations](../sp3_omega_coupled.py)
