# The two proposed level-36 instanton partners in SU(2)

In undeformed four-dimensional $\mathcal N=4$ SYM, Kim, Kim, Lee and Park proposed a possible nonperturbative pairing between the two level-36 operators $O'_1$ and $O''_1$ in their equations (3.4) and (3.5). Both survive in physical supercharge cohomology under the physical inputs stated below. The calculation exhausts all singlet sources of a differential into either operator, allowing the bonus charge to change. [Original operators and first differential](https://arxiv.org/html/2606.27955v1#S3.SS2)

The physical inputs are the protected threshold $O_0$ with its weak-coupling leading class, the superconformal generators commuting with the chosen supercharge and inducing the stated classical actions, and the cohomological primary product of these operators and the stress-tensor multiplet with its classical leading term. The perturbative product and supercharge are those of Li, Liu and Chang. The passage to a physical instanton differential also assumes a filtration whose first page is the full perturbative cohomology at the given charge. [Threshold index argument](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex#L291-L295) · [Independent threshold calculation](https://github.com/lbotao399-stack/quantum-lifts-results/blob/03fc9a5d8a68a1b70072242a0770335df1d862d3/algorithms/dual_singlets_20261118/su2_threshold_index/README.md) · [Perturbative product](https://arxiv.org/abs/2611.24160) · [Physical product and currents](https://arxiv.org/html/2509.16737#S4.SS3.SSS1)

## Physical sources

The two targets are singlets of dotted spin and $SU(3)_R$. Their incoming source must also be a singlet. At the three common index charges $R_i+J_L=6$, let $t$ be the number of derivatives in each dotted direction and $Y$ the number of letters. The three flavour degrees are each $6-t$, so $0\le t\le6$. The charge equation is $2J_L=18-t-Y$. The sources for $O'_1$ and $O''_1$ have $2J_L=10$ and $9$, respectively. The entries below are the tree cochain dimension minus the incoming and outgoing rank lower bounds. The two computations use primes $1000003$ and $1000033$.

| $t$ | $Y=8-t$, incoming to $O'_1$ | $Y=9-t$, incoming to $O''_1$ |
|---:|---:|---:|
| 0 | $6-0-6=0$ | $25-6-17=2$ |
| 1 | $163-22-141=0$ | $584-141-442=1$ |
| 2 | $902-74-828=0$ | $3604-828-2775=1$ |
| 3 | $1267-118-1149=0$ | $6505-1149-5298=58$ |
| 4 | $560-21-536=3$ | $3074-536-2538=0$ |
| 5 | $30-2-28=0$ | $386-28-352=6$ |
| 6 | $0$ | $6-0-6=0$ |

The three classes in the first column and the $1,58,6$ even-length classes in the second column are products of components of the bilinear stress-tensor multiplet. The verifier expands each bilinear in nonabelian trace words, checks classical closure, and tests the two flavour raising operations and the spin raising operation on each proposed product combination. On the defining Cartan $\mathrm{Tr}(H^2)=2$. Tree boundaries vanish there. The integer Cartan minors have orders $3,1,58,6$ and are all nonzero. Their values and coordinates are supplied. These exact lower bounds attain the cohomology upper bounds in the table. Each product admits a closed completion. [Products and determinants](verified.json) · [Public graded trace construction](https://github.com/lbotao399-stack/quantum-lifts-results/blob/03fc9a5d8a68a1b70072242a0770335df1d862d3/algorithms/dual_singlets_20261118/sp3_omega_product_replay.py)

## The seven-letter source invisible on the Cartan

The threshold has seven letters, flavour degrees $(4,4,4)$ and zero dotted degree. Its flavour and spin off-diagonal generators vanish in classical cohomology by the singlet condition. All dotted lowerings annihilate its derivative-free representative. The possible odd conformal lowering has target dimension $11-6-5=0$. Its five diagonal linear Taylor weights in dotted-first order are $(0,0,4,4,4)$. [Threshold ranks and lowerings](https://github.com/lbotao399-stack/quantum-lifts-results/blob/03fc9a5d8a68a1b70072242a0770335df1d862d3/algorithms/dual_singlets_20261118/su2_threshold_index/README.md)

Act first with two distinct odd flavour translations on $O_0$, then with two substitutions of the same flavour indices by dotted indices, then with two ordinary translations. Pair the first ordinary dotted index with the first substitution by the alternating spin tensor, and pair the second ordinary index with the second substitution in the same way. The operator order and the spin signs are given in `source()` in the supplied script. The result is a spin and flavour singlet at $(4,4,4;2,2)$, $Y=7$.

To return this expression to $O_0$, reverse each word, interchange the two directions of every linear substitution, and replace each translation by its conformal lowering. On a single divided jet of total degree $k$, the lowering removes the selected slot with coefficient $k-1$ and its exterior sign. The source convention gives an additional minus for an odd conformal lowering. There are two odd lowerings in every return word. The dotted translation carries its incremented multiplicity. Applying the returned sum in the superconformal enveloping algebra gives

$$
24883200[O_0]\ne0 .
$$

The calculation uses 24 nonzero ordered source words. Their individual return scores on the spin contraction are $-691200$ and $1382400$. The verifier checks 77,175 signed generator-bracket identities on divided jets. Commuting the lowered generators through the word then uses only the five diagonal weights and the annihilators above. This proves nonvanishing in classical cohomology. The descendant attains the upper bound one, and has a physical closed lift from the protected threshold. [Signed operations and return](descendant.py)

## The differential

At $t=0$, $Y=9$, KKLP's two independent classes are the heavy core $O_1$ and the dressed threshold $O'_1$. They found $d_1[O_1]\ne0$ and $d_1[O'_1]=0$. Their equation (3.4) gives $O'_1$ as a product of a bilinear graviton with a twice-translated $O_0$, contracted to a singlet. Equation (3.5) gives $O''_1$ as a differentiated bilinear graviton times $O_0$. Both have closed completions under the stated product and symmetry inputs. KKLP verified that $O''_1$ is nonzero modulo the complete incoming first image, even after enlarging the allowed leading source to all singlets at $t=0$. The resulting 930-dimensional span excludes $O''_1$. Higher perturbative orders entering its $t=1$ grade would require negative dotted degrees. The $t=0$ grade of $O'_1$ excludes every incoming positive perturbative order.

In the $2J_L=10$ source, three closed product classes exhaust the possible incoming domain. In the $2J_L=9$ source the one known first image leaves

$$
1+1+1+58+6=67
$$

independent closed leading classes before any incoming quotient. These span every surviving graded piece by the table. Subtracting their closed lifts successively in this finite charge filtration gives a spanning set with full closed representatives for the perturbative domain. Any incoming quotient preserves that property. A first nonzero instanton differential into $O'_1$ or $O''_1$ would act on one of these domains and therefore has rank zero. This is a protection argument with explicit physical hypotheses. The microscopic instanton integral is unevaluated.

## Replay

Use Python 3.13, NumPy, SciPy, SymPy, python-flint and a C++17 compiler. Set `REPO` to a checkout at `03fc9a5d8a68a1b70072242a0770335df1d862d3`. The full command rebuilds the native engine and recomputes all 28 finite-field blocks. The integer product determinants and the return calculation are recomputed in both modes.

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify.py --repo "$REPO" --rebuild --workers 2 | diff -u verified.json -
~~~
