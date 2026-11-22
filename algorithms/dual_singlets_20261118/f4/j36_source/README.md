# A classical preimage of the supplied $F_4$ source at level 36

Dot's `F4-j36-d2-verification.zip` supplies the 13-letter source $O_0=P^6(Na)$ in the full relative holomorphic jet complex. In its native divided jets and metric $B=\operatorname{Tr}_{26}/6$, the 12-letter invariant $X$ in this directory satisfies

$$
Q_0 X=O_0,\qquad d_2[O_0]=0.
$$

This determines the differential on the supplied source. Whether a different source hits the sixth glueball power $[S^6]$ remains open. [Native preimage](tree_preimage.json.gz) · [Full coefficient replay](verify.py)

## Trace identity and primitive

The archive contains 6,861 rational terms, with two-, three- and five-letter colour traces. On its four Cartan generators all 35 quartic coefficients of $12\operatorname{Tr}_{26}(X^4)-(\operatorname{Tr}_{26}X^2)^2$ vanish. Here $X$ in the trace identity ranges over colour matrices. Density of regular semisimple elements extends the identity to $\mathfrak f_4$.

For five independently labelled matrices, polarize the quartic polynomial on three matrices and the commutator of the other two. Use the coefficient of the product of the four polarizing parameters. The ten possible commutator pairs give ten relations on 22 trace products modulo cyclicity and orthogonal reversal. Their rank is six. The difference between the five-letter replacement in the replay and the original trace equals $-1/576$ times the six relations whose commutator uses only the first four labels. Multilinearity supplies the graded identities by restoring the sign of each permutation of odd coefficients. Two- and three-letter traces are retained. The native metric is one sixth of a two-letter trace. [Quartic and five-letter verification](verify.py)

This gives 10,401 ordered trace words. A tree split preserves the two dotted degrees and three flavour degrees in each trace. Grouping by the unordered list of these five degrees for the trace factors gives 870 blocks. Exact elimination in each block constructs the 6,417 terms of the stored preimage. Its only trace lengths are two and four. Direct differentiation followed by subtraction leaves zero on the full rational polynomial. Each summand is a product of colour traces and hence gauge invariant. [Block solver](construct.py) · [Exact counts](source-receipt.json)

The archive uses transposed dual generators. Transposition reverses the sign of the matrix commutator. The native replay checks all 1,378 unordered colour pairs, including equal pairs, with both output colour slots raised by the inverse metric. A separate component differentiation acts on the original 52-colour variables before multiplying the matrices. In the archive's exterior order it gives $2232,-24,5382$ on its three supplied probes for the matter fermion, differentiated scalar and second derivative of the ghost. [Native colour replay](native.py)

## The first correction and the numerical pairing

Choose $O_1=Q_1X$. The largest number of flavour slots on three letters is five for $X$ and four for $O_0$. The bare $Q_2$ requires six, so $Q_2X=Q_2O_0=0$. With the archive's supercharge identities,

$$
Q_0O_1=-Q_1O_0,\qquad
Q_1O_1=Q_1^2X=-Q_0Q_2X-Q_2Q_0X=0 .
$$

For any other first correction satisfying the same descent equation, $O_1-Q_1X$ is tree closed. Its first image accounts for the entire possible change in the second residual. The residual therefore represents zero on the second page.

A8 job `5693514` applies the native two-input $Q_1$ directly to the 684 matter-fermion terms of the collected $O_0$. It contracts both adjoint brackets, the inverse metric and the exterior spectators. The other two source profiles have fewer than three flavour slots on any two inputs. The two archived component sets contain 490 and 2,491 distinct entries, with 2,865 nonzero entries in their union. All 2,981 match the archived anchor divided by six with a single overall minus sign in the stated plus-first, positive local-kernel orientation. The archived scalar calculation gives $142560/49$. Its nonzero value on this tree boundary prevents its use as a functional on the second page. The test below isolates a discrepancy in the archived endpoint functional on a tree image. [Direct contraction](first_action.py) · [A8 result](a8-first-action.json)

## A boundary that tests the endpoint functional directly

Use the archive's $S$, its divided first ghost jets $\lambda_{\dot\alpha}$ and its first differentiated scalar $D_{\dot\alpha}\phi^I$. The dotted alternating tensor has $\epsilon^{+-}=1$, and the flavour tensor has $\epsilon_{123}=1$. In this basis take the gauge-invariant cochain

$$
Y=S^3\epsilon_{IJK}\epsilon^{\dot\alpha\dot\beta}\epsilon^{\dot\gamma\dot\delta}
\operatorname{Tr}_{26}\!\left((D_{\dot\alpha}\phi^I)\phi^J\lambda_{\dot\beta}
\phi^K\lambda_{\dot\gamma}\lambda_{\dot\delta}\right).
$$

Each of its twelve letters carries at most one flavour slot. The two inputs of $Q_1$ require three in total, so $Q_1Y=0$ as a cochain. Write $\lambda_6$ for the archive's classical target functional and $\eta_0,\eta_K$ for its two numerical endpoint tables. The asserted first adjoint identity $\lambda_6 Q_1=(\eta_0+\eta_K)Q_0$ would give zero on $Y$.

Only the differentiated scalar contributes to $Q_0Y$. With the archived transposed generators its ordered tree split inserts $\phi^I\lambda_{\dot\alpha}-\lambda_{\dot\alpha}\phi^I$. Applying the two endpoint tables to the plus-first coefficients of this tree image gives

$$
\eta_0(Q_0Y)=-\frac{98441322408000}{7}.
$$

| Prime | $\eta_0(Q_0Y)$ | $\eta_K(Q_0Y)$ | $(\eta_0+\eta_K)(Q_0Y)$ |
|---:|---:|---:|---:|
| $1000003$ | $559299$ | $592975$ | $152271$ |
| $1000033$ | $292644$ | $926374$ | $218985$ |

The first replay evaluates the commutator and the three glueball factors directly. The second differentiates the trace by the native jet rules, producing 36 canonical trace words and 144 ordered colour assignments. The two evaluations agree on all 2,981 archived component requests. The 7,276 unfolded queries for $\eta_0$ reproduce the first residue. The 40-coordinate and 477-query evaluations of $\eta_K$ agree at both primes. Substituting the archive's own input reproduces its endpoint scores. [Exact values and component digest](ward-boundary.json) · [Endpoint evaluation](ward_boundary.py) · [Independent trace calculation](ward_boundary_words.py)

Both nonzero residues disprove the indicated adjoint identity for the endpoint tables as applied to physical cochains. This test uses the local flavour rule and the tree action. A different $F_4$ source, or a proof that $[S^6]$ survives the incoming $d_1$ using a corrected dual, requires a separate calculation.

## The two missing metric factors

The classical target routine in the archive evaluates each harmonic branch at its highest exterior vector $H$ and opposite vector $L$ with the induced exterior metric. In the root order used by that routine the six metric entries multiply to $-4$ for the branch of dimension $7113106$ and to $2$ for the branch of dimension $1850212$. The exporter `export_eta0.py` multiplies by the branch dimension, harmonic coefficient and ordering sign. It omits this exterior metric in both branches.

For example, the archive's classical mixed-boundary row has entries $-836/1365$ and $736/1155$. Its pairing with the classical target coefficients $(2392,2299)$ vanishes. Replacing the two extremal values by the values used in the exporter changes those coefficients to $(-598,2299/2)$, and the same row gives $38456/35$.

We reconstructed both rational contributions of the scalar lift from the $52$-colour brackets and $B=\operatorname{Tr}_{26}/6$. With the factors omitted, the two contributions reproduce all $7276$ unfolded weights at the first prime and all $2491$ folded weights over the rationals. Four additional raw components enter when the metric factors are restored. The native ordered one-loop action and the archive's literal trace contraction give the same four values in the positive raw-anchor convention,

$$
(-51840,\ 51840,\ 38880,\ -38880).
$$

In the sign convention of the archive's final score, all entries in the following row are to be divided by $7$.

| First branch before its metric factor | Second branch before its metric factor | Corrected $\eta_0$ | Archived $\eta_K$ | Sum |
|---:|---:|---:|---:|---:|
| $-56702201707008000$ | $28351100853504000$ | $283511008535040000$ | $-283511008535040000$ | $0$ |

Thus the corrected pairing on the supplied $O_0$ vanishes over $\mathbb Q$, in agreement with $O_0=Q_0X$. On the independent test $Q_0Y$, the corrected $\eta_0$ is $984413224080000/7$. Adding the archived $\eta_K$ gives zero at both verification primes. The replay records the opposite overall sign for the source row because it reports the raw anchor divided by six, before the minus sign in the archive’s scoring formula. [Metric, raw components and exact arithmetic](endpoint-metric.json) · [Reconstruct both branch weights](endpoint_metric.py)

## Replay

Supply the original zip with SHA-256 `0a06b306a12abbe77cadbcfab616c08c579d82086e9798aece27b8c83a0f7b52`. The first two commands use the Python standard library. The block solver also uses the repository trace runtime with NumPy, SciPy and python-flint.

~~~sh
python3 verify.py /path/to/F4-j36-d2-verification.zip
python3 ward_boundary.py /path/to/F4-j36-d2-verification.zip
python3 ward_boundary_words.py /path/to/F4-j36-d2-verification.zip
python3 endpoint_metric.py /path/to/F4-j36-d2-verification.zip
python3 first_action.py /path/to/F4-j36-d2-verification.zip --workers 4
python3 construct.py /path/to/F4-j36-d2-verification.zip
~~~
