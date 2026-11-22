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

A8 job `5693514` applies the native two-input $Q_1$ directly to the 684 matter-fermion terms of the collected $O_0$. It contracts both adjoint brackets, the inverse metric and the exterior spectators. The other two source profiles have fewer than three flavour slots on any two inputs. The two archived component sets contain 490 and 2,491 distinct entries, with 2,865 nonzero entries in their union. All 2,981 match the archived anchor divided by six with a single overall minus sign in the stated plus-first, positive local-kernel orientation. The archived scalar calculation gives $142560/49$. Its nonzero value on this tree boundary prevents its use as a functional on the second page. The precise failing step in the archived dual construction or its transport is unresolved. [Direct contraction](first_action.py) · [A8 result](a8-first-action.json)

## Replay

Supply the original zip with SHA-256 `0a06b306a12abbe77cadbcfab616c08c579d82086e9798aece27b8c83a0f7b52`. The first two commands use the Python standard library. The block solver also uses the repository trace runtime with NumPy, SciPy and python-flint.

~~~sh
python3 verify.py /path/to/F4-j36-d2-verification.zip
python3 first_action.py /path/to/F4-j36-d2-verification.zip --workers 4
python3 construct.py /path/to/F4-j36-d2-verification.zip
~~~
