# The classical split for the sixth glueball of F₄

In the relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ Yang–Mills theory, the classical supercharge $Q_0$ has rational rank 703 on the thirteen-letter source in Taylor diagram $(4,4,3,3)$. This is the candidate source for the second differential into the sixth power of the glueball. At lengths twelve, thirteen and fourteen the three cochain dimensions are

$$
(10007,3720,704),\qquad
\dim\ker Q_0|_{13}=3017,\qquad
\dim\operatorname{coker}Q_0|_{13}=1 .
$$

The kernel here still contains the length-twelve classical boundaries. The candidate quantum obstruction is a functional on this kernel. [Source selection, first correction and the criterion for a lift](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md#one-row-decides-the-second-differential)

## The upper bound

Polarize in four ordinary scalar variables $X_1,X_2,X_3,X_4$. The transposed highest weight is $(4,4,4,2)$. On the first $n$ variables use the Gram matrix with entries $\operatorname{Tr}_{26}(X_iX_j)$, where the trace is in the 26-dimensional representation. The product of its determinants for $n=3,4$ has the specified highest weight. On the four native Cartan generators the determinants are $1728$ and $5184$, so their product is $8957952$. Every classical image vanishes on commuting scalars because $Q_0$ splits the single second jet into a commutator. The terminal multiplicity 704 therefore gives the upper bound 703. The verifier computes the three cochain dimensions by the Weyl and symmetric-group character formulas. [Published character calculation](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/source_symmetry.json)

## The lower bound

The rows are scalar polynomials evaluated on the same 703 tuples of integral Lie algebra elements. Each of the four index columns, of heights $4,4,3,3$, is alternated. The first 542 rows are the published defining- and adjoint-trace rows. Those containing an explicit scalar commutator have the corresponding second-jet predecessor. For an ordinary product of at most three traces, two entries of the first four-slot alternation occur in the same trace. When the scalars commute that trace is symmetric in those entries. The whole alternation vanishes in the free commuting trace presentation. Reversing the adjacent swaps expresses it as a sum with an explicit commutator in each term. [Original trace minor and its construction](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page/trace_obstruction.md)

The next 158 rows attach the four-scalar alternations to two four-slot invariant tensors. Each endpoint contracts two Albert cubics through the 26-dimensional summand of the endomorphisms of the defining representation. At each end one designated pair of scalars enters as a commutator. The three-scalar alternation on that side acts on the three specified tensor slots. The final three rows insert another pair of Albert cubics between the two halves. Replacing the designated commutator at either endpoint with the corresponding second jet gives an explicit thirteen-letter predecessor. The native normalization of the cubic contractions agrees with the preceding fifth-power calculation. [Integer tensor construction and normalization](https://github.com/lbotao399-stack/quantum-lifts-results/tree/2a19364c8b5c32da2f713268b19a3dcc2580c6ed/reviews/independent-f4/fifth-glueball-partners)

| Prime | Trace rows | Four-cubic rows | Six-cubic rows | Determinant |
|---:|---:|---:|---:|---:|
| 1009 | 542 | 158 | 3 | 1005 |
| 1013 | 542 | 158 | 3 | 734 |

The rows, points and ordering are identical at both primes. Each row specializes a rational invariant tensor with regular reduction at these primes. The first nonzero determinant already proves the rational lower bound.

## Replay

Use Python with NumPy and SciPy. The first command checks the two matrices with independent modular elimination, recalculates the rational cochain dimensions and evaluates the Cartan determinants. The second also regenerates three complete columns at both primes. The third regenerates every entry. The two matrix files contain row-major unsigned 16-bit little-endian residues compressed with gzip.

~~~sh
python3 -B verify_classical.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B verify_classical.py --samples 0 351 702
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 -B verify_classical.py --full
~~~

The reference algorithms and the native representation data are from the pinned trace calculation. The selected six-cubic recipes numbered 0, 1 and 3 are recorded explicitly by the layout ordering in the rebuilder. The Albert contractions in the rebuilder use signed integer sparse-matrix actions. The recorded result includes both determinants and the Cartan values. [Reference data and code](https://github.com/lbotao399-stack/quantum-lifts-results/tree/e834b91eb860bd52c209201fd1e17a0ced5072dd/algorithms/dual_singlets_20261118/f4/j36_source/first_page)
