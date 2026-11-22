# Two independent first images at SU(2) levels 48 and 50

In the derivative-free sector of undeformed four-dimensional $\mathcal N=4$ SYM, the outgoing first supercharge differential has rational rank two on the level-48 flavour singlet and on the level-50 flavour fundamental. These exhaust both source spaces. The level-50 flavour representation with Dynkin labels $(0,2)$ has zero classical cohomology at every length. The statements concern these three highest-weight columns.

Kim, Kim, Lee and Park established the first lift of the heavy core $O_2$ at level 48. Their equation (2.15) supplies both heavy cores used below, and (2.14) supplies the graviton $(w_2)^k$. The first kernel is the ordered Budzik–Kulp kernel in the public Fortuity convention. The independent result here is the second direction, the rational rank-two certificate, and its propagation to level 50. [Heavy cores and the earlier lift](https://arxiv.org/html/2606.27955v1#S3.SS3) · [Public ordered-pair algorithm](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/code/fortuity/full/q1.py) · [All-order supercharge and marked product](https://arxiv.org/abs/2611.24160)

## Classical source

Let $D^i$ be odd translation in the $i$th flavour direction and $S_i$ its classical lowering in the divided-jet convention. For the ordered basis at level 48 take

$$
O_2,\qquad
\epsilon_{ijk}(w_2)^kD^i D^j O_1 .
$$

Write $B$ for the second polynomial in this calculation. It is the descendant of a product with the level-36 heavy core. The graviton obeys $D^i(w_2)^k=0$. Fix the defining two-dimensional trace. In the standard $E,H,F$ basis of $\mathfrak{sl}_2$, the tree routine uses $Q\psi_i=[\phi^{i+1},\phi^{i+2}]$ with cyclic flavour labels and $Qf=[\phi^i,\psi_i]$. The program evaluates the physical polynomial (2.15) by substituting $f_{\rm KKLP}=-f/2$, with unchanged $\phi,\psi$. The native three-flavour superfield coefficient is $-f$, and its two-flavour coefficients for the complementary pairs are $(\psi_1,-\psi_2,\psi_3)$. This gives exactly (2.3) with $Q_{\rm KKLP}=-2iQ$. The overall ordered-pair constant is one in the cited public algorithm. The two $D$ actions on $B$ have the same sign in the frame and native bases.

The three scalar vectors form a flavour frame on the dense set of invertible scalar matrices. At fixed undotted fermion and field-strength content the scalar degree is fixed. Flavour covariance and homogeneity therefore determine a polynomial from its value on the unit scalar frame, with the auxiliary test vector retained for a nonsinglet. The first fundamental theorem for $SL(3)$ gives pairings and determinants of three flavour vectors or covectors. Contracting two opposite determinants reduces them to pairings. The frame code evaluates these generators, takes independent rational rows, and applies the full tree differential there. The residual colour group is $SL(2)$. Its singlet dimension, on both cochains and images, is the weight-zero dimension minus the weight-one dimension in adjoint weight units. An independent Weyl character checks the cochain dimension at every length.

At length eleven the unreduced cochain dimension, incoming tree rank, outgoing tree rank and cohomology dimension are

| Level | Flavour degrees | Cochains | Incoming | Outgoing | Cohomology |
|---:|---|---:|---:|---:|---:|
| 48 | $(8,8,8)$ | 28 | 6 | 20 | 2 |
| 50 | $(9,8,8)$ | 42 | 13 | 27 | 2 |
| 50 | $(9,9,7)$ | 21 | 3 | 18 | 0 |

All other lengths in the three columns have zero cohomology. The row reductions on the complete incoming scalar-frame spaces give the following equalities in classical cohomology. The flavour contractions use the fundamental and its dual.

$$
S_i[102O_2-B]=0,\qquad
D^i S_i[O_2]=-\tfrac13[B],\qquad
D^i S_i[B]=-34[B].
$$

## First target

Set all undifferentiated bosons to zero in the spin-singlet first target of flavour $(7,7,7)$. The specialization commutes with the tree differential since the tree variation of each omitted boson vanishes there. A term containing one second dotted jet belongs to the spin triplet. Every surviving spin-singlet term contains two first dotted jets and all nine undifferentiated fermion components. Its two flavour subsets are disjoint and exhaustive. Ordering these pairs and the three colour coordinates gives twelve rows.

There are 21 possible predecessor monomials of the correct Cartan weights on this specialization, even allowing coefficients that do not extend to invariant polynomials. Their integer tree matrix has rank eleven. In the sorted native basis the covector $\lambda$ equals

$$
(1,\tfrac12,1,\;-1,-\tfrac12,-1,\;1,\tfrac12,1,\;1,\tfrac12,1).
$$

It annihilates all 21 columns. Apply the ordered kernel using the inverse of the defining trace form on $E,H,F$. Direct exterior multiplication gives

$$
\lambda Q_1O_2=-9600,\qquad \lambda Q_1B=-172800 .
$$

On derivative-free words, both $D^i$ and $S_i$ anticommute with $Q_1$. The replay checks every flavour and ordered pair. Their tree anticommutators on the dotted target vanish on invariant traces, including the inner gauge term for $D^i$. Thus $\lambda D^i S_i$ gives a second functional on tree cohomology. Their evaluations on the ordered sources form the matrix

$$
\begin{pmatrix}
-9600&-172800\\
57600&5875200
\end{pmatrix},
\qquad
\det=-46448640000 .
$$

At level 50 the two highest-weight vectors occur among $D^iO_2,D^iB$. On a singlet source the contracted anticommutator is $\{S_i,D^i\}=-48$. The eigenvalues of $S_i D^i$ on the two-dimensional level-48 space are $-48,-14$. If a linear combination of the two translated flavour multiplets had zero first image, applying $S_i$ and contracting would give zero first image on the corresponding combination after this invertible map. The nonzero determinant excludes it. Consequently $E_2$ is zero in the two specified derivative-free source columns at both levels.

## Replay

Python 3.13 and NumPy suffice. The verification rebuilds the characters, integer boundary spaces, the two polynomials, all 21 specialized predecessors, and both ordered-kernel scores.

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify.py | diff -u verified.json -
~~~
