# First supercharge on the two SU(2) heavy-core singlets at levels 60 and 72

In the relative Taylor-jet complex of undeformed four-dimensional Yang–Mills theory with $\mathcal N=4$, let $O_n$ be the heavy core of Kim, Kim, Lee and Park, equation (2.15). The graviton $(w_2)^k$ is their equation (2.14), and $D^i$ is odd flavour translation. At $n=3,4$ the two columns of the first-supercharge matrix, in order, are the classes of

$$
O_n,\qquad \epsilon_{ijk}(w_2)^k D^i D^j O_{n-1}.
$$

Here $i,j,k$ are the three flavour components and $\epsilon_{ijk}$ is their alternating form. The level is $j=24+12n$ and the letter length is $7+2n$. The defining two-dimensional trace and the ordered-pair normalization are inherited from the level-48 calculation. Its colour coordinates have $f_{\rm KKLP}=-f/2$ and $Q_{\rm KKLP}=-2iQ_0$. [KKLP heavy cores and previous individual first lifts](https://arxiv.org/html/2606.27955v1#S3.SS3) · [Ordered first-supercharge kernel](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/code/fortuity/full/q1.py) · [Level-48 rank-two control and normalization](https://github.com/lbotao399-stack/quantum-lifts-results/tree/3843a42d05f064b9fb855571c03917707dee9a20/reviews/independent-su2/levels-48-50)

At both levels the source has $28-6-20=2$ tree classes. The two rows at each level are integer covectors on the one-plus, one-minus dotted target. Applying them to every incoming tree boundary gives zero. Applying them to the two $Q_1$ images gives

$$
j=60:\quad
\begin{pmatrix}79920933120&1278734929920\\30964076640&503391859200\end{pmatrix},
\qquad \det=636700739987747635200,
$$

$$
j=72:\quad
\begin{pmatrix}1321463185920&19821947788800\\-829440&-14976000\end{pmatrix},
\qquad \det=-3349116298395648000.
$$

The first-source kernel has dimension zero in these two singlet columns. The result adds exact tests of both directions at these levels. The operator family and earlier first lifts are due to KKLP. Page two has no source left in these columns. These ranks count supercharge differentials.

The verifier uses the invertible frame of the three scalars. At fixed fermion and field-strength content the scalar degree is fixed, so homogeneity and flavour covariance determine each coefficient on the dense invertible locus. The remaining colour action is the adjoint action of $SL(2)$. Restricting the field strength to the diagonal generator $H$ and using the raising and lowering generators recovers its transverse derivatives. On the zero-weight space the operator formed by lowering after raising has eigenvalue $\ell(\ell+1)$ on the integer-spin $\ell$ representation. A homogeneous polynomial of degree $m$ in the adjoint field contains spins $m,m-2,\ldots$. The projector in the verifier retains precisely these spins. Homogeneity fixes the radial derivative. Spin-singlet projection discards an isolated two-dotted jet, which is symmetric in its spin indices. A separate Weyl-character computation counts all predecessor cochains. Their numbers are 1,195 and 1,210. Colour, flavour and spin symmetry are reductive in characteristic zero, so invariant predecessors span the invariant boundaries. Each row of the displayed matrices is evaluated on the unreduced two-polynomial image in the same frame.

The marked quantum product induces the derivation rule on each page. A product whose factor is removed by $d_1$ requires its own cohomology calculation at the next page. [Quantum supercharge and marked product](https://arxiv.org/abs/2611.24160)

## Replay

Python 3.13 and NumPy suffice. The command rebuilds the colour and flavour singlets, checks tree nilpotence on every boundary column, checks all dual pairings, and compares the three matrices with the receipt. The level-48 calculation is a control.

~~~sh
shasum -a 256 -c SHA256SUMS
python3 -B verify_levels.py --all-columns 2 3 4 | diff -u verified.json -
~~~
