# Local contacts on the SU(3) level-30 threshold

The level-30 $SU(3)$ source gives a sharp test of proposed instanton contacts. At the conformal origin, three active input letters give zero on the component of Dot's Cartan target with five undifferentiated scalars, two matter fermions and two pure dotted gaugino jets. A four-input ansatz fits all 120 coordinates in the combined support. Four classical boundaries then exclude all 86 pure-fermion colour patterns together with the one mixed tensor used in the fit. The exact rank rises from seven to eight upon adjoining the required values. The renormalized one-instanton Ward coefficient remains open. [Dot's 6,139-term source](https://github.com/lbotao399-stack/quantum-lifts-results/tree/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-exact-input) · [Dot's 5,080-term target](https://github.com/lbotao399-stack/quantum-lifts-results/tree/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-t2-primary) · [Independent radial boundary projection](https://github.com/lbotao399-stack/quantum-lifts-results/blob/d28dfd7147da392772dd6dd292b921b479532fa7/reviews/independent-g2-and-r007/r007/README.md)

## Number of active inputs

An active input is a letter differentiated in the local contact. Spectators keep their fields and their dotted Taylor indices. The source has ten letters and odd Taylor degrees $(7,4,4)$. The target has nine letters and degrees $(5,2,2)$, with two Taylor indices of each dotted sign. Removing six odd Taylor indices takes at least three source letters. With three inputs there must be one matter fermion of each flavour. Their replacement has two pure dotted jets. Dotted rotational invariance pairs two copies of the same spin, so both jets have dotted degree two. The numbers of target trace words at the respective minimum input numbers $3,4,5,6,7$, computed by matching spectator species to every source content, are

$$
(69,\ 2007,\ 2877,\ 126,\ 1).
$$

This count permits arbitrary reconnections of the colour traces. The 69 words consist of 45 with the two opposite extremal jets and 24 with the two middle jets. They give 240 monomials on the rank-two Cartan subalgebra before selection of a polarization. [Endpoint reconstruction](cartan_contacts.py) · [Counts](cartan_contacts.json)

## The cubic colour test

Take the three scalar matrices diagonal and keep all eight independent colour entries of each odd matter fermion. With the three diagonals written as first entry, second entry and minus their sum, ordinary polynomial multiplication and exterior multiplication give the following numbers of nonzero source monomials by off-diagonal fermion degree $0,\ldots,5$.

$$
(0,\ 0,\ 3986,\ 0,\ 0,\ 0).
$$

Torus invariance pairs the two off-diagonal entries on opposite roots. As an independent colour check, at the integer specialization in the verifier multiply every off-diagonal fermion entry by $t$. The coefficient of its five exterior generators is $-77559704t^2/3$. Matrix multiplication at $t=0,1,-2$ gives $0,-77559704/3,-310238816/3$. [Exact expansion and separate matrix evaluation](cartan_contacts.py)

Use an ordinary test matrix $X$ to polarize the three fermion derivative colours. Contraction with the flavour alternating tensor and permutation of the three odd derivatives make their colour slots symmetric. Let $U,V$ be the two output jet colours. The spin singlet of two spin-one jets is symmetric in spin and the jets are odd, so the colour tensor is antisymmetric in $U,V$. Both endpoint classes are odd under colour conjugation, hence their coupling preserves colour conjugation. The five-index colour invariants of traceless matrices are spanned by traces of length five and products of lengths two and three. Paired colour epsilon tensors reduce to Kronecker deltas. Reversal with the five transpose signs removes the odd colour tensor. The traceless three-colour Cayley–Hamilton identity reduces the term with $X^3[U,V]$. The two remaining tensors are

$$
\operatorname{Tr}(X^2)\operatorname{Tr}(X[U,V]),\qquad
\operatorname{Tr}(X^2UXV-X^2VXU).
$$

Set $U=\operatorname{diag}(1,0,-1)$ and $V=\operatorname{diag}(0,1,-1)$. These span the Cartan. The first tensor vanishes and the second gives $3(X_{12}X_{23}X_{31}-X_{13}X_{32}X_{21})$ in ordinary matrix-entry notation. Both terms require three off-diagonal fermions in the source. The polynomial identity above kills them before any dynamical coefficient enters. [Cubic trace expansion](cartan_contacts.py)

## A functional on the four-input contact

In the same polarization, set the middle dotted jet to zero and attach independent exterior variables to $U,V$. For the four-fermion input, polarize three derivative colours on $X$ and label the fourth $Y$. Label the output matter-fermion colour $Z$. Flavour invariance pairs its flavour with the fourth derivative and antisymmetrizes the other three. Permutations of which derivative is distinguished are included by polarizing the colour tensors. Products of traces span the invariants of the seven colour slots. Reversal of all colour traces and exchange of $U,V$ each give a minus sign. The 337 cyclic trace products yield 86 projected trace patterns. There are 288 tests with two opposite root inputs, two Cartan inputs and a Cartan output $Z$. These exhaust the support by the quadratic source identity. Exactly three colour rows are independent. [Colour enumeration and rational elimination](four_fermion_contacts.py)

Use the three independent rows in the order of their unprojected trace words

$$
(\operatorname{Tr}(UX)\operatorname{Tr}(VXXYZ),\
 \operatorname{Tr}(UX)\operatorname{Tr}(VXXZY),\
 \operatorname{Tr}(UXVXXYZ)).
$$

The verifier sums the three-input polarizations and both signed projections with the defining trace. The projected target has 80 nonzero coordinates. On the four coordinates specified by scalar powers and fermion slots in the receipt, the three colour rows and the target give the columns of the matrix

$$
\begin{pmatrix}
-2592&0&3456&-1344/5\\
2400&912&-2016&1344/5\\
-1872&-384&1664&-1344/5\\
-528&1056&1408&1152/5
\end{pmatrix},
\qquad \det=311115644928.
$$

The row $(-11,-36,-36,18)/7104$ sends this matrix to $(0,0,0,1)$. Every classical boundary restricts to zero on the Cartan, where the tree supercharge is zero. This functional detects the target independently of classical boundary subtractions. Its value on other physical primaries requires separate projection. [Four coordinates and conventions](four_fermion_contacts.json)

For the mixed contact keep the three fermion derivative colours on $X$ and use $Y$ for the scalar derivative colour. Pair its flavour with that of the output scalar $Z$. Take the explicit tensor

$$
\operatorname{Tr}\!\left(U[X,Y]\right)\operatorname{Tr}(VX)\operatorname{Tr}(ZX)
-\operatorname{Tr}\!\left(V[X,Y]\right)\operatorname{Tr}(UX)\operatorname{Tr}(ZX).
$$

The off-diagonal scalar and fermion derivatives are on opposite roots and lie in one source trace. Summing their placements gives $(5184,-4128,2064,6528)$ on the four coordinates, with functional value $702/37$. At five scalars of the first flavour with diagonal $(2,-3,1)$ a direct unsplit substitution in all 49 contributing words gives $217728$. With the mixed tensor as fourth column the determinant is $5902788722688$. Solving on the four displayed coordinates gives

$$
\left(\frac{62453}{347490},\
\frac{8071}{694980},\
-\frac{10213}{463320},\
\frac{37}{702}\right).
$$

Applying these weights to the four contact polynomials agrees with the specified component of Dot's target on all 120 coordinates in their combined support. This establishes that four active inputs suffice for the projected polynomial matching problem. A test of a physical instanton requires the coefficients of the matched supercharge on perturbative cohomology. [Exact fit and independent substitution](mixed_contact.py) · [Rational values](mixed_contact.json)

## Classical boundary test

For a correction that induces a map on classical cohomology, its value on a tree boundary must be a tree boundary. Every such output has zero Cartan restriction. The first test word uses $f$, the undifferentiated triple odd Taylor jet. In the public BMN matter-fermion conventions its tree variation is the negative of the flavour contraction of $[\phi^i,\psi_i]$. The replay obtains all six ordered terms from Dot's native letter differential. With the standard alternating tensor on the three matter flavours, the first invariant predecessor is

$$
\epsilon^{ijk}\operatorname{Tr}\!\left((\phi^1)^3 f\,\psi_i\psi_j\psi_k\psi_l\phi^l\right).
$$

Tree variation on a matter fermion leaves $f$ as a spectator, which the specified output projection removes. Tree variation on $f$ gives the surviving words. On this boundary the four fitted contact tensors have functional values $(757/1184,11/32,-591/592,-201/1184)$. The fitted combination gives $6032431/45714240$. Thus agreement on all target coordinates allowed a false positive.

For all four predecessors, the receipt records cyclic words with $0$ denoting $\phi^1$, $a,b,c$ the three epsilon-contracted fermions, and $d,e$ the flavour-contracted matter fermion and scalar. A slash starts a new trace. Impose the four displayed target coordinates on R-007, and zero on the functional of each exact boundary. The resulting eight-by-87 rational matrix has rank seven. Its integer left-null relation annihilates every column and takes the value $22052448777162432/5$ on the required column. Consequently the tested family cannot give any nonzero multiple of the specified target projection while also descending on these four boundaries. The boundary test uses all admissible root configurations of each predecessor, including three and four off-diagonal derivatives. Further mixed tensors, terms acting on $f$, and corrections from matching the perturbative representatives lie outside this family. [All matrix entries, words and integer witness](contact_ward.json) · [Native tree kernel and transposed colour replay](contact_ward.py)

Run from this directory in a checkout with the pinned Dot input commit.

~~~sh
python3 cartan_contacts.py | diff -u cartan_contacts.json -
python3 four_fermion_contacts.py | diff -u four_fermion_contacts.json -
python3 mixed_contact.py | diff -u mixed_contact.json -
python3 contact_ward.py | diff -u contact_ward.json -
~~~
