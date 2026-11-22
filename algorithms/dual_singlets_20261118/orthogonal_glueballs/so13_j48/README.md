# The eighth glueball power in SO(13)

In the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, the first incoming supercharge differential that kills $S^8$ at $j=48$ has order three. Here $S$ is Dot's defining-trace glueball, with the dotted order used in the $SO(13)$ certificate. Its eighth power survives both earlier differentials. The source is the product of $S$ with Dot's order-three source for $S^7$. [Signed source and normalization](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md)

At the pure-gaugino weight $(8,8)$ the stable classical quotient has just the product of eight bilinears. Sorting each even trace kills a repeated odd letter, and an odd-length trace is classically exact by reversal. The full finite-colour relation counts are zero at $(8,8)$ and at the first dual weight $(7,7;1,1,1)$, where the three last entries count the scalar flavours. The functional extracting $S^8$ has a first correction supported on 13 ordered words. Its transposed first equation holds on all 33 potentially contributing columns. It descends to thirteen colours, so $[S^8]$ is nonzero on the second page. All positive-order charges annihilate the pure representative. [First dual](../../orthogonal_glueball_second_pair.py) · [Finite-colour characters](../../orthogonal_glueball_relation_span.py)

Write $x_0,x_1,x_2$ for Dot's source and its first two corrections, and $x_3$ for the last pure-gaugino primitive in that certificate. With its loop parameter $\hbar$ and full supercharge $Q$, the chain satisfies

$$
Q(x_0+\hbar x_1+\hbar^2x_2+\hbar^3x_3)=-921600\,\hbar^3 S^7 .
$$

For each of these cochains, multiplication by the even spectator $S$ commutes with $Q$. The classical charge is a derivation and kills both spectator letters. The first correction selects two letters with at least three flavour slots in total. Each spectator gaugino has zero, and every letter of $x_0,x_1,x_2,x_3$ has at most two. Thus every active pair is internal to the original cochain. For order $r\ge2$, the $r+1$ selected letters have at most $r+2$ flavour slots because each cochain has at most one second jet. The required number is $3r>r+2$. Hence

$$
Q\bigl(Sx_0+\hbar Sx_1+\hbar^2 Sx_2+\hbar^3 Sx_3\bigr)
=-921600\,\hbar^3 S^8 .
$$

## Survival through the second differential

The stable second dual has 251 ordered-word coefficients. It obeys the second transposed equation on all 1,266 potentially contributing columns. Inverse reconnection checks the first and second equations on 39 and 1,362 columns including extra zero columns. At the second-dual weight $(6,6;2,2,2)$ the orthogonal invariant-tensor character gives 301 finite-colour relations with balanced dotted weight and singlet flavour. This includes all 81 dotted-spin singlets. The flavour projection commutes with the supercharge, and all $SO(13)$ adjoint invariants are $O(13)$ invariant.

For any relation $R$, write $c_1(R),c_2(R)$ for its coefficients on entries 1 and 2 in the zero-based distinct-plane Cartan basis supplied by the verifier. Restriction to commuting Cartan matrices kills classical boundaries. Define $\zeta(R)$ to be the stable second dual evaluated on the ordered trace polynomial $R$. On every one of the 301 supplied polarized Pfaffians the following equality holds over the rationals,

$$
\zeta(R)+\frac{9}{1280}c_1(R)+\frac{1}{80}c_2(R)=0.
$$

For example, row 41 gives $18432-1474560/80=0$, and row 64 gives $13824-9(1966080)/1280=0$. The 301-by-301 ordered-word minor of these Pfaffians has determinant $427177$ modulo the prime $1000003$. Every row specializes to zero in thirteen colours, its polarized Pfaffian degree exceeding 13. The nonzero minor proves independence over the rationals, so the character count proves completeness. Adding the two Cartan functionals to the stable dual preserves its transposed equation and makes it vanish on the entire relation space. The corrected dual therefore exists in thirteen colours and annihilates every incoming second differential to this target. [Basis, probes, minor and rational scores](relations.json.gz) · [Rational and modular replay](relations.py)

The product chain consequently gives

$$
d_3[Sx_0]=-921600[S^8]\ne0 .
$$

Dot's uniform pure-glueball threshold predicts this product case. A separate calculation reduces the finite identities to eight path types and obtains the same two Cartan coefficients on 56 polarizations. The ordered-word minor above tests the whole 301-dimensional relation space before that reduction. [Path-reduction calculation](https://github.com/lbotao399-stack/quantum-lifts-results/tree/f05b289fc69fb5f65bd8286c09060156a554f23d/reviews/orthogonal-mixed-lifts)

Run the following from this directory with the repository's trace runtime and the directory containing Dot's three encoded archives. The script checks the first-page statement, both versions of Dot's third coefficient, all spectator pairs on the five first and twelve second jets, and the full first action on six nonzero trace inputs.

~~~sh
python3 spectator.py /path/to/dot-orthogonal-d3
python3 relations.py --rebuild --workers 4
~~~
