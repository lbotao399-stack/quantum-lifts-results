# The eighth glueball power in SO(13)

In the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, the first incoming supercharge differential that kills $S^8$ at $j=48$ has order two or three. Here $S$ is Dot's defining-trace glueball, with the dotted order used in the $SO(13)$ certificate. The result uses Dot's order-three source for $S^7$. [Signed source and normalization](https://github.com/lbotao399-stack/quantum-lifts-results/blob/e1d184d388c588c0a5a825a19ecbca4e24427242/reviews/dot-orthogonal-d3/REPLAY.md)

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

If $[S^8]$ survives the second differential, the displayed cochains give $d_3[Sx_0]=-921600[S^8]$. Otherwise the target is an incoming $d_2$ boundary. The coefficient alone permits both possibilities. The second-dual weight has 301 balanced flavour-singlet finite-colour relations, of which 81 are dotted-spin singlets. These are dimensions of the relation spaces.

Run the following from this directory with the repository's trace runtime and the directory containing Dot's three encoded archives. The script checks the first-page statement, both versions of Dot's third coefficient, all spectator pairs on the five first and twelve second jets, and the full first action on six nonzero trace inputs.

~~~sh
python3 spectator.py /path/to/dot-orthogonal-d3
~~~
