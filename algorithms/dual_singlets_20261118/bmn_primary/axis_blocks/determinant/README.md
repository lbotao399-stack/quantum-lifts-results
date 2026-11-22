# One possible second target in the fifteen-letter G₂ singlet sector

For the level-36 flavour singlet of $G_2$, a fifteen-letter source can have at most **one** classical target for the second supercharge obstruction under filtered Taylor equivariance. The target Taylor diagram is $(5,5,3,3)$ at length fourteen. The high pair of axes must carry $\det^5$ and the transverse pair $\det^3$. The latter contains one copy each of the colour singlet and $(4,2),(6,0),(8,0)$ in short-root-first Dynkin labels. The high pair has zero singlets. At this total length the earlier pure-axis filtration excludes every occupation with a jet connecting the two pairs. Keeping both axes of each pair together gives a filtration invariant under both changes of Taylor coordinates. The selected restriction of $(5,5,3,3)$ occurs once. [Occupations and transverse cycles](../README.md)

| High-pair colour | Chain multiplicities at lengths 7, 8, 9 | Rank from 8 to 7, lower bound | Rank from 9 to 8, lower bound | Rational homology at length 8 |
|---|---|---:|---:|---:|
| $(4,2)$ | $396,411,194$ | 253 | 158 | 0 |
| $(6,0)$ | $649,543,231$ | 355 | 188 | 0 |
| $(8,0)$ | $91,103,54$ | 60 | 42 | 1 |

The last entry is exact in the **high-pair** complex. Tensoring it with the transverse copy leaves one direction on the filtration page. The boundary joining the two pairs, earlier physical images and the actual corrected second pairing can still remove it. The original index allows source lengths 11, 13 and 15. [A8 rank files](5701203/receipt.json) · [Exact rational weight calculation](q-8-0.json)

## Taylor projection before taking ranks

Use the exterior chains of positive polynomial currents in two ordinary Taylor variables $x,y$. The classical boundary brackets two colours and multiplies the Taylor monomials. Fix the colour, total Taylor degree and exterior occupation of each set of identical letters. Within that subspace the lowering $y\partial_x$ on $x^a y^{d-a}$ has coefficient $a$. The colour and exterior multiplicity are preserved. Quotient the weight $(5,5)$ by its image from $(6,4)$. In a rational irreducible $GL(2)$ representation this quotient keeps exactly the determinant $\det^5$. The same finite calculation of the quotient is reused whenever the jet degrees and exterior occupations agree. Each projected boundary has integer entries after multiplying by $10!$. Both adjacent matrices compose to zero before reduction modulo any prime. [Exterior chain bases](currents.py) · [Taylor projection](determinant.py)

For a dominant colour weight, quotient the projected chains by the sum of the two simple colour lowering images into that weight. Complete reducibility in characteristic zero identifies the quotient with the multiplicity space of that highest weight. The Weyl character formula gives its dimension independently at all three lengths. The boundary commutes with the lowering maps. Its rank on the quotient equals the rank after adjoining both lowering images to its columns, minus the dimension of the colour relations in the target. The latter is the full target weight dimension minus its rational highest-weight multiplicity. A modular rank of the augmented columns therefore gives a lower bound on the rational boundary rank. The calculation requires no semisimplicity assumption in either finite field. [Augmented matrices and Weyl dimensions](colour_quotient.py)

For example, at $(6,0)$ the two relation dimensions are $14668-649=14019$ and $17376-543=16833$. Both primes give augmented ranks $14374,17021$. Thus the rational homology dimension is bounded by
$$
543-(14374-14019)-(17021-16833)=0 .
$$

A second calculation uses integer row elimination on the weight-space boundaries followed by the rational Weyl character formula. It gives multiplicities zero and one in $(4,2)$ and $(8,0)$. The alternating sum for $(8,0)$ is $63-35-49+3+19=1$. Every other term vanishes. [Rational ranks for (4,2)](q-4-2.json) · [Rational ranks for (8,0)](q-8-0.json)

## Replay

Python 3.10 or later suffices for the two rational replays. The augmented route additionally uses NumPy and a CBLAS library. Its matrix panels use primes 127 and 131. Every intermediate in a BLAS dot product has absolute value at most $128(131-1)^2+131=2163331<2^{24}$, so binary32 represents all integers encountered exactly. The scalar elimination and the reduction of each panel use integer modular arithmetic. [Panel implementation](panel_rank.cpp) · [Sparse rational elimination](exact_rank.py)

~~~sh
python3 verify.py rational 8 0 | diff -u q-8-0.json -
python3 verify.py rational 4 2 | diff -u q-4-2.json -
c++ -std=c++17 -O3 -shared -fPIC panel_rank.cpp -lopenblas -o panel_rank.so
python3 verify.py augmented 6 0 | diff -u 5701203/result-60.json -
python3 verify.py augmented 4 2 | diff -u 5701203/result-42.json -
python3 verify.py augmented 8 0 | diff -u 5701203/result-80.json -
~~~

On macOS the compile command can use `-framework Accelerate` in place of `-lopenblas`.
