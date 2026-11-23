# Classical products in the second differential into the sixth glueball

In the full relative holomorphic jet complex of $F_4$ at $j=36$, the second supercharge differential into $S^6$ vanishes on the part of its source spanned by products of positive-degree classical cohomology classes. Here $S$ is the native glueball and the classical page is $E_1=H(Q_0)$. Thus a nonzero source must survive the quotient of classical cohomology by these products. A cycle expressed as a sum of products of cochains whose factors fail to close can survive this quotient. The existence of such a source with nonzero second pairing remains open. The argument uses the native Taylor Ward identities in this sector. The comparison with physical local operators uses the all-loop manuscript's Appendix B. [Sector and Ward identities](README.md#a-single-taylor-representation-can-hit-the-sixth-glueball)

## Which products can contribute

The source has fourteen Taylor slots in thirteen letters. Its only contributing Taylor diagram is $(4,4,3,3)$. Exactly one factor in a classical product has a second jet. Every other factor contains only first jets. Such a factor is closed under the first supercharge correction since a pair of first jets carries at most two flavour slots and the correction requires three. Combine these factors into one terminal spectator. The spectator has at least two slots since the simple algebra has no invariant linear form.

If the spectator has at least six slots, the factor with the second jet has at most eight. The first differential shifts its Taylor highest weight by $(1,1\mid-1,-1,-1)$. A polynomial target requires three positive odd coordinates before the shift. Its source diagram must therefore contain a square of three rows and three columns. Eight slots cannot contain this square. Both factors represent cycles on $E_2$, allowing a zero class if a factor is a first boundary.

For spectators of degrees two through five the defining traces span the invariant tensors. Modulo classical boundaries the first-jet matrices graded commute, since their graded commutators are the variations of the corresponding second jets. Orthogonal reversal kills each odd-length trace in this quotient. The polarized identity $12\operatorname{Tr}_{26}(X^4)=(\operatorname{Tr}_{26}(X^2))^2$ expresses every degree-four trace as a product of quadratics. It suffices to test a quadratic spectator and any classical cycle in the remaining twelve slots. This includes an arbitrary polarization of the quadratic graviton. [Native invariant dimensions and quartic identity](../mixed_colours.md)

The two-slot quadratic has Taylor diagram $(1,1)$. Removing a vertical two-strip from $(4,4,3,3)$ gives precisely the following shorter diagrams. An independent symmetric-group character sum gives multiplicity one in each case.

| Shorter diagram | Eleven-letter cochains | First target diagram | Classical obstruction to a first image |
| --- | ---: | --- | --- |
| $(4,4,2,2)$ | 750 | empty | The third odd coordinate becomes negative. |
| $(3,3,3,3)$ | 141 | $(4,4,3)$ | The shortened source has dimension $141-80-61=0$ in cohomology. |
| $(4,3,3,2)$ | 791 | $(5,4,2)$ | The terminal target has zero classical cohomology by the calculation below. |

Every shorter component that can reach the full diagram is therefore a first-differential cycle. Projection to the Taylor diagram commutes with the induced differential up to its stated Berezinian twist. Every product under consideration is consequently a sum of projected products on $E_2$. [The two ranks at $(3,3,3,3)$](../../../bmn_primary/f4_scalar/README.md#application-at-level-24) · [Character sum and all three shifts](glueball_factor_ward.json)

On $E_2$ the differential is a derivation. For a product to hit the zero-flavour, zero-excess target, its unchanged positive-degree factor must already have these two degrees. The invariant two-gaugino ring then makes that factor a multiple of $S^k$ for $k\ge1$. A nonzero result would require a second image proportional to $S^{6-k}$. Multiplication by the closed $S^{k-1}$ would make $S^5$ a second boundary, contradicting its perturbative survival. This also applies after taking weight components and their linear combinations. [Survival of the fifth power and the ring argument](README.md#one-row-decides-the-second-differential)

## A rational eight-dimensional test

The terminal cochain multiplicity for $(5,4,2)$ is 55. Its highest Taylor weight has five positive dotted gauginos, four negative dotted gauginos and two antisymmetrized scalars. A functional on the classical quotient kills every bracket of two first jets. Killing the same-species brackets separately places the three colour legs in the top-Casimir components of exterior degrees five, four and two. In Dynkin labels the first leg has $(0,1,0,3)$ and $(0,0,3,0)$, the second $(0,0,2,1)$ and the third $(0,1,0,0)$. The multiplicities of invariant trilinears are respectively five and three. [Native exterior modules](dotted_filtration.py) · [Weyl tensor multiplicities](intermediate_rank.py)

Fix the highest wedge of the first leg and the lowest wedge of the second. Weight conservation leaves 11 and 14 possible scalar pairs for the two choices of highest wedge. These values determine an invariant functional. Indeed opposite highest and lowest flags form an open dense diagonal orbit in the product of the two flag varieties. The translates of either extremal vector span its irreducible module. The replay also reconstructs the functional by raising the first wedge and lowering the second with the native brackets. For a weight $w$ on a degree-$s$ exterior leg the two Casimir denominators are $18s-(w,w+2\rho)$ and $18s-(w,w-2\rho)$, where $2\rho$ is the sum of the positive roots and the pairing uses the archived $B=\operatorname{Tr}_{26}/6$. Moving the outer root operator to the other tensor legs contributes the minus sign in the Ward identity. The code checks the extremal weight multiplicities directly in the full exterior powers.

The scalar bracket and the common stabilizer impose 17 independent equations on the 25 coordinates. Their rational kernel has dimension eight, equal to the independent Weyl-character count. In the free-coordinate basis printed in the certificate, two negative-gaugino–scalar brackets, five positive-gaugino–scalar brackets and one bracket of opposite gauginos give an eight-by-eight determinant

$$
\frac{523032987223594133}{125439321093750}.
$$

The unreduced 25-by-25 determinant is $523032987223594133/83626214062500$, with residues $806$ modulo $1009$ and $936$ modulo $1013$. Every functional annihilating classical boundaries therefore vanishes. Schur decomposition gives classical rank $55$ on the selected terminal component. The same component is absent from the perturbative terminal operators with the physical charge of $S^5$, where the dotted and flavour degrees are $(4,4\mid1,1,1)$. [All rational rows, basis and minors](glueball_factor_ward.json) · [Native Ward replay](glueball_factor_ward.py)

An independent forward calculation polarizes eleven commuting-variable slots down three columns of heights $5,4,2$. At odd total length defining traces are classical boundaries by reversal. A pair of Albert multiplication-operator projections on words of lengths five and six is also a boundary. Reversing the length-five word changes its sign since its distinguished multiplication operator is self-adjoint and each Lie input is skew. Subtracting the reversed word expands into brackets with second-jet primitives. The 51 independent defining traces and two additional Albert contractions attain rank 53 directly. The local and A8 computations used different primes and integer assignments. [Defining-trace recipes](glueball_factor_traces.json) · [Albert replay](glueball_factor_jordan.py) · [Forward and cluster results](glueball_factor_receipts.json)

Run from this directory with Python, NumPy, SciPy and python-flint.

~~~sh
python3 glueball_factor_ward.py --depth 0 --out /tmp/f4-ward.json
cmp glueball_factor_ward.json /tmp/f4-ward.json
python3 glueball_factor_test.py --prime 1013 --replay glueball_factor_traces.json --out /tmp/f4-traces.json
python3 glueball_factor_jordan.py --traces glueball_factor_traces.json --prime 1013 --sample-start 852179 --seeds 96 --points 62 --out /tmp/f4-albert.json
~~~
