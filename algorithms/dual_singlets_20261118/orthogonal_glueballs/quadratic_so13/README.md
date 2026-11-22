# Both SO(13) quadratic targets survive the second differential

In the full relative holomorphic jet complex of undeformed $\mathcal N=4$ SYM, the incoming second differential at dotted and flavour weight $(7,7;2,0,0)$ has rank zero. At level $j=46$ the classical target has dimension three, the incoming first image has dimension one, and both quotient classes survive to $E_3$. The order counts supercharge obstructions. The known third differential on the scalar multiple of the seventh glueball power still applies. [First differential and the strict substitution](../scalar_multipliers/README.md) · [Full jet complex and physical comparison](https://arxiv.org/abs/2611.24160)

Use the defining trace $S=\operatorname{Tr}_{13}(\bar\lambda_{\dot+}\bar\lambda_{\dot-})$ and put $A=\phi^1$. The preceding certificate orders the three trace coordinates as

$$
S^7\operatorname{Tr}(A^2),\qquad
S^6\operatorname{Tr}(\bar\lambda_{\dot-}A)\operatorname{Tr}(\bar\lambda_{\dot+}A),\qquad
S^6\operatorname{Tr}(\bar\lambda_{\dot+}\bar\lambda_{\dot-}A^2).
$$

The first image is the line $(1,-2,-18)$. Two strict odd substitutions give the covector $(2,1,0)$ times $S^8$, known to survive to the third page. The independent covector $(18,0,1)$ annihilates the first image. Its first transpose has 71 terms and its second transpose has 2,809. The latter satisfies the required equation on all 15,928 potentially active native columns. These are first-jet endpoints, and an incoming cochain at the bare second-supercharge order has one second jet. Three selected letters supply at most four of the six required flavour slots. The bare second term therefore vanishes. [Transpose construction](../scalar_multipliers/second.py) · [Pure eighth-power certificate](../so13_j48/README.md)

## Finite determinant identities at the second transpose

The ordered transpose must vanish on finite-colour relations after subtraction of a functional on the stable graded-commuting quotient. At second order its weight is $(5,5;4,2,2)$. The orthogonal second fundamental theorem supplies determinants with fourteen left and fourteen right endpoints among 36 indices. At most four further contraction edges remain. The free endpoints of a determinant can be labelled by two auxiliary directions and polarized by a Pfaffian. The fixed numbers of the two odd first jets are five each. [Orthogonal invariant tensor relations](https://arxiv.org/abs/1102.3221) · [Polarized determinant coefficient](../../orthogonal_glueball_relation_span.py)

Sorting the letters along an open path or closed trace changes the relation by $Q_0$ of a finite-colour identity. The second row on this change equals the first row on its $Q_1$ image, up to the fixed descent sign. That image is another finite-colour identity, and the first row descends by the first-page calculation. A sorted path contains at most one gaugino of either sign. A path with $l$ letters costs $l-1$ ordinary edges if open and $l$ if closed. A closed trace has even length. Transposition leaves one identity auxiliary matrix at even open length or three traceless axes at odd length. At most three single open paths of each dotted sign survive exterior antisymmetry.

The resulting 46 types consist of 33 with no closed trace, seven with a closed two-gaugino spectator, four with a closed gaugino–scalar spectator, one with a closed scalar bilinear, and one with two closed two-gaugino spectators. In the last three cases the spectator must be multiplied into the ordered relation before its flavour projection. The verifier strips its specified factors from the ordered probes with their Koszul signs, and distributes the spectator among labelled planes for the commuting test.

The target selects the flavour representation with diagram $(4,2,2)$. For the three-by-three scalar marker matrix this is two determinants and a symmetric pair of derivatives of the highest flavour. The scalar singleton matrices can be chosen as a basis of the three-dimensional space of traceless auxiliary matrices. Invertible frames are dense, and the marker polynomial recovers singular choices by polarization. Give each internal or closed scalar slot its own covector on that frame. Let each traceless open-path matrix and each of the two remaining flavour derivatives vary over all three axes. This generates the required flavour projection with arbitrary internal scalar insertions. Repeated paths commute or anticommute according to their physical parity.

There are 1,848 auxiliary weight-zero components. Simultaneous conjugation of all auxiliary matrices preserves the relation. With the basis $H,E,F$ of $\mathfrak{sl}_2$, its raising action sends $H$ to $-2E$ and $F$ to $H$, with the negative transpose on the scalar covectors. At each type the verifier forms the full raising matrix from weight $-1$ to zero, in units of a root. It quotients by its image using rational row reduction. This leaves 208 representatives. The replay can also evaluate discarded components directly. [Weight matrices](auxiliary.py) · [Closed and open factors](relations.py)

## Exact descent

There are 1,327 sorted trace monomials. Coefficients on distinct two-dimensional Cartan planes vanish on every $Q_0$ boundary. All finite relations vanish on six or fewer planes. The 953 remaining plane coordinates include an invertible 20-column minor on the determinant rows. Appending the ordered second pairing leaves the rank equal to 20. The program solves the square minor and checks the resulting rational functional on all 208 rows. Subtracting it from the ordered second transpose annihilates the full required finite-colour subspace. Applying the flavour projection to the dual commutes with the two supercharges. Both independent target covectors now annihilate $d_2$. Positive-order supercharges annihilate these target representatives by their flavour support. [Rational rows](rows.json.gz) · [Functional and pivot indices](result.json)

As an independent check, project the three previously transported families to the same flavour component using the quadratic $\mathfrak{gl}_3$ Casimir. Their unprojected trace rank is 23. After the projection they have rank 15. The direct determinants and the projected transported families together have rank 20 on all 953 plane coordinates. [Independent transport test](transport_check.py)

The faster trace-log calculation fixes one physical slot on a cyclic trace. For a cyclic tiling with $t$ matrix tiles, the $t$ possible marked tile boundaries cancel the $1/t$ in the logarithm. Divide by the rotational stabilizer of the colour word and include its reversed necklace when distinct. The native replay uses the original sum over every rotation and reversal. It independently reproduces all 20 pivot rows and the same rational functional on one discarded auxiliary component in each of the 46 path types. [Rooted trace calculation](rooted_trace.py) · [Native trace calculation](../../orthogonal_glueball_second_sources.py)

## Antisymmetric scalar partners

The same determinant test gives twelve nonzero first pairings in the highest dotted-spin-one, antisymmetric-flavour component. For $SO(N)$ use $S$ in the defining trace and set the terminal power to $p$. The unique stable coordinate contains $p$ bilinears and one trace of each of $\bar\lambda_{\dot+}\phi^1,\bar\lambda_{\dot+}\phi^2$. Its level is $6p+10$. The replay uses the native sorted trace convention. For the colour range in the table there are too few ordinary edges in the classical relation to accommodate the excess positive gauginos, so the coordinate is nonzero already at finite colour.

The source relation has $p-2$ marked two-gaugino open paths, three distinct positive single axes, one negative single axis and five scalars. Project the latter with a scalar-marker determinant and antisymmetrize the two remaining derivatives over flavours 1 and 2. Its Pfaffian degree is $p+7>N$. The total trace length $2p+5$ is odd, so its stable commuting section vanishes. Its graded-reversal primitive becomes closed in $N$ colours. The coefficient obtained from that primitive by the ordered first supercharge follows. Each transposed equation covers 69 native columns. [Source and replay](triplet.py) · [All coefficients](triplet.json)

| $N$ | $p=N-6$ | $p=N-5$ |
|---:|---:|---:|
| 7 | | 12 |
| 8 | 0 | -36 |
| 9 | -12 | 96 |
| 10 | 48 | -240 |
| 11 | -144 | 576 |
| 12 | 384 | -1344 |
| 13 | -960 | 3072 |

## Two infinite first-differential families

Take $m=2$ or $3$ distinct scalar flavours, each in its own trace with a positive dotted gaugino, multiplied by $S^p$. Here $k=p+m$ is the pure glueball power obtained by replacing each scalar by a negative dotted gaugino. Denote the first-image coefficient in the unique sorted target coordinate by $c$. In the source choose three distinct positive singleton axes, $3-m$ distinct negative singleton axes, $k-4$ labelled double paths, and the scalar projector specified below. For $k\ge4$ and $N=k+3,k+4$ the coefficient satisfies the exact identity

$$
c^2=
\begin{cases}
[12\,2^{k-5}(2N-3k-4)]^2,&m=2,\\
[72\,2^{k-5}(2N-3k-4)]^2,&m=3 .
\end{cases}
$$

To see this, apply one strict flavour-to-negative-dotted substitution for each scalar flavour. Each acts on its unique target trace. The output is $S^k$ with coefficient of absolute value one. The substitutions anticommute with the tree charge. At these grades the only active first-order pair is a matter second jet and its complementary first scalar. Substitution on that pair or on its two outputs vanishes; the spectator terms anticommute by parity. The signed product of the substitutions therefore intertwines the first obstruction with the pure one. Since $k<N-2$, the classical glueball theorem gives $[S^k]\ne0$ and certifies finite-colour survival of the mixed classical target. [Strict substitutions and checks](marker_norms.py) · [Classical glueball theorem](https://arxiv.org/abs/math/0312452)

On the finite identity, the same substitutions shift the negative single-ghost markers by the corresponding rows of the scalar marker matrix. Antisymmetry forces the new axes to fill exactly the $m$ vacant slots. The scalar polynomial supplied by this wedge is the minor of rows 1 and 2 and columns $E,F$ for $m=2$, or the full determinant for $m=3$. Write $D$ for the scalar-marker determinant and $M$ for this minor. Projection of the residual three scalar slots to the flavour singlet leaves only $D$. With the ordinary polynomial-differentiation pairing, distinct monomials are orthogonal and a monomial has norm squared equal to the product of the factorials of its multiplicities. The ratio of the new singlet extraction to the old one is $\lVert D M\rVert^2/\lVert D\rVert^2$. Direct expansion gives $\lVert D\rVert^2=6$, $\lVert D M\rVert^2=24$ for $m=2$, and $\lVert D^2\rVert^2=144$ for $m=3$. Labelling the $k-4$ identical double paths supplies $(k-4)!$. Their reversal and moving exterior slots contribute only signs. This proves the factor $4(k-4)!$ or $24(k-4)!$ in absolute value.

Two stable primitives with the same classical derivative differ by a stable cocycle. Its ordered first pairing with the pure target is polynomial in $N$, since the elementary index reconnections only insert factors of the defining dimension. For every sufficiently large $N$ the pure power is nonzero on the commuting Cartan and cannot be a perturbative boundary by the Cartan-image theorem. Hence that polynomial vanishes identically. The first coefficient is therefore independent of which primitive of this derivative is used. Substitution of the pure coefficient $3\,2^{k-5}(2N-3k-4)/(k-4)!$ proves the squared formula. It eliminates both families as second-differential targets on these two bands for all $N\ge9$. [Pure first coefficient](../first_family.md) · [Cartan-image theorem](../../../perturbative_cartan_targets.md)

For three flavours the native trace calculation covers $7\le N\le16$ and $p=N-7,N-6$ whenever $p\ge1$. The level is $6p+15$ and the dotted-spin highest weight is $3/2$. Eighteen entries are nonzero. Together with the two-flavour spin-one rows above, this gives 30 explicitly evaluated nonzero target multiplets. The $SO(7)$ entries at levels 21 and 22 fall inside the earlier complete census and count as controls. The remaining 28 are new to the campaign ledger. The ledger counts only the explicitly tested multiplets. [Signed three-flavour coefficients](quartet.json) · [Earlier census](../../../../data/sdual_SO7_Sp3.txt)

## Replay

With the repository's Python trace dependencies, run from this directory.

~~~sh
python3 replay.py --rebuild --workers 4
python3 replay.py --native --components sample --workers 4
python3 transport_check.py
python3 triplet.py
python3 marker_norms.py
~~~
