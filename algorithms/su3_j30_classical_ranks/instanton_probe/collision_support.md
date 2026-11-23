# Scalar support of the level-30 SU(3) class

On the reduced Cartan of the three undifferentiated scalars, the classical R-007 class has support precisely where two **joint eigenvalue vectors** coincide. Here a background belongs to the support if every polynomial in scalar graviton traces that annihilates the class vanishes at that background. Six explicit scalar discriminants detect the complement. Every choice of five directions misses a background with distinct joint eigenvalues. The charge-one supercurrent matrix element remains open.

## Annihilators in every scalar direction

Let $F$ denote the public 6,139-term source. It has flavour weight $(7,4,4)$. The exact graded trace actions of the raising operators and of the two shears on the last two flavours make $F$ the highest vector of a symmetric cubic for $SL_3$. Consequently a scalar direction $v$ determines a cubic $F(v)$ in its three coordinates, with $F(e_1)=F$ for the original first-flavour direction $e_1$. Let $\Delta(v)$ be the discriminant of the same combination of the three scalar matrices. Flavour covariance of the polynomial classical primitive gives

$$
\Delta(v)^6[F(v)]=0 .
$$

The checks on the full source use both the scalar and the dual matter-fermion flavour actions, including the graded cyclic signs. [Classical primitive](discriminant_ward.md) · [Exact flavour actions](discriminant_pencil.py)

To keep the external operator fixed as the probe changes, substitute $v+t e_1$. Write $g=\Delta(v)$ and $z=\Delta(v+t e_1)-g$, which is divisible by $t$. For a positive integer $m$,

$$
(g+z)^m
\left[g^3-mg^2z+\binom{m+1}{2}gz^2-\binom{m+2}{3}z^3\right]
=g^{m+3}\pmod {t^4}.
$$

The coefficient of $t^3$ in $F(v+t e_1)$ equals $F$. Taking $m=6$ gives $\Delta(v)^9[F]=0$ for every $v$, including singular directions, without dividing by a scalar trace.

For an $SL_3$-equivariant perturbative completion, the existing single-direction exponent is $70$. Such a completion exists by projecting any homogeneous completion onto the symmetric-cubic isotypic summand in each finite charge space. The perturbative differential commutes with this projection. Apply the same identity in the quotient by the subalgebra of completed gravitons, regarded as a module for the scalar graviton dressings using the marked product. The coefficient of $t^3$ gives $[\Delta(v)]^{73}[\widehat F]$ in that subalgebra for every $v$. The exponents are sufficient. [Single-direction exponent and assumptions](https://github.com/lbotao399-stack/quantum-lifts-results/blob/fb9d3c6/reviews/independent-g2-and-r007/np-discriminant/README.md)

## Six projected discriminants detect joint collisions

Take $v=(1,t,t^2)$ at $t=0,1,2,3,4,6$. For three distinct joint eigenvalue vectors, each of their three pair differences pairs with $v$ to give a nonzero polynomial of degree at most two. If the six discriminants all vanished, each polynomial would have two of these roots, and the three root pairs would be disjoint. The coefficient row for a pair $u,w$ would be a nonzero multiple of $(uw,-u-w,1)$. The oriented eigenvalue differences sum to zero, so the determinant of the three rows would vanish. The exact determinants for all fifteen pairings, in the replay order, are in the certificate. Each is nonzero. Conversely a joint collision makes every projected discriminant vanish. [Fifteen determinants and traceless five-probe control](collision_support.json)

To see that five probes can never suffice over the complex numbers, regard their directions as points in the projective plane. If no three are collinear, draw the lines through two disjoint pairs and join their intersection to the fifth point. If three or more are collinear, take their line, one line for the remaining points, and a third line through the intersection. If all are collinear, choose any two further lines through a point on that line. These constructions cover the five points by three distinct concurrent lines. Their nonzero normals admit a linear relation with all three coefficients nonzero. After rescaling, their sum is zero. They are therefore the three nonzero successive differences of a traceless triple of eigenvalue vectors. All five projected discriminants vanish on this triple.

## Detection on the collision locus

On the pair-collision block replace each scalar by $\operatorname{diag}(X^i+h_i\mathbf 1_2,-2h_i)$, with $X^i$ traceless, and restrict the fluctuating fields to $SU(2)$. Central shifts commute with the block. Restriction intertwines the classical supercharge, which raises fluctuation length by one. Use Dot's rational functional in the balanced $(4,4,4)$ block at fluctuation length seven. Its pairing on every one of the 2,808 incoming columns is zero. The parent has flavour $(7,4,4)$ and length ten. Only removal of three first-flavour scalar slots reaches this functional. Hence Dot's one-variable block evaluation gives the entire polynomial on arbitrary $h$ in these coordinates,

$$
\left\langle\ell,F\big|_{\rm block}\right\rangle
=-\frac{5184}{5}\,i\,h_1^3\,(-30)=31104\,i\,h_1^3 .
$$

Here $\ell$ is the archived functional, its value on the 117-term $SU(2)$ class is $-30$, and the colour coordinates are the original unscaled Pauli matrices. The replay regenerates the block coefficient and checks all incoming columns. [Dot's block and normalization](https://github.com/lbotao399-stack/quantum-lifts-results/tree/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-exact-input) · [Replay](collision_support.py)

Suppose a scalar polynomial annihilates $[F]$ and has nonzero value on this background. Its restricted series is invertible in the completion by $SU(2)$ fluctuation length. Its inverse is classically closed. Multiplying the proposed boundary by that inverse makes the restricted $F$ a boundary. The displayed functional excludes this for $h_1\ne0$. Every annihilating scalar polynomial thus vanishes on the dense open $h_1\ne0$ of the collision component. Polynomiality extends vanishing to the entire component. Permuting the three colour eigenvalues covers all pairs. Together with the six annihilators this proves equality of reduced supports.

For a proposed perturbative relation taking a scalar multiple of $\widehat F$ into the completed graviton subalgebra, specialization at zero loop parameter gives a classical graviton polynomial modulo a boundary. Each symmetrized single-trace colour tensor on the block with fixed central scalar slots is an invariant symmetric tensor of the fluctuating $SU(2)$ colours. Such tensors have even degree, since the symmetric invariant algebra of the adjoint is generated by the quadratic. Taylor derivatives preserve colour arity. Every restricted graviton polynomial, scalar dressing, and the inverse series of a scalar dressing with nonzero constant term thus has only even fluctuation lengths. The length-seven functional still kills its graviton part. This gives the same lower support bound for equivariant completions modulo gravitons over the formal power series ring, using scalar dressings with these classical leading representatives. It does not compute a nonperturbative differential.

Run the following in this directory from a repository with the cited Git objects.

~~~sh
python3 -B collision_support.py | diff -u collision_support.json -
~~~
