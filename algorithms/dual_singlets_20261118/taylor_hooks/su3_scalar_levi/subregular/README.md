# Scalar charts for the $SU(3)$ supercharge

For the archived levels 32 and 34, the first supercharge equation is locally soluble on a neighbourhood of every pair of independent square-zero scalars with a common image line or a common kernel plane. The remaining possible support of either first obstruction consists of triples proportional to one rank-one nilpotent matrix, including the zero triple. The same bound applies to the second residual at level 34 once a global first correction has been found. The primitives on these open sets belong to the localized matrix-coordinate complex. The global equation requires a gauge-invariant polynomial primitive. [Fixed-scalar calculation](../square_zero/README.md) · [Simple-eigenvalue calculation](../README.md#arbitrary-matrices-above-a-scalar-with-a-simple-eigenvalue)

## Two independent scalars

In the Taylor superfield write $\Psi=A+\theta^1 B+\theta^2 C+\theta^1\theta^2 D$. The remaining odd coordinate and the two spacetime directions are free. The fields $B,C$ are even, $A,D$ are odd and $A(0)=0$ in the remaining directions. The source convention $Q_0\Psi=\Psi^2$ gives

$$
Q_0 A=A^2,\qquad Q_0 B=[A,B],\qquad Q_0 C=[A,C],\qquad Q_0D=AD+DA+[B,C].
$$

Filter by the number of $A$ factors. The filtration is finite at fixed remaining Taylor degrees since $A(0)=0$. Its first differential sends $D$ to $[B,C]$ and kills the other three fields. On the chart where $T=(v,Bv,Cv)$ has invertible constant determinant for a fixed vector $v$, conjugate $B,C,D$ by $T$ and choose coordinates with $v$ first. The first columns of the transformed $B,C$ are the second and third basis columns. On their commuting locus their full matrices, with $p,q,r,t$ the entries in positions $22,32,23$ of $B$ and $23$ of $C$, are

$$
B=\begin{pmatrix}0&2p^2+2qr&pr+qt\\1&p&r\\0&q&-p\end{pmatrix},\qquad
C=\begin{pmatrix}0&pr+qt&2r^2-2pt\\0&r&t\\1&-p&-r\end{pmatrix}.
$$

Off this locus six further coordinates are the entries $11,21,31,22,32,23$ of $[B,C]$. Substitution gives a polynomial inverse for these six coordinates and sets the entire commutator to zero when they vanish. The two remaining independent entries of $D$ can be replaced by $\operatorname{tr}(DB),\operatorname{tr}(DC)$, which are already closed for the full $Q_0$ by graded cyclicity. Their coefficients on the two unused entries $D^1{}_2,D^1{}_3$ form the identity. The first differential therefore splits off six coordinate pairs at every Taylor coefficient. Each pair has scalar cohomology by contraction with its degree-counting operator. [Symbolic inverse, trace Jacobian and all 180 source signs](charts.json) · [Symbolic verification](../../../su3_scalar_charts.py)

For the induced differential write $\gamma$ for the terms containing one additional $A$. In the moving frame the odd matrix $T^{-1}AT-T^{-1}\gamma T$ has columns $w,Bw,Cw$, where $w$ is its first column. Since $B,C$ commute and $v,Bv,Cv$ span, this matrix is a linear combination of $1,B,C$ with odd coefficients. It commutes with $B,C$, so $\gamma$ kills their four displayed coordinates. Its square and its $\gamma$ variation both vanish. Its scalar coefficient is fixed by tracelessness of the original $A$ and equals $-\operatorname{tr}(T^{-1}\gamma T)/3$. The two nonscalar coefficients and the six positive-degree entries of $\gamma T$ reconstruct the original $A$ degree by degree. After removing the pairs formed by the positive-degree entries of $T$ and $\gamma T$, this page is a free supercommutative algebra over the domain of the six constant frame coordinates and $p(0),q(0),r(0),t(0)$ with $\det T(0)$ inverted.

There are no further differentials. To prove this, first pass to the fraction field of that domain. Every base coordinate has a full closed lift since $A(0)=0$. The four closed traces $\operatorname{tr}B^2,\operatorname{tr}B^3,\operatorname{tr}(BC),\operatorname{tr}(B^2C)$ have Jacobian determinant $6\operatorname{disc}(B)$ in the displayed coordinate order, where the discriminant is that of the characteristic polynomial. At $(p,q,r,t)=(0,1,1,0)$ it equals $648$. Recursion in the positive Taylor degrees therefore lifts every commuting-pair coordinate using those four traces. The frame $(v,Bv,B^2v)$ is generically invertible: its determinant relative to $T$ is $q$. In that frame the transformed ghost, using the full $Q_0$ transformation of the frame, is an odd polynomial in the companion matrix. Its two nonscalar coefficients are closed. On the page they give the two ghost coordinates by changing the centralizer basis from $1,B,C$ to $1,B,B^2$. The other two generators lift by $\operatorname{tr}(DB),\operatorname{tr}(DC)$. Thus every generator survives over the fraction field. Freeness implies successively that each differential vanishes over the original domain. [Single-matrix companion identity](../square_zero/README.md#a-cyclic-scalar)

Multiplication by the nonzero polynomial $\operatorname{disc}(B(0))$ is injective on the resulting associated graded cohomology, hence on the localized classical cohomology by the finite filtration. Every Cartan-zero trace cocycle is killed in classical cohomology by a power of this polynomial. It is consequently exact on the entire chart, including the locus where the discriminant vanishes. This last implication uses the scalar discriminant theorem. [Discriminant identity and its use for page residuals](../../../../perturbative_cartan_targets.md)

For two independent rank-one square-zero matrices with common kernel plane, choose $v$ outside the plane. Their two values on $v$ span that plane. Negative transposition gives the chart on a common image line. The single-scalar companion chart covers a rank-two nilpotent scalar. The existing two-colour calculation handles a simple scalar eigenvalue for the specified residuals. What remains has the form $S_i=x_i N$ with $N^2=0$ and $\operatorname{rank}N=1$, or all three $S_i=0$. A nonzero projective line of such matrices is determined by its image line inside its kernel plane. These flags have dimension $2+1=3$. For each flag the three scalar entries take values in the one-dimensional space of maps from the quotient by the kernel plane to the image. The remaining locus thus has dimension $3+3=6$ in the 24-dimensional space of traceless scalar triples. [Classification of rank-one pencils](../square_zero/README.md#rank-one-scalar-pencils)

## The chart at a nonzero rank-one scalar

Now write $\Psi=A+\theta^1 B$, retaining both other odd coordinates in the Taylor series. Choose component coordinates in which the constant entry $B^2{}_1(0)$ is invertible. Take the frame with columns $v,Bv,e_3-kv$, where $v=e_1$ and $k$ is the entry in row two, column three of $B$ in the preliminary frame $(v,Bv,e_3)$. The resulting scalar matrix and the traceless part of $T^{-1}AT-T^{-1}Q_0T$ are

$$
B=\begin{pmatrix}0&p&u\\1&q&0\\0&r&-q\end{pmatrix},\qquad
A=\begin{pmatrix}
a&dp+gu&du-2fq\\
d&a+dq&f\\
g&dr-gq&-2a-dq
\end{pmatrix}.
$$

Here $p,q,u,r$ are the four indicated even matrix entries and $a,d,f,g$ are the odd entries in positions $11,21,23,31$ of $A$. The form of $A$ follows by setting entries $11,21,31,23$ of $[A,B]$ to zero; the linear system for the other four entries has determinant of absolute value one. The central odd trace removed from the moving ghost is fixed by $T^{-1}Q_0T$, has zero variation and drops out of both equations. The three entries of $Bv$ and $k$ recover the frame. Their four positive Taylor series pair with their own classical variations. Conversely these variations and the four displayed odd fields recover the original traceless ghost by the inverse frame change. The remaining four constant frame coordinates are closed. This gives a polynomial reduced differential on the chart.

| Field | $Q_0$ variation | Field | $Q_0$ variation |
|---|---|---|---|
| $p$ | $-2q(fr-gu)$ | $a$ | $-2qfg$ |
| $q$ | $fr-gu$ | $d$ | $fg$ |
| $u$ | $3au-(p-2q^2)f$ | $f$ | $3af$ |
| $r$ | $(p-2q^2)g-3ar$ | $g$ | $-3ag$ |

Substitution into $Q_0B=[A,B]$ and $Q_0 A=A^2$ checks every entry, with all six independent quadratic exterior products of the four odd coordinates. At spacetime bidegrees through $(1,1)$ the 124 independent Taylor generators square to zero under $Q_0$. The images of all 31 original letters agree coefficientwise in all 279 matrix entries, comprising 3504 nonzero polynomial coefficients. [Taylor coefficients and signs](chain.json) · [Independent comparison with the original letter rules](../../../su3_subregular_differential.py)

## First corrections at the rank-one point

Set the moving frame to the identity. At the expansion point the first scalar has its only nonzero entry in row two, column one, and the other letters vanish. Give each of the 124 remaining Taylor generators polynomial degree one. This degree measures the order of vanishing at that point. A lower bound for a trace is the minimum total degree over its closed colour paths. For a product of traces add the individual bounds. Expanding every retained product over the integer polynomial and exterior coordinates gives the following exact counts.

| Physical level | Full first-image words | Words of minimum allowed degree 8 | Nonzero degree-8 monomials | Nonzero degree-9 monomials | Terms of the degree-7 primitive |
|---:|---:|---:|---:|---:|---:|
| 32 | 331201 | 153915 | 397 | 39588 | 248 |
| 34 | 517468 | 200409 | 1255 | 124633 | 573 |

All lower degrees vanish word by word. The archives use $+Q_1\mathcal O_0$ at level 32 and $-Q_1\mathcal O_0$ at level 34. The primitive replay restores the sign at each level and cancels every term through degree eight in $Q_1\mathcal O_0+Q_0\mathcal O_1$. A6 job 7657430 independently expanded both complete input streams through degree nine. Both uncompressed polynomial hashes agree with the local computation. [Level-32 polynomial](j32_through9.json.gz) · [Level-34 polynomial](j34_through9.json.gz) · [A6 receipt](a6.json) · [Level-32 correction](j32_degree8.json) · [Level-34 correction](j34_degree8.json)

From the repository root, run the following commands with Python 3.10 or later and SymPy. The polynomial replay uses the standard library.

```sh
d=algorithms/dual_singlets_20261118
p="$d/taylor_hooks/su3_scalar_levi/subregular"
python3 "$d/su3_scalar_charts.py" --out /tmp/su3-charts.json
cmp /tmp/su3-charts.json "$p/charts.json"
python3 "$d/su3_subregular_differential.py" --out /tmp/su3-subregular-chain.json
cmp /tmp/su3-subregular-chain.json "$p/chain.json"
for j in 32 34; do
  python3 "$d/su3_subregular_image.py" prepare --level "$j" --degree 9 --out "/tmp/su3-packet-$j.json.gz"
  python3 "$d/su3_subregular_image.py" expand --input "/tmp/su3-packet-$j.json.gz" --workers 4 --out "/tmp/su3-image-$j.json.gz"
  python3 -c 'import gzip,sys; assert gzip.open(sys.argv[1],"rb").read()==gzip.open(sys.argv[2],"rb").read()' "/tmp/su3-image-$j.json.gz" "$p/j${j}_through9.json.gz" || exit 1
  python3 "$d/su3_subregular_correction.py" --degree 8 --image "$p/j${j}_through9.json.gz" --recipe "$p/j${j}_degree8_primitive.json.gz" --out "/tmp/su3-primitive-$j.json"
  cmp "/tmp/su3-primitive-$j.json" "$p/j${j}_degree8.json" || exit 1
done
```
