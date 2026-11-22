# The first differential into two powers of the orthogonal glueball

Let $N$ be the defining dimension of $SO(N)$ and $S=\operatorname{Tr}_N(\lambda_+\lambda_-)$. For $p\ge4$, take the existing Pfaffian source with one scalar-marker determinant, six single ghost paths, and $p-4$ double ghost paths. Its classical derivative is a polarized Pfaffian coefficient of degree $3+6+(p-4)=p+5$. For $N<p+5$ this derivative vanishes after specialization to $N$ colours. The first differential on the resulting source $O_0$ in the full relative holomorphic jet complex is

$$
d_1[O_0]=\frac{3\,2^{p-5}}{(p-4)!}(2N-3p-4)[S^p].
$$

The normalized odd-trace reversal specifies $O_0$ in the earlier source prescription. At the output weight every word has $2p$ first gaugino jets. A classical predecessor has exactly one second dotted jet, whose images generate the three classical gaugino anticommutator relations. The classical chiral-ring theorem for $\mathfrak{so}_N$ gives $[S^p]\ne0$ for $p<N-2$ and $[S^{N-2}]=0$. This justifies detecting the target after finite-colour specialization. [Source, orientation and trace conventions](README.md) · [Classical chiral-ring theorem](https://arxiv.org/abs/math/0312452) · [Known $SO(7)$ first pairing and higher-rank question](https://arxiv.org/html/2610.22332v3#S5)

For $N=p+4$, the coefficient is $-3\,2^{p-5}/(p-5)!$ for $p\ge5$, and zero for $p=4$. For $N=p+3$, it is $-3(p-2)2^{p-5}/(p-4)!$ for every $p\ge4$. Thus for every $N\ge9$ both $S^{N-4}$ and $S^{N-3}$ are already boundaries on the second page. The one-source zero at $(N,p)=(8,4)$ gives no exclusion by itself.

## Coefficient calculation

In the classical quotient an even trace with a repeated odd letter is a boundary, and a trace of odd length is a boundary by graded reversal. The sole surviving ordered monomial is the product of $p$ bilinears. Joining different traces under $Q_1$ would require $p+1$ traces in a source of length $2p$, with at least two letters per trace. Thus only a single three-letter active trace and twelve five-letter orientations enter the calculation. Their reversal primitives give two-letter and four-letter inputs to $Q_1$.

Use the marker order $+H,+E,+F,-H,-E,-F$ of the published Pfaffian. In the table each mask is the sum of the six bit values $1,2,4,8,16,32$ for the selected odd markers. The last column multiplies one copy of the commuting double-path marker. Entries already include scalar-determinant extraction for the active traces.

| Trace contribution | No odd markers or double path | Mask 9 | Mask 20 | Mask 34 | One double path |
|---|---:|---:|---:|---:|---:|
| One spectator $\operatorname{Tr}(-+)$ | 0 | 1 | $1/2$ | $1/2$ | $-2$ |
| Three-letter active trace after $Q_1$ | $-12(N-2)$ | 0 | 0 | 0 | 0 |
| Five-letter active traces after $Q_1$ | 0 | 24 | 12 | 12 | $-36$ |

Literal Pfaffian minors on disjoint colour cycles give the table's marker entries. The contraction of two elementary skew generators gives the last two rows by the defining-representation commutator. A one-bilinear output permits at most one empty index loop, so the colour dependence is affine in $N$. A two-bilinear output uses both traces and permits no empty loop. Exact calculations at $N=7,11$ therefore fix the entire colour dependence. The verifier checks every active orientation against both implementations. [Direct contractions and polynomial check](first_family.py)

The product of the three disjoint two-marker terms in a spectator cube gives the six-marker coefficient $3!(1)(1/2)(1/2)=3/2$ in the stated order. An identical spectator repeated $k$ times enters the Pfaffian exponential with denominator $k!$. Collecting the two possible trace outputs, with the global odd-marker sign and the reversal $\operatorname{Tr}(-+)=-S$ restored, gives the following two coefficients of $[S^p]$,

$$
\begin{aligned}
&\frac{2^{p-4}}{(p-4)!}\frac{12(N-2)}{3!}\frac32
  =\frac{2^{p-4}}{(p-4)!}\,3(N-2),\\
&\frac{2^{p-4}}{(p-4)!}
  \left(-\frac{24}{2!}\frac32-\frac{36}{2\cdot3!}\frac32(p-4)\right)
  =\frac{2^{p-4}}{(p-4)!}\left(-18-\frac92(p-4)\right).
\end{aligned}
$$

The term containing $p-4$ chooses one double path from the five-letter active trace. Summing the two rows gives the stated coefficient, including the $p=4$ endpoint. The verifier also evaluates the original ordered-trace formula on both diagonals for every $4\le p\le14$.
