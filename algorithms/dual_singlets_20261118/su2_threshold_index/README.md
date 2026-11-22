# The index partners of the SU(2) threshold at level 24

The $SU(2)$ threshold $\mathcal O_0$ has $(E,J_L,J_R)=(19/2,5/2,0)$ and $(R_1,R_2,R_3)=(3/2,3/2,3/2)$. Its proposed use as a nonperturbative lifting target or source is excluded by the weak-coupling index and the unitary superconformal multiplet decomposition. Fortuity had already given this index argument. The calculation here evaluates the full nonabelian classical complex on **every** opposite-statistics rung and builds separate integer lowering matrices for the only surviving blocks. [Prior argument and its classification inputs](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex#L291-L295) · [Index and multiplet argument](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex#L785-L796)

The three charges kept by the index are $R_i+J_L=4$. The source is a flavour and dotted-spin singlet. In the full classical jet complex let $t$ be the number of each dotted derivative and $Y$ the letter number. The three flavour degrees are all $4-t$, so $0\leq t\leq4$. The letter charge formula gives

$$
2J_L=12-t-Y,\qquad 5-2J_L=Y+t-7,\qquad 1\leq Y\leq12-t .
$$

Odd $Y+t-7$ gives opposite statistics. The resulting 26 singlet blocks exhaust the tower at that parity. At the threshold itself the zero-jet singlet cohomology has classical dimension $11-4-6=1$ with the known explicit nonboundary. The only possibly nonzero classical superconformal lowering on this singlet removes one flavour slot. Its highest-weight target at flavour $(4,4,3)$ has dimension $11-6-5=0$. Dotted lowerings annihilate the zero-jet representative. Thus the threshold is a lowest-weight class. [Explicit threshold and rational classical pairing](../taylor_hooks/su3_bmn_sources/quartic.json) · [Classical generator action](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/scripts/oneloop/su-threshold-classes/abelian_lw.py)

In 23 opposite-statistics blocks the incoming and outgoing classical ranks saturate the full singlet dimension. In the other three the Cartan restriction is injective on classical cohomology. The numbers below are equal at primes 1000003 and 2147483647.

| $t$ | $Y$ | $5-2J_L$ | cochain dimension minus incoming and outgoing ranks | Cartan rank | Cartan monomials before imposing the lowest-weight equations |
|---:|---:|---:|---|---:|---:|
| 2 | 4 | $-1$ | $86-7-76=3$ | 3 | 737 |
| 2 | 6 | $+1$ | $258-147-100=11$ | 11 | 711 |
| 4 | 4 | $+1$ | $8-3-4=1$ | 1 | 11 |

For the Cartan kernel over the rationals the upper bound is the cochain dimension minus the incoming rank and the rank of the joint map given by the outgoing differential and Cartan evaluation. Both modular ranks are lower bounds on the rational ranks. In each block this upper bound is zero. Hence every possible partner has a nonzero image in the rank-one abelian jet algebra.

Cartan evaluation intertwines the classical superconformal generators. Impose the two flavour raising operators, one dotted-spin operator, two conformal lowerings, three odd conformal lowerings, and six odd linear Taylor generators on that abelian algebra. The independently assembled 14 simultaneous equations have matrices with 4465, 4696, and 191 rows in the three cases. In the order of the table the square minors have determinants

$$
130460,\qquad 623010,\qquad 993091 \pmod{1000003}.
$$

All three have full column rank over the rationals. There are zero opposite-statistics lowest-weight multiplets of the threshold's commutant weight. In a unitary short spectrum each removed multiplet changes this weight's signed index multiplicity unless an opposite-statistics multiplet leaves with it. The superconformal index is independent of the exactly marginal coupling. Under the usual weak-coupling identification of the classical supercharge complex and the physical multiplet spectrum, the threshold therefore cannot acquire an instanton lift. This argument does not evaluate an instanton insertion. The result serves as a vanishing control for a projected Ward calculation.

With Python 3.12 or newer, python-flint and a C++17 compiler, run from this directory

~~~sh
python3 verify.py --out /tmp/su2-threshold-index.json
cmp /tmp/su2-threshold-index.json result.json
~~~

The replay builds the pinned full-jet engine, recomputes the five singlet columns and the lower representation at both primes, and regenerates the three simultaneous lowering matrices. [Raw ranks and determinants](result.json) · [Rebuild](verify.py) · [Integer Cartan maps](abelian.py)
