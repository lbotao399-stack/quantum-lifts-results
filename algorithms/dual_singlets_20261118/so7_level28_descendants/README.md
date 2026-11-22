# First-supercharge images of the dressed level-26 SO(7) multiplet

In the full relative holomorphic jet complex of undeformed four-dimensional $\mathcal N=4$ SYM, the two ten-letter $SO(7)$ BMN classes outside the multi-graviton span at $j=28$ both have nonzero first supercharge obstruction. All other BMN classes at that level have graviton representatives. Consequently the outgoing second differential from the entire derivative-free $j=28$ source sector has rank zero. The full descendant calculation also gives two first-image flavour multiplets at $j=30$ and one at $j=32$. The finite-colour supercharge uses the normalization of the [level-26 parent](../so7_level26/README.md). Fortuity reported the two non-graviton classes at level 28 and their possible descendant charges. [Fortuity census, §5](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex#L593)

## Return to the known first image

Let $F$ be the archived ten-letter highest-weight source of flavour $(5,4,4)$. Its class is nonzero and its first variation $Q_1F$ increases the rank of the classical boundary matrix from 12929 to 12930. The flavour multiplet of this image is a triplet. Let $D_i$ be the outer odd Taylor derivative in flavour direction $i$ and $S_i$ the classical superconformal lowering operator on the divided Taylor jets. Their action in the native zero-based flavour basis is the one in Dot's [jet replay](https://github.com/lbotao399-stack/quantum-lifts-results/blob/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-t2-primary/replay.py#L73).

On a letter of total Taylor degree $d$ with alternating flavour set $I$, $D_i$ appends an absent $i$ with the sign for placing it in sorted order. For $i\in I$, $S_i$ removes it with that same reordering sign multiplied by $-(d-1)$. The other cases give zero. On a trace polynomial of total Taylor degree $n$ they satisfy

$$
\{S_i,D_j\}=E_{ji}-\delta_{ij}n,\qquad \{D_i,D_j\}=0 .
$$

Here $E_{ji}$ replaces a flavour $i$ by $j$ in every slot, with $E_{ii}$ counting occurrences of $i$. Both odd actions use the graded product rule. On the jet letters with at most one dotted derivative of each sign, $\{Q_0,S_i\}=0$ and $\{Q_0,D_i\}$ is the commutator with the scalar of flavour $i$, with the scalar on the right in the first term. The latter vanishes on invariant trace polynomials. These are all the letters in a possible classical primitive of any returned first variation, since a first variation of a BMN state has dotted degrees $(1,1)$. On every ordered pair of BMN input letters the rational first kernel obeys $\{Q_1,D_i\}=\{Q_1,S_i\}=0$. The exact checks give 186 classical identities and $7^2\cdot3=147$ ordered-pair identities for each odd operation. The pair test precedes colour contraction, so it applies to every trace splitting and joining. Spectator terms cancel by the odd product rule. Successive raising operations preserve the BMN input condition. [Rational verifier](verify.py) · [Native first kernel](../../su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-protected/q1d_exact.py)

The three lowerings of the triplet take values in the singlet and octet at level 24 and length ten. At highest weights $(4,4,4)$ and $(5,4,3)$ the respective classical dimension bounds are

$$
502-313-189=0,\qquad 2432-1546-873=13 .
$$

Integer colour characters fix the cochain dimensions. Two native finite-field computations give the adjacent rank lower bounds. Thirteen explicit highest-weight products of closed single-graviton factors have Cartan determinant $-7409643204723800691383785960243200$. The defining traces of positive-degree even powers use the three opposite pairs of Cartan eigenvalues in $SO(7)$; the seventh eigenvalue is zero. The verifier reconstructs the nonabelian closed factors, the highest-weight identities and the integer determinant. Thus all lowerings of $F$ have closed perturbative representatives, and every $S_i$ kills the classical class of $Q_1F$ and of its flavour partners. [Singlet](native-444.json) · [Octet](native-543.json) · [Graviton recipes](gravitons.json.gz)

Write $F_1=E_{10}F$, so that $E_{01}F_1=F$. The two highest vectors $D_0F$ and $D_1F-D_0F_1$ have respectively 5111 and 10110 rational trace terms. Lower them with $S_0$ and $S_1$. After taking the first variation the terms containing $D_iS_jF$ or $D_iS_jF_1$ are classical boundaries. With brackets denoting classical cohomology the signed identities are

$$
[S_0 Q_1 D_0F]=8[Q_1F],\qquad
[S_1 Q_1(D_1F-D_0F_1)]=10[Q_1F].
$$

The opposite signs $-8,-10$ apply when $D$ acts on $Q_1F$ instead. Either displayed equality detects a nonzero class since each $S_i$ takes classical boundaries to boundaries.

## All exterior degrees

On the image triplet $Q_1F$ the three flavour weights are the permutations of $(4,3,3)$ and $n=12$. Apply a subset of the three raising operators in increasing index order in the written product. For every subset $i_1<\cdots<i_p$ of the same size, apply $S_{i_1}$ first and then the remaining lowerings in order. Project the result onto each of the original three components. The return matrix acts on $\wedge^p\mathbf3\otimes\mathbf3$. Commuting each $S_i$ past the first $D_j$ reduces to the displayed anticommutator and a word with one fewer raising. In the resulting exterior word $E_{ji}$ acts on every exterior slot and on the last triplet, with $3\delta_{ij}$ from the determinant weight. The total Taylor degree there is twelve plus the number of remaining raising operations.

| Number of raising operations | Level | Flavour Dynkin labels | Return eigenvalues with multiplicity |
|---:|---:|---|---|
| 0 | 26 | $[1,0]$ | $1$ (3) |
| 1 | 28 | $[2,0],[0,1]$ | $-8$ (6), $-10$ (3) |
| 2 | 30 | $[1,1],[0,0]$ | $80$ (8), $110$ (1) |
| 3 | 32 | $[1,0]$ | $-880$ (3) |

Every return matrix is invertible, so the 24 first-image components are linearly independent within their indicated grades. Anticommuting $D$ through $Q_1$ multiplies each degree-$p$ matrix by $(-1)^p$ for the first variations of the specified descendant sources. The five positive-degree highest-weight multiplets are new to the campaign ledger. Each has a ten-letter BMN source; the previously counted $j=30$ singlet into $SF$ has dotted source degrees $(3,3)$. The three BMN raising supercharges generate all these components from the original flavour triplet. [All four integer matrices by recursion](verify.py) · [Recorded ranks and characteristic polynomials](result.json)

## Exhaustion at level 28

For each dominant BMN flavour partition of fourteen, compare the complete integer cochain dimensions with the modular ranks recorded by Fortuity. On a given middle space, their code calculates the rank of incoming boundaries and graviton columns together; adjoining the outgoing rank gives a rational lower bound for the number of directions accounted for. The 145 nonempty highest-weight blocks agree with the independent character calculation. In 143 blocks the combined ranks saturate the cochain dimension. The two remaining upper bounds on the quotient by gravitons are

$$
13755-5412-37-8305=1,\qquad
10217-4057-15-6144=1.
$$

The two displayed first variations attain these bounds because graviton products admit closed quantum representatives. Thus every BMN leading class at this level surviving $d_1$ has a closed graviton representative. This proves the second-order source exclusion without any quantum Taylor selection assumption. [Upstream computation of the augmented rank](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/code/fortuity/hwcoh/sectorpar.py#L243) · [Pinned public columns](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/g2/so7_hw_cohomology_j28.jsonl) · [Perturbative product theorem](https://arxiv.org/abs/2611.24160)

Run in this directory with Python 3.12 or later, NumPy, SciPy, python-flint and a C++17 compiler.

~~~sh
c++ -std=c++17 -O3 -pthread ../taylor_hooks/native/engine.cpp -o /tmp/so7-jet-engine
python3 verify.py /tmp/so7-jet-engine | diff -u result.json -
PAR_THREADS=4 A_CHECK=1 /tmp/so7-jet-engine chain SO7 4 4 4 0 0 10 10 772801 1000003
PAR_THREADS=4 A_CHECK=1 /tmp/so7-jet-engine chain SO7 5 4 3 0 0 10 10 772801 1000003
~~~
