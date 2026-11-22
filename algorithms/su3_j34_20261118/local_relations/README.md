# Formal perturbative survival of the bosonic level-34 SU(3) seed

The charge-conjugation-even representative of Kim, Lee, Lee and Song's $\mathcal O_{34}^{(3)}$ admits a closed lift over $\mathbb Q[[u]]$ in the full relative holomorphic Taylor-jet complex of undeformed $\mathcal N=4$ SYM. Here $u$ counts supercharge corrections in the Li–Liu–Chang convention. The resulting class is free of $u$-torsion. Its physical instanton differential is open. We adapt Dot's level-30 local-relation argument to the different relation $R_{20}^{3}[1,0]$. [Published seed and Appendix A](https://arxiv.org/html/2610.20340v2#S4.SS1) · [Dot's argument](https://github.com/lbotao399-stack/quantum-lifts-results/blob/2e22109bc140b10e6a9bf70892dcf6b586510475/nonperturbative-lift/contributions/dot-r007-all-order/ALL_ORDER_PROOF.md)

## The classical identities

Only $R_{14}^{2}[1,3]$, $R_{16}^{2}[2,3]$ and $R_{20}^{3}[1,0]$ occur. For their highest components, expanding $Q_0r-R$ on traceless three-colour matrices gives zero on 76, 102 and 17992 distinct monomials, respectively. The three $R$ controls have 32, 44 and 5484 nonzero coefficients on the same specializations. One scalar is diagonal, the other matrices are generic, and invariance extends each identity from the dense regular-semisimple locus. Flavour covariance gives the other components. [Finite-colour replay](verify_local.py)

Expand the paper's 23 contractions in independent graviton symbols, with just their linear symmetries and trace subtractions. Their 5369 separate monomials cancel exactly. Replacing $R$ by its displayed primitive $r$ and applying the sign for the fermionic prefactor gives the archived 5815-term polynomial. Its even projection agrees coefficient by coefficient with all 6270 terms. Project each local primitive to the parity of its $R$ to obtain this even source. [Contraction and preimage replay](syzygy.py) · [Archived representative and projection](../README.md)

Replace graviton products by closed marked products in one common normalization. A change of bracketing or a graded interchange of closed factors is a supercharge boundary with primitive starting at order $u$. The graded free cancellation therefore gives an exact relation between the marked trees up to that boundary. These steps use the marked Leibniz identity, homotopy and graviton theorem in Li–Liu–Chang, equations (6.13), (6.26), (8.3) and Theorem 9.1. [Marked identities](https://arxiv.org/pdf/2611.24160v1)

## The three local corrections

Order $k$ removes $k$ from every flavour degree and inserts dotted degrees $(k,k)$. The highest flavour degrees of the first two relations are $(4,3,0)$ and $(5,3,0)$. Only $R_{20}^{3}$, with degrees $(4,3,3)$, can have a local quantum defect. This relation has colour parity $-1$. In the table, the last two columns give the rank of the outgoing $Q_0$ and the additional rank of Cartan restriction on its evaluated kernel, measured on common cochain columns. Each is used through the lower bound on their stacked rank.

| Order | Flavour | Length | Exact dimension | Incoming rank bound | Outgoing rank | Additional Cartan rank |
|---:|---|---:|---:|---:|---:|---:|
| 1 | $(3,2,2)$ | 6 | 761 | 312 | 426 | 23 |
| 2 | $(2,1,1)$ | 5 | 430 | 146 | 270 | 14 |
| 3 | $(1,0,0)$ | 4 | 30 | 7 | 21 | 2 |

Both primes give these lower bounds on the rational maps. Independent integer characters give the three dimensions. Every incoming image has zero outgoing differential and zero restriction. The three sums saturate, so every classical cocycle with zero restriction in each row is a boundary. At even length, charge conjugation acts as $+1$ on the Cartan. The first defect is consequently a boundary, and the same argument applies to the third after the second has been repaired. Choose the primitives with the required flavour representation and colour parity by equivariance. [Characters and rebuild at both primes](ranks.py)

For the second defect, choose local primitives through order one and subtract the marked reassociation homotopy from the global source. Its order-two residual is the sum of the three contractions numbered 37, 38 and 40, with the components of the local defect in place of $R_{20}^{3}$. Its Cartan restriction vanishes by the all-page Cartan-image theorem applied to this partial global lift. That theorem follows from the classical discriminant lemma, the marked graviton lifts and the page Leibniz rule. [Cartan-image proof](https://github.com/lbotao399-stack/quantum-lifts-results/blob/2e22109bc140b10e6a9bf70892dcf6b586510475/nonperturbative-lift/contributions/dot-cartan-image/PROOF.md)

Use the paper's flavour convention here with highest flavour $1$ and lowest flavour $3$. Denote the Cartan restrictions of the three upper-flavour components of the order-two local defect by $e^i$. The Cartan multiplication tested is

$$
199920\,v^{1}{}_{i}v^{i1}{}_{3} e^1
-97104\,v^{1}{}_{i}v^{11}{}_{3}e^i
+237048\,v^{1}{}_{3}v^{11}{}_{i}e^i .
$$

The second graviton in each product is the three-letter graviton. The even cyclic flavour permutation sends $e^1$ to $e^2$ and $e^3$ with the exterior signs on the Taylor letters. The 4678 raw degree-five monomials on two independent Cartan colours split into 20 blocks by the multiset of dotted and flavour degrees of each letter. The two simple Weyl invariance equations and the three highest-weight equations have a joint kernel of dimension 53 over $\mathbb Q$. Independently, the identity, transposition and three-cycle have highest-weight traces $356,0,-19$, giving $(356+3(0)+2(-19))/6=53$. Each multiplication adds three undifferentiated scalars and two undifferentiated matter fermions, so outputs from different blocks have disjoint supports. Integer kernel and rank calculations give image rank 53. The three prefactors were also reconstructed directly from the paper-field traces and have 36, 14 and 4 terms. Thus $e^1=0$ and all three components vanish. The middle row of the table supplies the missing primitive. [Exact multiplication and symmetries](cartan_multiplier.py) · [Paper-trace comparison](verify_local.py)

After the third correction no higher local defect has nonnegative flavour degrees. Substituting these corrected local primitives in the marked syzygy gives a closed global lift. In each allowed finite charge and loop block the normalized differential has rational coefficients. Solvability over an extension gives solvability over $\mathbb Q$ by linear elimination. Finally the previously constructed mixed $SU(2)$ block functional on the zero-jet projection evaluates to $120294720$ on the source and kills all classical boundaries. Higher supercharge coefficients increase dotted degree, so the coefficientwise functional kills every formal boundary. Its nonzero constant pairing excludes $u$-torsion. [Exact block functional and rational primitive](../match_threshold.py)

## The independent level-36 seed and the bosonic tower

The repaired marked $R_{20}^{3}$ is a local tensor identity and can be used again with the same marked trees. The level-36 syzygy uses three additional relations of degree 12. For $[0,3]$ and $[2,2]$ the smallest highest flavour is zero. For $[1,1]$ the only possible defect has dotted degrees $(1,1)$, flavour $(2,1,0)$, length four and $C=-1$. Its Cartan restriction vanishes by parity. The exact cochain dimension and joint bounds give $12-4-4-4=0$. Thus all six local relations admit compatible repairs. [Published contractions](level36_terms.json) · [Finite-colour identities](verify_local.py)

All 15259 separate monomials in the 52-term level-36 graviton identity cancel freely. Substitution of the six classical preimages gives 11686 terms. The $C=-1$ projection has 12106. On the mixed coefficient of four abelian scalars of flavour 1 and one abelian fermion of flavour 3 the threshold functional gives $19780588800$. Its value on the $SU(2)$ threshold is $-20$. The ratio is $-989029440$, and multiplication by the flavour dimension 35 gives the paper's singlet coefficient $-34616030400$. Subtracting this multiple of the threshold gives exactly $Q_0$ of a 30-term primitive. The marked syzygy therefore gives a non-torsion formal class at $j=36$. [Exact source, pairing and reduction](level36.py)

For either seed, let $n=3$ or $4$ be its abelian scalar degree. The highest quadratic scalar graviton restricts to the sum of its $SU(2)$ trace and $6(\hat\phi^1)^2$. In each seed the coefficients of $(\hat\phi^1)^{n+2}\hat\psi_3$ and $(\hat\phi^1)^{n+4}\hat\psi_3$ both vanish identically on the $SU(2)$ block. Higher relevant powers are forbidden by flavour. After multiplication by $k$ quadratic gravitons the mixed coefficient at degree $n+2k$ therefore pairs to precisely $6^k$ times the original nonzero answer. Its leading term is a highest component in $[n+2k,1]$. Fixed-bracketing marked products give a closed formal representative. This supplies at least one non-torsion class at every $j=28+2n$ in the published bosonic $B_1$ tower, $n\ge3$. Its instanton differential still requires a physical matrix element. [Higher-power checks](level36.py) · [Original classical tower](https://arxiv.org/html/2610.20340v2#S4.SS1)

## Grade-zero instanton partners at level 34

At zero dotted degree, the even bosonic seed has eleven letters. Opposite-statistics partners of the same charges have ten letters on the incoming side and twelve on the outgoing side. The exact charge-conjugation-even cohomology dimensions there are 10 and 6. All sixteen directions have independent Cartan-visible graviton products. In this grade there is no additional even fortuitous partner. A physical instanton map into or out of a formally protected product requires separate physical information.

The public finite-field ranks give upper bounds on the rational total cohomology dimensions at lengths $9,10,11,12,13$ of $1,11,18,6,1$. The remaining lengths have upper bound zero. In the same order, exact rational Cartan minors exhibit $1,10,16,6,1$ graviton products, all with $C=(-1)^Y$ at letter length $Y$. The parity-resolved cochain characters give Euler characteristics $14$ and $-17$ using the length sign. There is one assignment consistent with these bounds and Euler characteristics.

| Length | 9 | 10 | 11 | 12 | 13 |
|---|---:|---:|---:|---:|---:|
| $C=+1$ | 0 | 10 | 2 | 6 | 0 |
| $C=-1$ | 1 | 1 | 16 | 0 | 1 |

The verifier checks every product determinant in both cyclic and ordinary traceless colour coordinates, recomputes the integer cochain dimensions, and solves the two Euler equations. This count uses no filtered Taylor Ward condition on an instanton. Positive dotted grades and the physical matrix elements are open. [Partner calculation](bmn_partners.py) · [Public complete column](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/hw_cohomology_8_5_4.json)

## Replay

Use Python 3.13 with SymPy 1.14, python-flint, NumPy and SciPy, and a C++17 compiler. Extract the frozen R-007 packet and use its directory containing `fortuity/bmn/dmk.py` as the vendor argument. Its two source hashes are verified before import. From this directory run

~~~sh
python3 syzygy.py /path/to/q0-certificate/vendor free
python3 syzygy.py /path/to/q0-certificate/vendor source
python3 verify_local.py /path/to/q0-certificate/vendor
python3 cartan_multiplier.py
python3 ranks.py
python3 ../match_threshold.py
python3 level36.py /path/to/q0-certificate/vendor free
python3 level36.py /path/to/q0-certificate/vendor source
python3 bmn_partners.py
~~~

[Frozen public vendor](https://github.com/lbotao399-stack/quantum-lifts-results/tree/60282a6b37e16fddf80990d6d7eb97a908082067/nonperturbative-lift/contributions/dot-r007-exact-input)
