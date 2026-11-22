# Flavour support of the charge-one instanton in the SU(3) level-30 Ward channel

At the conformal origin the classical instanton profile of the outgoing class of R-007 vanishes at a single insertion in either orientation. This holds for every local operator with its R charges, including all perturbative mixtures with those charges. The bare R-007 source also vanishes separately. The conclusions concern evaluation on the classical one-instanton supermoduli. A physical lift requires the projected quantum Ward matrix. [Pinned source and target](README.md) · [Green and Kovacs, §§3.1–3.3 and Appendix C](https://arxiv.org/pdf/hep-th/0212332) · [Alday and Korchemsky, §2](https://arxiv.org/pdf/1605.06346)

## A weight test without the operator

For a charge-one $SU(N)$ instanton the fermionic coordinates of each $SU(4)_R$ flavour consist of four exact modes and $2(N-2)$ further modes. All carry the same fundamental flavour index, including both additional modes usually denoted by $\nu,\bar\nu$. Use the physical charges $R_I$ of the jet dictionary. In the orientation with chiral matter zero modes, the three matter-fermion weights are $(1/2,1/2,1/2)-e_I$ and the fourth weight is $(-1/2,-1/2,-1/2)$. A monomial with flavour occupations $n_1,n_2,n_3,n_4$ therefore satisfies

$$
R_i-R_j=n_j-n_i,\qquad R_i+R_j=n_k-n_4
\qquad (\{i,j,k\}=\{1,2,3\}).
$$

Every $n_A$ lies between zero and $2N$. The opposite orientation reverses the weights. Thus a necessary condition on the classical profile at a single insertion is

$$
\max_{1\le i<j\le3} (|R_i|+|R_j|)\le 2N .
$$

The level-30 outgoing operator has charges $(5,2,2)$, so $5+2=7>6$. Its radial conjugate has the same absolute charges. Equivalently, its four occupations in the two orientations must be common integer translates of $(5,2,2,9)$ and $(4,7,7,0)$. Each tuple has spread seven. The collective-coordinate exterior algebra contains neither weight. [Physical charge dictionary](../../../supplement/supplement.tex) · [Outgoing charges and representative](https://github.com/lbotao399-stack/quantum-lifts-results/blob/8ba864d4eb9413131922b75aa653ab119d10a04a/nonperturbative-lift/contributions/dot-r007-t2-primary/PROOF.md)

## The bare source at one point

Every one of the 6,139 source terms contains five undifferentiated scalars and five chiral matter fermions. In the chiral-zero-mode orientation a scalar of flavour $i$ takes the two complementary matter flavours at leading order, and a matter fermion of flavour $i$ takes flavour $i$. Directly counting the actual words gives $(3,6,6,0)$ on every word. A higher term with the same Cartan weight adds the same integer to every entry. Flavours two and three saturate six and flavour four already has zero. Only total degree 15 is possible, so every contributing elementary profile would have to be its lowest component.

At this order the ADHM equations (C.3) and (C.5) involve the projected matrices $\bar U\mathcal M^A$ and $\bar{\mathcal M}^A U$. Here $U$ spans the null space of $\bar\Delta(x)$, the one-instanton ADHM matrix. In (C.14) one can subtract $4\Delta(x)\bar\xi^A$ from $\mathcal M^A$ and the corresponding term from its barred projection. Both projections thus depend on the two components of $\eta^A-(x-x_0)\bar\xi^A$, and on $\nu^A,\bar\nu^A$. The other scalar term in (C.5) contains the antisymmetric flavour part of $\bar{\mathcal M}^A\mathcal M^B$. The exact spinor contractions are symmetric in $A,B$. The mixed exact and lifted terms vanish by $\bar w\nu=\bar\nu w=0$, with $w$ the bosonic orientation matrix. Hence this extra term uses only lifted modes. At $N=3$ the lowest profiles at one point depend on at most $2+2(N-2)=4$ independent odd variables of each flavour. Six factors of flavour two already vanish. This argument only uses the lowest scalar and chiral gaugino fields. [ADHM constraints and projected fields, equations (3.3), (C.3)–(C.8) and (C.13)–(C.17)](https://arxiv.org/pdf/hep-th/0212332)

In the opposite orientation the lowest component of each source letter occupies the complementary flavours out of four. The source minimum is $(7,4,4,10)$, which exceeds six in flavours one and four. The higher elementary terms add multiples of the four-flavour determinant. Thus the bare source has zero profile in this orientation as well.

## Free contractions between the two endpoints

In the separated three-point calculation take a classical supercurrent and allow Gaussian propagators only between the bare source and the radial target, in the bosonic instanton background. The minimum current occupations required by its R weight are $(1,1,1,0)$ and $(0,0,0,1)$. Adding source, target and current gives $(9,9,9,9)$ and $(11,11,11,11)$. A direct scalar propagator or conjugate fermion propagator removes one of each flavour from these lower bounds. A covariant derivative on the radial field preserves the bound. Consequently the two orientations require at least **three** and **five** such direct propagators respectively. Current contractions, higher interaction vertices, terms in the matched operators and contact terms require their own calculation. The lower bounds supply no value for an allowed coefficient.

Run the following from this directory in a checkout containing the pinned input commit. The verifier reconstructs both endpoints, checks all 6,139 and 5,080 terms, and tests both sets of integer Cartan occupations.

~~~sh
python3 flavour_modes.py > /tmp/r007-flavour.json
cmp /tmp/r007-flavour.json flavour_modes.json
~~~
