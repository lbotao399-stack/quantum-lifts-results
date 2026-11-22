# Two degree-17 relations and one surviving third-page source in SO(13)

Use the full relative Taylor-jet complex, the defining trace, and the ordered-pair one-loop kernel of the [public runtime](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7). The derivative-free highest flavour weight is $(7,5,5)$ and the leading source has sixteen letters. The first jet of flavour $i$ is $X_i$ in the code. The matter letter $P_{ij}$ has $Q_0 P_{ij}=[X_i,X_j]$ for $i<j$ in the runtime's ordered Taylor basis. Write $I_Y$ for the kernel of specialization of the stable orthogonal trace space at length $Y$, in this flavour representation. The symmetric-group character gives

$$
(\dim I_{13},\dim I_{14},\dim I_{15},\dim I_{16},\dim I_{17})=(0,9,35,69,42).
$$

## The classical quotient

The C++ program polarizes the orthogonal determinant on two auxiliary indices. The coefficient of Pfaffian degree greater than thirteen vanishes at thirteen colours. An open path containing more than one adjoint replaces several singletons by one matrix in this determinant. A closed trace multiplies the determinant. The scalar singleton couplings are traceless two-by-two auxiliary matrices. Coefficients of the stable expansion follow from one half the auxiliary and colour trace of the matrix logarithm. [Determinant construction and conventions](https://github.com/lbotao399-stack/quantum-lifts-results/blob/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/orthogonal_glueballs/README.md)

The probe bank antisymmetrizes five disjoint triples of flavour slots and fixes the other two to flavour zero. In the one-matter-letter bank its two slots are fused with the wedge sign. These tensors select the highest-weight multiplicity space. For a chosen cyclic trace word of length $l$, the connected coefficient on labelled slots is its collected trace-log coefficient times $2l$ divided by the size of its dihedral orbit. For a product of traces one further multiplies by the factorials of identical factors. These are the cyclic, reversal, and factor multiplicities in the multilinear pairing. The one-half trace logarithm on a circle can thus be computed by placing the tile containing site zero, then summing linear tilings of its complement with the alternating sign. A closed spectator is paired with all oriented cyclic cuts, so repeated letters and repeated factors have their full multiplicities.

An upper bound on the image uses twenty-nine independent closed identities at length sixteen. They come from the following three constructions.

1. Vary a determinant containing the three-flavour letter, or two independent odd matter letters, by $Q_0$. In the latter case the two differentiated polarizations carry opposite signs. These images are closed by $Q_0^2=0$.
2. On a pure-scalar determinant replace one scalar in every possible slot by $P_{vi}$, with $Q_0 P_{vi}=[X_v,X_i]$. The total tree variation is a simultaneous infinitesimal conjugation of an invariant.
3. Multiply a degree-fourteen scalar determinant by $\operatorname{Tr}(P_{uv}X_u)$. Its variation is $\operatorname{Tr}([X_u,X_v]X_u)=0$.

The sampled image of $Q_0:I_{16}\to I_{17}$ has rank forty. Each sampled column is the actual variation of a finite identity. The two product rows use seeds zero and one. Their leading relations are a scalar bilinear multiplying the degree-fifteen, all-singleton Pfaffian. The matrices and both bilinear factors are in the parameter file. Both determinants vanish at thirteen colours.

| Prime | 29 closed identities | 40 image rows | 40 image rows and both products |
|---|---:|---:|---:|
| 1000003 | 6477 | 326095 | 660854 |
| 1000033 | 790759 | 470068 | 174056 |

Each entry is the determinant of the indicated square minor. The complete row and column selections are saved. Rank-nullity gives $40\leq\operatorname{rank}_{\mathbb Q} Q_0\leq69-29=40$. The two products form a basis of $I_{17}/Q_0 I_{16}$ over the rationals. At odd terminal length seventeen every trace monomial contains an odd trace. Reversing this trace by successive interchanges gives its negative, and every interchange difference is a tree variation. Dividing by two produces a stable primitive of each product. The connecting map for the stable-to-finite trace quotient therefore realizes both products as finite classical cycles modulo stable cycles. Each may also be represented by the scalar bilinear times the primitive of its degree-fifteen relation. Stable cohomology in this matrix family consists of multi-gravitons and these admit closed quantum completions. [Orthogonal stable primitive](https://github.com/lbotao399-stack/quantum-lifts-results/blob/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118/orthogonal_odd_trace_primitive.py) · [Quantum product and multi-graviton theorem](https://arxiv.org/abs/2611.24160)

## First differential and the second page

The first target has dotted degrees $(1,1)$ and flavour degrees $(6,4,4)$ at length sixteen. On the two scalar product sources the ring derivation gives the scalar bilinear times the first variation of the primitive of the singleton Pfaffian. In the stable commuting quotient only one odd trace of the scalar determinant contributes. Terms with three or more odd traces retain an odd trace after the pair contraction. Contracting between the active trace and another trace gives zero in the commuting quotient since the trace of a generator times a commutator with either output gaugino vanishes there. Thus the first action is the single-odd-trace variation of the logarithm times its commuting exponential. The latter restricts on labelled planes to the product of $1+\frac12\operatorname{Tr}_2(A^2)$, where $A$ is the linear combination of the three singleton couplings on that plane. Seven- and eight-plane coordinates at this length require only active odd trace lengths three, five and seven. The implementation computes the tree primitive and the native $Q_1$ action for each and includes repeated-factor and exterior signs.

For the full flavour-highest and dotted-singlet first target, replace the two gauginos by one auxiliary scalar with two independent odd coefficients and extract their ordered product. This polarizes the symmetric colour pair in the dotted singlet. Signed-site orbit counting and the flavour and dotted Weyl differences give fifteen higher-plane coordinates modulo adjoint symmetry lowering. Fifteen independent multilinear contractions use four alternating flavour triples, two fixed highest-flavour slots and the two distinguished slots. They have support on 744 site monomials. The scalar prefactors include the factorials for repeated plane coordinates.

The finite determinant has at most two pairings left after removing the fourteen antisymmetrized endpoints on each side. Its commuting path types are singletons, one double path, two double paths, one triple path, and a closed bilinear times singletons. Substitution on an individual plane gives $1+q$, where $q=\operatorname{Tr}_2(A^2)/2$. An even double path enters the scalar part of the two-by-two auxiliary matrix. Thus the one-plane factor for two marked double paths $B,C$ is $(1-uB-vC)^2+q$. The mixed term on a single plane is $2uvBC$ and the term on separate planes is $4uvB C$. A triple path with auxiliary matrix $M$ contributes a multiple of $\operatorname{Tr}_2(AM)$ times its three scalar factors. The closed bilinear multiplies the product of plane factors. Fifteen parameter choices per type already supply fourteen independent relations. Appending the two first actions increases this rank to fifteen. A second computation on the six substituted components tests the same relation without presupposing that these fifteen samples exhaust the finite relations.

For that check, substituting all six distinct flavour-to-gaugino odd generators gives the intermediate weight $(4,4;4,2,2)$. The 76 weight-allowed determinant polarizations exhaust the high-plane image modulo internal tree commutators by the path reduction used for the eighth glueball power. Here the extra two scalar slots transform as the symmetric square. The three first-jet flavour axes are tested modulo the auxiliary maximal-torus weight; every zero-weight choice is included. The previous identity space at the preceding weight vanishes. The supplied second-dual equations were verified on 314 and 16,978 native columns. The rank of the 346-coordinate plane rows is five. Including both simple flavour lowerings has rank 320, of which 315 comes from lowering. This gives a functional on actual finite classes. On the two product sources its values are $218103808$ and $93972234240000$. A direct, one-pair calculation has respective evaluations $872415232$ and $375888936960000$. In particular the first map has rank exactly one. Its target has at most one excess direction after the fourteen finite identities, and the Cartan-visible directions cannot be first images by the discriminant theorem. [Eight-path method](https://github.com/lbotao399-stack/quantum-lifts-results/tree/f05b289fc69fb5f65bd8286c09060156a554f23d/reviews/orthogonal-mixed-lifts) · [First-order Cartan-image theorem](https://arxiv.org/abs/2610.22332)

Write $O_0^{(0)},O_0^{(1)}$ for the two product primitives after projecting to the highest component. Use the same projection on both. Dividing the evaluations by their greatest common divisor gives

$$
358475625\,[O_0^{(0)}]-832\,[O_0^{(1)}]\in\ker d_1 .
$$

The class is nonzero by the 42-row minor. There is no incoming positive-order map into its zero-derivative grading. Its second target has length fifteen and consists entirely of first jets. Every formal trace monomial has an odd trace and vanishes modulo the classical graded commutators. Hence the displayed class survives to $E_3$. The third target has dotted degrees $(3,3)$ and flavour $(4,2,2)$. A determination of its third differential requires the corrected third pairing.

## Cartan-visible targets at all perturbative orders

In the classical multi-trace sector where Kibe, van Leuven and Mathieson prove their discriminant identity, let $\Delta$ be the characteristic discriminant of one scalar and $E_1=H(Q_0)$. Their Eq. (110) says that for any class $x$ there is a power $k$ for which $\Delta^k x$ is a polynomial multi-graviton. Thus $E_1[\Delta^{-1}]$ is generated by multi-gravitons and $\Delta^{-1}$. Their Eq. (112) implies that a class with nonzero Cartan restriction has nonzero localization. [Classical input, §6.2](https://arxiv.org/html/2610.22332v3#S6.SS2)

The quantum ring theorem supplies closed lifts of all these graviton polynomials, including $\Delta$, and gives the Leibniz rule on each page. Assume inductively that $E_r[\Delta^{-1}]=E_1[\Delta^{-1}]$ via the preceding pages. The localized differential is a derivation, and vanishes on each generating graviton and on $\Delta^{-1}$ since $d_r\Delta=0$. It therefore vanishes identically. Localization is exact, hence

$$
E_{r+1}[\Delta^{-1}]
=H(E_r[\Delta^{-1}],d_r)=E_r[\Delta^{-1}].
$$

Induction proves that all perturbative differential images localize to zero. A target class with nonzero classical Cartan restriction cannot occur in such an image. At a later page this assertion concerns a representative that has survived the earlier pages, and restriction of earlier boundaries is zero by the same induction. A Cartan-visible source may have a nonzero unlocalized differential. The argument uses pagewise localization and requires no interchange of localization and formal completion. [Filtered quantum product and graviton closure](https://arxiv.org/abs/2611.24160)

## Replay

Use Python 3.13, NumPy, SciPy, python-flint and a C++17 compiler. From this directory set REPO to a checkout of the public trace runtime at commit ee9e06494a926fb88baeb6fd875cee36762f63a7. The rank verification compiles the C++ program in a temporary directory, regenerates every entry over both primes and re-evaluates the four character dimensions.

```sh
python3 verify_bmn_relations.py --repo "$REPO" --workers 2 --out replay-ranks.json
cmp bmn13-ranks.json replay-ranks.json
c++ -O3 -std=c++17 bmn_relation_probe.cpp -o /tmp/so13-bmn-probe
/tmp/so13-bmn-probe --product-parameters 0 8 > replay-parameters.dat
cmp bmn13-parameters.dat replay-parameters.dat
python3 mixed_bmn_d1.py --repo "$REPO" --relations mixed13-rows.json --parameters bmn13-parameters.dat --out mixed-replay.json --count 4
cmp bmn13-first-products.json mixed-replay.json
python3 bmn_one_pair.py --repo "$REPO" --parameters bmn13-parameters.dat --out paired-replay.json
cmp bmn13-one-pair.json paired-replay.json
```
