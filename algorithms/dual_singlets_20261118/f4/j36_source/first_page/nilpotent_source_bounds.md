# Classical source tests on singular (F_4) tuples

For the possible second supercharge differential into the sixth glueball at level 36, the ordinary source polarization has weight $(4,4,4,2)$ and contains one marked second jet. On the 30 Weyl orbits of quadruples of distinct commuting positive-root generators, every global classical source cycle is locally a classical boundary through scalar deformation degree six. At degree seven the bounds on the possible leading symbols are one and two for the two ordered cases below. The global corrected quantum pairing remains open. [Complete rational certificates](nilpotent_source_bounds.json.gz) · [Global pairing criterion](README.md#one-row-decides-the-second-differential)

### Local classical equations

Fix a commuting tuple $X_i$, with $i=0,1,2,3$ the four Taylor directions. The differential on adjoint-valued alternating forms is the alternating insertion of $[X_i,\cdot]$. For a scalar variation $u$ and marked second jet $Y$ its first two maps are

$$
(du)_{ij}=[X_i,u_j]-[X_j,u_i],\qquad
(dY)_{ijk}=[X_i,Y_{jk}]-[X_j,Y_{ik}]+[X_k,Y_{ij}] .
$$

The Jacobi identity and commutation of the $X_i$ give $d^2=0$. The adjoint invariant form pairs the degree-$k$ colour complex with degree $4-k$ at the complementary weight. Both identities are verified on all 30 orbits. The degree-one cohomology $H^1$ describes scalar variations modulo infinitesimal gauge transformations. The degree-two cohomology $H^2$ describes marks satisfying the linear Bianchi equation modulo linearized commuting equations. To take a transverse slice, choose a complement to the orbit tangent and use the group action on this complement. Its derivative, with the orbit parameters, is invertible at the chosen point. The formal change of coordinates gives the same leading complex. [Rational weight blocks and duality checks](nilpotent_tangent.py)

On this slice, the linear differential splits as contractible pairs and cohomology generators. At scalar deformation degree $k$, leading one-mark symbols therefore lie in $\operatorname{Sym}^k(H^1)^*\otimes(H^2)^*$. For root colours $X_i$ of weights $\alpha_i$, the weight of a deformation of colour $\beta$ on a set of Taylor axes is $\beta$ minus the sum of the corresponding $\alpha_i$. The weights in a nonzero evaluation must add to the negative of the weighted sum of the four root colours in $(4,4,4,2)$. The first such degree over the 120 short-axis choices is six. It occurs only on $(8,19,27,47)$ with colour 19 on the short axis. Its adjacent differential ranks are $41,123,123,41$, and its cohomology dimensions are $11,44,66,44,11$. Root numbers label the native basis in the archived colour file. [Root data](inputs.json.gz) · [Weight search](nilpotent_tangent.py)

The quadratic maps on this linear cohomology are the projections of the colour bracket with the exterior product of Taylor axes. The map out of a one-mark symbol replaces the mark by the quadratic scalar commutator. Incoming maps split two marks by the same commutator or split a third jet by the mixed scalar–mark bracket. The verifier computes both maps over the rationals and checks their product is zero. [Quadratic calculation](nilpotent_quadratic.py)

### An exact four-row test at degree six

A global source cochain vanishes when the scalar tuple commutes, for arbitrary marked colour, by the reduced-commuting-variety argument. Evaluate it on a polynomial family of commuting scalars with a marked two-form satisfying the Bianchi equation along the family. These evaluations also annihilate every classical boundary. Consequently they constrain the first non-exact homogeneous coefficient even after earlier boundaries are subtracted. A vanishing homogeneous class can be removed at that order by a formal primitive. [Global vanishing argument](regular_source.md#vanishing-on-the-reduced-commuting-variety) · [Polynomial tests](nilpotent_arcs.py)

In the five ordered coordinates in the archive the outgoing quadratic map vanishes and the incoming image is spanned by $(0,2,-2,0,0)$. The four exact evaluations give

$$
\begin{pmatrix}
0&0&0&0&1\\
0&0&0&1&0\\
1&0&0&0&0\\
0&1&1&0&0
\end{pmatrix}.
$$

Their joint kernel equals the incoming image. In the rational tangent basis the only possibly nonzero curvature on the six tangent directions used is twice the difference of the products of coordinates $(27,30)$ and $(11,40)$, in colour 49 on axes $(0,1)$. The four test families make these products equal. Direct substitution at $t=0,1,-1$ checks the commuting and marked equations of degree at most two. This proves the degree-six exclusion. A separate elimination of the cubic equation gives three nonzero individual columns and no obstruction after the 84 allowed first corrections. [Exact test families and pivots](nilpotent_source_bounds.json.gz) · [Cubic replay](nilpotent_cubic.py)

### Degree seven

Of 120 labelled cases, twelve have the required weight at degree seven. The table adds ranks across these twelve local tests. The entries in the last column bound possible leading coefficients in the separate completions.

| Scalar degree | Linear symbols | Incoming quadratic rank | Outgoing quadratic rank | After quadratic equations | After straight families | After polynomial families | After finite stabilizers |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 5 | 1 | 0 | 4 | 0 | 0 | 0 |
| 7 | 259 | 51 | 147 | 61 | 35 | 5 | 3 |

For a polynomial scalar family and mark of degree at most three, evaluate their original commutators and Bianchi equations at seven distinct integers. This proves the identities over the rationals. A finite Weyl stabilizer exchanging long Taylor axes acts trivially on the prescribed highest line under the fourth power of the permutation sign. The lifts, axis scalings and character of every stabilizer used are checked in the native Lie algebra. Finite-group averaging permits invariant representatives. The two remaining bounds are $1$ on $(8,19,25,46)$ and $2$ on $(8,19,27,47)$, both with colour 19 on the short axis. Attainment by global cycles and their corrected second pairings remain open. [Polynomial curves](nilpotent_curved_arcs.py) · [Stabilizer actions](nilpotent_stabilizer.py) · [Complete test](nilpotent_source_bounds.py)

### Two proper Levi subgroups

For the ordinary diagram $(4,4,2,2)$ the physical polarization $(4,4,0;2,2)$ occurs with multiplicity one. The integer cochain dimensions come from the group and Taylor characters. At both primes 1000003 and 1000033, the incoming and outgoing classical minors attain these dimensions together. Nilpotence gives the corresponding upper bound over the rationals.

| Gauge group | Parent dimension | Source dimension | Terminal dimension | Incoming rank | Outgoing rank | Rational source cohomology |
|---|---:|---:|---:|---:|---:|---:|
| $SO(7)$ | 16347 | 7266 | 1435 | 5837 | 1429 | 0 |
| $USp(6)$ | 16227 | 7206 | 1422 | 5790 | 1416 | 0 |

The four A8 jobs were 5706542, 5707146, 5707154 and 5707188. For the remaining proper rank-three Levi diagram $(4,4,4,2)$ the corresponding source dimensions are 18894 and 18425. The subgroup calculations constrain the support of a global source. [Scheduler receipts and pivot indices](levi_alternate_minors.json.gz) · [Verifier](levi_alternate_minors.py) · [Character reconstruction](levi_relation_spaces.json)

From this directory the first two commands replay the local bounds and inspect the scheduler-backed minors. The last two regenerate any selected native minor on an available compute node.

~~~sh
python3 nilpotent_source_bounds.py
python3 levi_alternate_minors.py
c++ -O3 -std=c++17 -pthread -I ../../../taylor_hooks/native ../../../taylor_hooks/su2_cartan_kernel.cpp -o /tmp/f4-levi-minor
python3 levi_alternate_minors.py --engine /tmp/f4-levi-minor --case SO7-44022-1000003.json --threads 8
~~~
