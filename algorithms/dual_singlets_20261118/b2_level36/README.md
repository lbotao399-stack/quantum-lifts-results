# The level-36 BMN restriction to $USp(4)$

The restriction of the published $SU(4)$ class $\mathcal O_{36}^{(4)}$ has a nonzero $SU(2)$ threshold component inside $USp(4)$. Its BMN charges are $(2J;q_1,q_2,q_3)=(5;10,4,4)$, its flavour representation is $[6,0]$ and its letter number is $Y=13$. The input from Kim, Lee, Lee and Song is the nonzero mixed coefficient in their equation (4.15). Their section 3 describes the finite-field checks used for large sectors. The calculation here reconstructs the ratio of two restriction coefficients over the rationals. The full parent trace polynomial has not been independently reconstructed. [Published parent and restriction](https://arxiv.org/html/2610.20340v2#S4.SS3) · [Integer calculation](wall_restriction.json)

## The restriction coefficient

Set the spectator fermions and field strengths to zero and polarize both spectator scalars along the highest flavour. In equations (2.15) and (4.14), write $h=\hat\phi^1$ and $t=\tilde\phi^1$. The scalar background is

$$
\operatorname{diag}(h+t,h+t,-2h+t,-3t).
$$

The first two colours carry the active $SU(2)$. Extract the coefficient of its seven-letter threshold class $\mathcal O_{24}^{(2)}$. The coefficient is homogeneous of degree six. The classical supercharge $Q_0$ acts by commutators, so restriction and coefficient extraction commute with $Q_0$. The threshold pairs to $-20$ with the saved rational functional, which annihilates the variations of all 197 spanning predecessor trace words. [Core and functional](../../su3_j34_20261118/threshold_identity.json) · [Replay](wall_restriction.py)

For outer scalar power $k=6,5,4$, the intermediate $SU(3)$ operator has level $36-2k$, length $13-k$, and charge-conjugation parity $(-1)^k$. Every graviton at that length has parity $(-1)^{13-k}$. The complete spin-five census at levels $24,26,28$ gives respectively $11=11$, $31=31$, $60=60$ cohomology and graviton multiplicities, summed over highest weights. All three opposite-parity coefficients vanish in cohomology. Thus the projected polynomial is divisible by $h^3$. At outer power three, the only excess is the $[3,0]$ threshold at level 30. The $[6,0]$ projection in the product of the two symmetric cubes has multiplicity one. Equation (4.15) gives a nonzero coefficient of $h^3t^3$. [Full $SU(3)$ columns](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su3/bmn_hw_cohomology_j30.json)

Exchanging colours three and four fixes the active doublet and takes $(h,t)$ to $((h+4t)/3,(2h-t)/3)$. Its invariance forces divisibility by $(h+4t)^3$ as well. If $c$ denotes the published mixed coefficient in the highest component, the full polynomial is

$$
\frac{c}{64}h^3(h+4t)^3.
$$

On the symplectic long-root block the spectator eigenvalues are $(z,-z)$ and the active offset vanishes. Hence $(h,t)=(-z/3,z/3)$ and the threshold coefficient is $-c z^6/1728$. The two constraints have rank six on the seven homogeneous coefficients. The integer direction, ordered by increasing power of $h$, is $(0,0,0,64,48,12,1)$. Odd-length gravitons vanish in $USp(4)$ since its invariant scalar polynomials have degrees two and four and their BMN descendants preserve length. [Exact constraint matrix](wall_restriction.json)

## Lower generators

In the rank-one cohomology through level 24 the sole odd-length highest weight is $\mathcal O_{24}^{(2)}$. In a product contributing to its wall coefficient, one factor must already give that class. Every remaining factor must give the identity times a spectator scalar polynomial. The spectator Weyl reflection requires even powers of $z$. The public $SO(5)=USp(4)$ ranks exhaust all 297 nonempty odd-length highest-weight blocks through level 28. Thus the only possible positive-degree dressing is a quadratic scalar multiplying an eleven-letter class at $(5;8,4,4)$, level 32. The highest weight is $[4,0]$ because four identical spectator scalars multiply the flavour singlet. Within the BMN sector, the three raising supercharges commute with extraction of the pure spectator scalar component. They annihilate a constant spectator background after its fermions and field strength are set to zero. A descendant contribution would require an odd-length active class below level 24. [Rank-one census](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/data/su2/bmn_hw_cohomology_j44.json) · [Rank-two predecessor audit](odd_through_j28.json.gz)

## Independent rank bounds

At $(10,4,4)$ the integer character gives 67,720 even-length and 67,698 odd-length cochains. Independent closed products contribute at least $3+19+2=24$ even-length classes. The total odd cohomology has dimension at least $24-(67720-67698)=2$. The nonzero minors in consecutive lengths have orders 768, 4,329 and 12,989. Consequently the individual classical dimensions at lengths 10, 11 and 12 are $771-768=3$, $5097-768-4329=0$ and $17337-4329-12989=19$. This lower bound of two covers the entire odd column. [Product determinants](rational_products.json) · [Length 10](length10.json.gz) · [Length 11](length11.json.gz) · [Length 12](length12.json.gz)

## Selection of higher-order targets

The index gives a generator-existence statement independent of the published $SU(4)$ restriction. The only previously counted $USp(4)$ core is the ten-letter level-30 singlet. The symplectic scalar invariants have degrees two and four. Products add letter lengths, and the BMN raising supersymmetries preserve length. These generators and their descendants therefore generate only even lengths. The odd cohomology lower bound of two forces at least one further BMN generator. The 297 exhausted odd-length blocks through level 28 bound the first such generator to $j\in\{30,32,34,36\}$. This argument does not determine its individual charges. [Rank-two census](odd_through_j28.json.gz) · [Exact product minor](rational_products.json)

For the published restricted class, 231 distinct scalar-insertion trial polynomials give classical-variation rank 230 modulo 1000003. Appending the long- and short-root evaluations, individually or together, still gives 230. This particular ansatz has no closed trial combination detected by either wall at that prime. [Trial polynomials and ranks](all_insertions.json.gz)

For the five Taylor slots take the two even spacetime directions and the three odd flavour directions. A BMN highest weight $(a,b,c)$ with $a\ge b\ge c\ge2$ occurs at zero spacetime degree in precisely the polynomial representation with column lengths $(a,b,c)$. With the even directions ordered first its highest weight is $(3,3;a-2,b-2,c-2)$. An order-$r$ supercharge preserves the special linear part and shifts the two even degrees by $r$ and the three odd degrees by $-r$. Under filtered superconformal equivariance its image must therefore have highest weight

$$
(3+r,3+r; a-2-r,b-2-r,c-2-r).
$$

Negative entries in the odd degrees exclude a polynomial target. At $(10,4,4)$ this removes $d_3$ and all higher outgoing orders under the stated Ward identity. At $d_2$ the Young diagram has rows $(5,5,1,1,1,1,1,1)$. Its twelve-letter multiplicity is 46,860. The original dotted-spin and flavour block has 833,176 cochains. At $(10,5,5)$ and $(11,5,5)$ the third-order diagrams exist. [Character calculation](schur_complex.json) · [Superconformal assumption and perturbative checks](https://github.com/davidsd/fortuity-ai/blob/4256f207986029540060fb44dd494746828f65ce/summary/fortuity-summary.tex)

The multiplicity calculation uses the exterior algebra of the odd superfield with positive symmetric Taylor jets. For each permutation of the letters it evaluates the trace on colour-invariant adjoint tensors by the $C_2$ Weyl formula. It then applies the symmetric-group character of the required diagram to the permutation of its Taylor slots. All 31 derivative-free multiplicities in three columns agree with the independent finite-colour characters. These are dimensions of cochains. The second obstruction requires the class of $Q_2\mathcal O_0+Q_1\mathcal O_1$ for the solutions of $Q_0\mathcal O_1=-Q_1\mathcal O_0$. [Algorithm](../b2_schur_complex.py)

From the repository root the two character and restriction calculations use the Python standard library.

```sh
d=algorithms/dual_singlets_20261118
python3 "$d/b2_level36/wall_restriction.py" > /tmp/b2-wall.json
cmp /tmp/b2-wall.json "$d/b2_level36/wall_restriction.json"
python3 "$d/b2_schur_complex.py" --out /tmp/b2-schur.json
cmp /tmp/b2-schur.json "$d/b2_level36/schur_complex.json"
```
