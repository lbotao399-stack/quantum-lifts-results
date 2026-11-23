# The lifted-colour functional for the fifth glueball

For the charge-one instanton embedded along the long root of $F_4$, the fifth glueball at level 30 has two terms on linear radial gauginos. In the previously fixed normalization, their fourth-flavour doublet factors are $\omega^5$ and $5S_{\rm trip}\omega^4$, multiplied by the fifth and fourth powers of $\det M(x)$ respectively. Here $S_{\rm trip}$ uses two exact fermion modes and $M(x)$ is the radial doublet profile. At the minimum allowed source-plus-current shift in either term, every remaining fourth-flavour insertion outside the quartic vertices is an exact mode. The all-doublet term then requires four quartic vertices, and the triplet term requires six. [Mode count and normalization](README.md#integrating-one-flavour-of-lifted-instanton-modes)

Write the archived native quartic as a term linear in each fourth-flavour variable. Its coefficient on each such variable, placed at the right, is a cubic Grassmann polynomial in the first three flavours. For one of the seven ordered doublet pairs with form coefficient $2b$, denote its two cubic coefficients here by $L_i,L_j$. Integration of this pair, in its indicated order, gives $2b-L_iL_j$. Thus the sextic from this pair is $-L_iL_j/(2b)$. Adding the seven sextics gives the polynomial whose complete sparse coefficients are archived. Multiplication of the seven pair integrals, with permutation to the native colour order, gives the Pfaffian $-1024$. If the quartic is taken to divided power $2r$ and the fourth-flavour form to power $7-r$, the integral before the remaining three flavours are integrated is

$$
(7-r)!(-1024)\,\frac{(\text{archived sextic})^r}{r!}.
$$

At $r=2$ this gives $-61440$ times the squared sextic. At $r=3$ it gives $-4096$ times its cube, or $-20480$ after the radial factor five. These coefficients use the native quartic polarization of $B(f,({\rm ad}\,v)^4 f)$ with $B=\operatorname{Tr}_{26}/6$. Physical matching of the coefficient of that quartic and the source–current insertion remains open. [Native quartic](fifth_quartic.json) · [Sparse sextic, stored with denominator 16](fifth_radial_sextic.json.gz)

Each of the first three flavours contributes an alternating pair of the fourteen lifted colours. The sextic has 9058 nonzero coordinates. Flattening any one pair against the other two gives rational rank $\binom{14}{2}=91$. A selected integer determinant for the stored scaled tensor proves the rank. Direct exterior integration of the original four-flavour quartic checks twelve coordinates at two, four and six vertices. At six vertices the four checked values are $-552960,-7077888,-1013760,-172032$. The last value uses an off-diagonal colour pattern. The rank and these numbers describe the colour functional at fixed bosonic moduli. The physical nonperturbative lift has no certified matrix element. [Exact minor and test inputs](fifth_radial_kernel.json) · [Independent exterior and factorized calculations](fifth_radial_kernel.py)

~~~sh
python3 fifth_radial_kernel.py --out /tmp/f4-radial-kernel.json
cmp fifth_radial_kernel.json /tmp/f4-radial-kernel.json
~~~
