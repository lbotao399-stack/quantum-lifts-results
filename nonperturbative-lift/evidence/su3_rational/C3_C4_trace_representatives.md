# Exact trace representatives for SU(3) $C_3$ and $C_4$

These formulas preserve the saved representatives, source trace-dual $b$, Chevalley $c$, and the original loop parameter. The tables contain trace monomials, not color-component expansions.

The source rotation generators $J_+,J_-$ act by the Leibniz rule. On a jet with nonnegative derivative orders $m,n$, their action is

$$
J_+(\partial_1^m\partial_2^n b)=n\partial_1^{m+1}\partial_2^{n-1}b,\qquad J_-(\partial_1^m\partial_2^n b)=m\partial_1^{m-1}\partial_2^{n+1}b,
$$

and identically on $c$. A zero prefactor means the corresponding term is zero. For each table, multiply each trace row by its coefficient, add the rows, then apply the displayed polynomial in $J_-J_+$. The factors remove every positive integer rotation spin present in that weight-zero input. This is the original `scalar` projection written explicitly.

## $C_3$

Here $[\hbar^j]$ extracts the coefficient of the indicated power of the unchanged loop parameter. The source expressions are $C_0=\operatorname{tr}(\partial_1c\partial_2c)$ and $a_1=-\frac12\operatorname{tr}(b\partial_1c\partial_2c-b\partial_2c\partial_1c)$.

$$
C_3=[\hbar^0]C_3+\hbar[\hbar]C_3+\hbar^2[\hbar^2]C_3.
$$

$$
[\hbar^0]C_3=-\frac{111}{112}\operatorname{tr}(b^2)a_1+\frac{57}{56}\operatorname{tr}(b^3)C_0
+\frac{15}{14}\left[\operatorname{tr}(b^2\partial_2c)\operatorname{tr}(b\partial_1c)-\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b\partial_2c)\right].
$$


### Coefficient of $\hbar^1$ in $C_3$

Apply

$$
\left(1-\frac{J_-J_+}{2}\right)\left(1-\frac{J_-J_+}{6}\right)
$$

to the coefficient-weighted sum of these 11 trace monomials, preserving the displayed factor order.

| Coefficient | Trace monomial |
|---:|:---|
| $-\frac{9}{2}$ | $\operatorname{tr}\!\left(b^{2} (\partial_2^{2}c) (\partial_1^{2}c)\right)$ |
| $-\frac{81}{16}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_2c) (\partial_1^{2}c)\right)$ |
| $\frac{15}{8}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_1^{2}c) (\partial_2c)\right)$ |
| $\frac{15}{8}$ | $\operatorname{tr}\!\left(b (\partial_2c) (\partial_1^{2}c) (\partial_2b)\right)$ |
| $\frac{81}{16}$ | $\operatorname{tr}\!\left(b (\partial_2c) (\partial_2b) (\partial_1^{2}c)\right)$ |
| $\frac{33}{16}$ | $\operatorname{tr}\!\left(b (\partial_1^{2}c) (\partial_2c) (\partial_2b)\right)$ |
| $-\frac{33}{16}$ | $\operatorname{tr}\!\left(b (\partial_1^{2}c) (\partial_2b) (\partial_2c)\right)$ |
| $\frac{57}{8}$ | $\operatorname{tr}\!\left(b (\partial_2^{2}b) (\partial_1c) (\partial_1c)\right)$ |
| $\frac{33}{4}$ | $\operatorname{tr}\!\left(b (\partial_1c) (\partial_2^{2}b) (\partial_1c)\right)$ |
| $-\frac{27}{2}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2b) (\partial_1c) (\partial_1c)\right)$ |
| $\frac{3}{2}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1b) (\partial_2c) (\partial_1c)\right)$ |

$$
[\hbar^2]C_3=\frac14\operatorname{tr}\!\left\{b\left[(\partial_2^3c)(\partial_1^3c)-3(\partial_1\partial_2^2c)(\partial_1^2\partial_2c)+3(\partial_1^2\partial_2c)(\partial_1\partial_2^2c)-(\partial_1^3c)(\partial_2^3c)\right]\right\}.
$$

## $C_4$

$$
C_4=[\hbar^0]C_4+\hbar[\hbar]C_4+\hbar^2[\hbar^2]C_4,
$$

$$
[\hbar^0]C_4=\operatorname{tr}(b^3)a_1-\frac32\operatorname{tr}(b^2\partial_1c\,b^2\partial_2c)-\frac32\operatorname{tr}(b^2\partial_1c)\operatorname{tr}(b^2\partial_2c).
$$

### Coefficient of $\hbar^1$ in $C_4$

Apply

$$
\left(1-\frac{J_-J_+}{2}\right)\left(1-\frac{J_-J_+}{6}\right)
$$

to the coefficient-weighted sum of these 42 trace monomials, preserving the displayed factor order.

| Coefficient | Trace monomial |
|---:|:---|
| $-\frac{14463008783824}{996131895}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_2c)\right)\operatorname{tr}\!\left(b^{2} (\partial_1\partial_2c)\right)$ |
| $\frac{414853627740}{66408793}$ | $\operatorname{tr}\!\left(b^{2} (\partial_1\partial_2c) b (\partial_1\partial_2c)\right)$ |
| $\frac{12616544145859}{996131895}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) (\partial_1c)\right)\operatorname{tr}\!\left(b^{2} (\partial_1c)\right)$ |
| $\frac{966060767062}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1b) (\partial_1c) b (\partial_2c)\right)$ |
| $\frac{37723309127}{66408793}$ | $\operatorname{tr}\!\left(b^{2}\right)\operatorname{tr}\!\left(b (\partial_2^{2}c) (\partial_1^{2}c)\right)$ |
| $\frac{2032861225165}{132817586}$ | $\operatorname{tr}\!\left(b (\partial_2^{2}b) (\partial_1c)\right)\operatorname{tr}\!\left(b (\partial_1c)\right)$ |
| $-\frac{25908519850867}{996131895}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2b) (\partial_1c)\right)\operatorname{tr}\!\left(b (\partial_1c)\right)$ |
| $\frac{16955324342201}{3984527580}$ | $\operatorname{tr}\!\left(b^{2}\right)\operatorname{tr}\!\left((\partial_1^{2}b) (\partial_2c) (\partial_2c)\right)$ |
| $-\frac{365131791832}{66408793}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_2c) (\partial_1b) (\partial_1c)\right)$ |
| $\frac{83616196435}{6037163}$ | $\operatorname{tr}\!\left(b (\partial_2b) b (\partial_2c) (\partial_1^{2}c)\right)$ |
| $-\frac{11908180583773}{2988395685}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1c) (\partial_1c)\right)\operatorname{tr}\!\left(b (\partial_2b)\right)$ |
| $\frac{562762484551}{265635172}$ | $\operatorname{tr}\!\left(b^{2}\right)\operatorname{tr}\!\left((\partial_1b) (\partial_2^{2}c) (\partial_1c)\right)$ |
| $\frac{2760726060305}{66408793}$ | $\operatorname{tr}\!\left(b (\partial_2b) b (\partial_1\partial_2c) (\partial_1c)\right)$ |
| $\frac{59240431399}{132817586}$ | $\operatorname{tr}\!\left(b^{2}\right)\operatorname{tr}\!\left((\partial_1b) (\partial_1c) (\partial_2^{2}c)\right)$ |
| $\frac{1516560170}{548833}$ | $\operatorname{tr}\!\left((\partial_2b) b (\partial_1c) (\partial_2b) (\partial_1c)\right)$ |
| $\frac{444281621753}{132817586}$ | $\operatorname{tr}\!\left((\partial_1b) b (\partial_1c)\right)\operatorname{tr}\!\left(b (\partial_2^{2}c)\right)$ |
| $-\frac{1192822079472}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2b) b (\partial_1c) (\partial_1c)\right)$ |
| $-\frac{47000697216538}{2988395685}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_2c)\right)\operatorname{tr}\!\left(b (\partial_1b) (\partial_2c)\right)$ |
| $-\frac{15128950725731}{996131895}$ | $\operatorname{tr}\!\left((\partial_1^{2}b) (\partial_2c) b^{2} (\partial_2c)\right)$ |
| $-\frac{5259830085757}{1992263790}$ | $\operatorname{tr}\!\left((\partial_1^{2}b) b (\partial_2c)\right)\operatorname{tr}\!\left(b (\partial_2c)\right)$ |
| $-\frac{38786081057}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2b)\right)\operatorname{tr}\!\left(b (\partial_1c) (\partial_1c)\right)$ |
| $\frac{14812285312223}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_1b)\right)\operatorname{tr}\!\left(b (\partial_2c) (\partial_1\partial_2c)\right)$ |
| $-\frac{5385039753947}{996131895}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_2^{2}c)\right)\operatorname{tr}\!\left(b^{2} (\partial_1c)\right)$ |
| $\frac{25477327814321}{1992263790}$ | $\operatorname{tr}\!\left(b (\partial_1^{2}b) b (\partial_2c) (\partial_2c)\right)$ |
| $-\frac{48976031844446}{2988395685}$ | $\operatorname{tr}\!\left(b (\partial_2b)\right)\operatorname{tr}\!\left((\partial_1b) (\partial_2c) (\partial_1c)\right)$ |
| $-\frac{1002743742401}{996131895}$ | $\operatorname{tr}\!\left((\partial_2b) b (\partial_2c) b (\partial_1^{2}c)\right)$ |
| $\frac{2078764226204}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1b) b (\partial_1c) (\partial_2c)\right)$ |
| $-\frac{1011194068210}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1c) b (\partial_2b) (\partial_1c)\right)$ |
| $-\frac{5517472065494}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2c) (\partial_1c) b (\partial_1b)\right)$ |
| $\frac{13232032358514}{332043965}$ | $\operatorname{tr}\!\left(b^{2} (\partial_1b)\right)\operatorname{tr}\!\left((\partial_2c) (\partial_1\partial_2c)\right)$ |
| $\frac{29047206252619}{996131895}$ | $\operatorname{tr}\!\left(b^{2} (\partial_2b) (\partial_1\partial_2c) (\partial_1c)\right)$ |
| $\frac{444281621753}{66408793}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_1c) b (\partial_1\partial_2c)\right)$ |
| $\frac{11850829540168}{2988395685}$ | $\operatorname{tr}\!\left((\partial_1b) b (\partial_2c)\right)\operatorname{tr}\!\left((\partial_1b) (\partial_2c)\right)$ |
| $-\frac{23371221243386}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_2^{2}b) (\partial_1c) b (\partial_1c)\right)$ |
| $\frac{121188684971}{996131895}$ | $\operatorname{tr}\!\left((\partial_1b) b^{2} (\partial_2c) (\partial_1\partial_2c)\right)$ |
| $\frac{22475161760617}{1992263790}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_1^{2}c) b (\partial_2c)\right)$ |
| $-\frac{16341298110424}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_2b)\right)\operatorname{tr}\!\left((\partial_2c) b (\partial_1^{2}c)\right)$ |
| $-\frac{198086336691}{132817586}$ | $\operatorname{tr}\!\left((\partial_2b) b (\partial_1^{2}c)\right)\operatorname{tr}\!\left(b (\partial_2c)\right)$ |
| $\frac{25333224931111}{2988395685}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2c)\right)\operatorname{tr}\!\left((\partial_1b) b (\partial_1c)\right)$ |
| $-\frac{12896972854556}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_1\partial_2b)\right)\operatorname{tr}\!\left((\partial_2c) b (\partial_1c)\right)$ |
| $\frac{182565895916}{66408793}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_2c) (\partial_1c) (\partial_1b)\right)$ |
| $-\frac{20835815350171}{2988395685}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_2c)\right)\operatorname{tr}\!\left((\partial_1b) (\partial_1c)\right)$ |

### Coefficient of $\hbar^2$ in $C_4$

Apply

$$
\left(1-\frac{J_-J_+}{2}\right)\left(1-\frac{J_-J_+}{6}\right)\left(1-\frac{J_-J_+}{12}\right)
$$

to the coefficient-weighted sum of these 49 trace monomials, preserving the displayed factor order.

| Coefficient | Trace monomial |
|---:|:---|
| $\frac{105395696441}{199226379}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_1b)\right)\operatorname{tr}\!\left((\partial_2c) (\partial_1\partial_2^{2}c)\right)$ |
| $-\frac{6418825920644}{996131895}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}\partial_2c)\right)\operatorname{tr}\!\left(b (\partial_1\partial_2c)\right)$ |
| $\frac{8691631667194}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2c) (\partial_1^{2}\partial_2c) (\partial_1b)\right)$ |
| $-\frac{1915187807358}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}c)\right)\operatorname{tr}\!\left((\partial_2^{2}b) (\partial_1c)\right)$ |
| $-\frac{180409500189}{66408793}$ | $\operatorname{tr}\!\left((\partial_1^{2}b) (\partial_2^{2}c) b (\partial_1\partial_2c)\right)$ |
| $-\frac{87573313754}{199226379}$ | $\operatorname{tr}\!\left(b (\partial_2^{3}b)\right)\operatorname{tr}\!\left((\partial_1c) (\partial_1^{2}c)\right)$ |
| $-\frac{464869764342}{66408793}$ | $\operatorname{tr}\!\left((\partial_1\partial_2b) (\partial_2^{2}c)\right)\operatorname{tr}\!\left(b (\partial_1^{2}c)\right)$ |
| $-\frac{548409125078}{18111489}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_1\partial_2^{2}c)\right)\operatorname{tr}\!\left((\partial_1b) (\partial_2c)\right)$ |
| $-\frac{1978711683154}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1\partial_2c) (\partial_1b) (\partial_1\partial_2c)\right)$ |
| $\frac{882681246764}{597679137}$ | $\operatorname{tr}\!\left((\partial_2^{3}c) b^{2} (\partial_1^{3}c)\right)$ |
| $\frac{3531952147491}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}b) (\partial_1\partial_2c) (\partial_2c)\right)$ |
| $\frac{30558084450}{66408793}$ | $\operatorname{tr}\!\left((\partial_1\partial_2b) (\partial_1\partial_2^{2}c) b (\partial_1c)\right)$ |
| $-\frac{604466630774}{271672335}$ | $\operatorname{tr}\!\left((\partial_2b) b (\partial_1^{3}c) (\partial_2^{2}c)\right)$ |
| $-\frac{276183923886}{6037163}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1\partial_2c)\right)\operatorname{tr}\!\left((\partial_1\partial_2b) (\partial_1c)\right)$ |
| $-\frac{29479109157416}{2988395685}$ | $\operatorname{tr}\!\left((\partial_2^{3}b) (\partial_1c) b (\partial_1^{2}c)\right)$ |
| $-\frac{8174847382046}{996131895}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}c) b (\partial_1\partial_2^{2}c)\right)$ |
| $-\frac{1805046760108}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1\partial_2^{2}c) (\partial_1b) (\partial_1c)\right)$ |
| $-\frac{353568711105}{66408793}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) (\partial_1c)\right)\operatorname{tr}\!\left((\partial_1^{2}b) (\partial_2c)\right)$ |
| $\frac{26891024616}{6037163}$ | $\operatorname{tr}\!\left(b (\partial_1\partial_2b) (\partial_2c) (\partial_1^{2}\partial_2c)\right)$ |
| $-\frac{2855151673184}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1c) (\partial_2c) (\partial_1^{2}\partial_2b)\right)$ |
| $\frac{30558084450}{66408793}$ | $\operatorname{tr}\!\left(b (\partial_1\partial_2^{2}b) (\partial_1^{2}c) (\partial_2c)\right)$ |
| $-\frac{1900768372542}{66408793}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_2c) (\partial_1^{2}b) (\partial_2^{2}c)\right)$ |
| $\frac{852598523944}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1\partial_2c) (\partial_1^{2}c) (\partial_2b)\right)$ |
| $-\frac{989355841577}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2^{2}c) (\partial_1b) (\partial_1^{2}c)\right)$ |
| $-\frac{363074994730}{199226379}$ | $\operatorname{tr}\!\left(b (\partial_1\partial_2^{2}c) b (\partial_1^{2}\partial_2c)\right)$ |
| $-\frac{5709492261796}{2988395685}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_1^{2}c)\right)\operatorname{tr}\!\left(b (\partial_2^{3}c)\right)$ |
| $\frac{14268913753}{398452758}$ | $\operatorname{tr}\!\left(b^{2}\right)\operatorname{tr}\!\left((\partial_2^{3}c) (\partial_1^{3}c)\right)$ |
| $\frac{4290983068502}{271672335}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_1^{3}c) (\partial_2^{2}c)\right)$ |
| $\frac{353568711105}{66408793}$ | $\operatorname{tr}\!\left((\partial_1\partial_2b) (\partial_2c)\right)\operatorname{tr}\!\left((\partial_1^{2}b) (\partial_2c)\right)$ |
| $\frac{28269532613216}{2988395685}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) b (\partial_2c) (\partial_1^{3}c)\right)$ |
| $\frac{1167235907679}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}b)\right)\operatorname{tr}\!\left((\partial_2c) (\partial_1\partial_2c)\right)$ |
| $\frac{8174847382046}{2988395685}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{3}c) b (\partial_2^{2}c)\right)$ |
| $\frac{8691631667194}{597679137}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2c) (\partial_1^{3}b) (\partial_2c)\right)$ |
| $\frac{7107504205264}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_2b) (\partial_1^{2}c) (\partial_1\partial_2^{2}c)\right)$ |
| $-\frac{5836479994010}{199226379}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}\partial_2b) (\partial_1c) (\partial_2c)\right)$ |
| $\frac{1961910516539}{66408793}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) (\partial_1c) (\partial_1b) (\partial_1\partial_2c)\right)$ |
| $\frac{9322877660656}{2988395685}$ | $\operatorname{tr}\!\left(b (\partial_1b)\right)\operatorname{tr}\!\left((\partial_2^{3}c) (\partial_1^{2}c)\right)$ |
| $-\frac{26891024616}{6037163}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) b (\partial_1^{2}\partial_2c) (\partial_1c)\right)$ |
| $-\frac{18931281140306}{996131895}$ | $\operatorname{tr}\!\left(b (\partial_1\partial_2b) (\partial_1^{2}\partial_2c) (\partial_2c)\right)$ |
| $\frac{483136944516}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1^{2}c) (\partial_2^{2}b) (\partial_1c)\right)$ |
| $-\frac{725762664838}{66408793}$ | $\operatorname{tr}\!\left((\partial_1\partial_2^{2}b) (\partial_2c)\right)\operatorname{tr}\!\left(b (\partial_1^{2}c)\right)$ |
| $\frac{44400552291}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2^{2}b)\right)\operatorname{tr}\!\left((\partial_1c) (\partial_1^{2}c)\right)$ |
| $-\frac{1168354566305}{132817586}$ | $\operatorname{tr}\!\left((\partial_1\partial_2b) (\partial_1^{2}b) (\partial_2c) (\partial_2c)\right)$ |
| $\frac{785020815476}{66408793}$ | $\operatorname{tr}\!\left((\partial_1b) (\partial_2c) (\partial_2^{2}c) (\partial_1^{2}b)\right)$ |
| $-\frac{17721704596106}{996131895}$ | $\operatorname{tr}\!\left((\partial_1^{2}\partial_2b) (\partial_1\partial_2c) b (\partial_2c)\right)$ |
| $-\frac{1598884501582}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_1c) (\partial_2^{2}c) (\partial_1^{2}b)\right)$ |
| $-\frac{437591417792}{66408793}$ | $\operatorname{tr}\!\left((\partial_2b) (\partial_2c)\right)\operatorname{tr}\!\left((\partial_1^{3}b) (\partial_2c)\right)$ |
| $\frac{87573313754}{199226379}$ | $\operatorname{tr}\!\left((\partial_2^{3}b) (\partial_1^{2}c)\right)\operatorname{tr}\!\left(b (\partial_1c)\right)$ |
| $\frac{1659194256051}{132817586}$ | $\operatorname{tr}\!\left((\partial_2^{2}b) (\partial_1c) (\partial_2c) (\partial_1^{2}b)\right)$ |

The higher coefficients of these saved representatives vanish. Every table was reconstructed from the exported trace words, projected with the displayed polynomial, and compared coefficientwise over $\mathbb Q$ to the independently read-back full cocycle.

The scalar uniqueness statements are $26-18-7=1$ for $N_b+\text{loop order}=3$ and $129-97-31=1$ for $N_b+\text{loop order}=4$. They do not classify other homogeneous blocks or non-scalar lower-filtration sectors.
