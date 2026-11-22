<a id="np-section-19"></a>

## 21. $\mathcal N=2$ SQCD: exact periods, duality, and deformation-induced lifting

Four-dimensional $\mathcal N=2$ supersymmetry has eight supercharges.

Its gauge multiplet includes a complex adjoint scalar; a matter hypermultiplet consists of a fundamental and an antifundamental $\mathcal N=1$ chiral multiplet.

On the Coulomb branch, adjoint scalar expectation values leave massless Abelian gauge fields.

The auxiliary curves below encode their exact low-energy couplings; their coordinates $x,y$ are curve coordinates, and $u$ labels the vacuum.

<a id="np-6074772dd3c3803a-1"></a>

### Curves and singular loci

Use the Seiberg–Witten matter-paper normalization.

Here $x,y$ are curve coordinates, $u$ is the source Coulomb curve coordinate, and $\Lambda_{N_f}$ is the strong scale for each flavor theory.

| $N_f$ | Curve $y^2=$ | Cubic discriminant |
|---|---|---|
| $0$ | $x^3-ux^2+\Lambda_0^4x/4$ | $\Lambda_0^8(u^2-\Lambda_0^4)/16$ |
| $1$ | $x^2(x-u)-\Lambda_1^6/64$ | $-\Lambda_1^6(u^3+27\Lambda_1^6/256)/16$ |
| $2$ | $(x^2-\Lambda_2^4/64)(x-u)$ | $\Lambda_2^4(u^2-\Lambda_2^4/64)^2/16$ |
| $3$ | $x^2(x-u)-\Lambda_3^2(x-u)^2/64$ | $-\Lambda_3^2u^4(u-\Lambda_3^2/256)/16$ |

Scale matching gives $\Lambda_0^4=m\Lambda_1^3$ and $\Lambda_1^3=m_2\Lambda_2^2$.

[Seiberg–Witten, §§11–14](https://arxiv.org/abs/hep-th/9408099) In the Dorey–Khoze–Mattis Pauli–Villars prescription for $N_f=3$, the curve coordinate is $\widetilde u$, while the physical coordinate $u=\langle\operatorname{Tr}A^2\rangle$ satisfies

$$
\widetilde u=u+\frac{(\Lambda_3^{\mathrm{SW}})^2}{432},
\qquad
(\Lambda_{N_f}^{\mathrm{SW}})^{4-N_f}
=4(\Lambda_{N_f}^{\mathrm{PV}})^{4-N_f}.
$$

The two physical singular values are therefore

$$
u=-\frac{(\Lambda_3^{\mathrm{SW}})^2}{432},
\qquad
u=\frac{11}{6912}(\Lambda_3^{\mathrm{SW}})^2.
$$

This operator-coordinate matching preserves the geometry of the source curve.

[Dorey–Khoze–Mattis](https://arxiv.org/abs/hep-th/9611016)

The electric and magnetic periods are $a,a_D$, their charges are $n_e,n_m$, and $S_i$ and $m_i$ are the flavor charges and masses.

A period is the integral of the source Seiberg–Witten one-form around a closed curve.

The two independent electric and magnetic integration cycles give $a,a_D$.

Their charge-weighted combination $Z$ is the supersymmetry central charge.

A BPS particle saturates the mass bound fixed by $Z$, giving the equality below.

A degeneration is a point where an integration cycle shrinks; a particle carrying its charge can become massless there.

The source central charge is

$$
Z=n_ea+n_ma_D+
\sum_{\substack{i=1,\ldots,N_f\\\mathrm{flavor\ charges}}}
\frac{S_im_i}{\sqrt2},
\qquad M_{\mathrm{BPS}}=\sqrt2|Z|.
$$

The flavor-dependent residues and monodromies of the periods determine charged BPS data.

Existence of stable particles also depends on the chamber; alignment of central-charge phases permits wall-crossing.

<a id="np-9fe85cd156b64d38-1"></a>

### Exact instanton coefficients

Use the period normalization of the pure-gauge paper, with its strong scale set to $1$:

$$
y^2=(x^2-1)(x-u),
\qquad
a(u)=\frac{\sqrt2}{\pi}\int_{-1}^{1}
\frac{\sqrt{u-x}}{\sqrt{1-x^2}}\,dx.
$$

> **Proof.**
>
> For $|u|>1$, the convergent binomial expansion is
>
> $$
> a(u)=\sqrt{2u}
> \sum_{\substack{k\ge0\\\mathrm{binomial\ series}}}
> \binom{1/2}{2k}\frac{(2k)!}{2^{2k}(k!)^2}u^{-2k},
> $$
>
> using
>
> $$
> \frac1\pi\int_{-1}^{1}\frac{x^{2k}}{\sqrt{1-x^2}}\,dx
> =\frac{(2k)!}{2^{2k}(k!)^2}.
> $$
>
> Series inversion gives the exact coefficients
>
> $$
> \lim_{a\to+\infty}a^2(u-a^2/2)=\frac14,
> \qquad
> \lim_{a\to+\infty}a^6(u-a^2/2-1/(4a^2))=\frac5{64}.
> $$
>
> In the Matone relation $u=\pi i(F-a\partial_aF/2)$, $F$ is the prepotential, and $F_k$ is its $a^{2-4k}$ coefficient.
>
> Hence
>
> $$
> 2\pi iF_1=\frac14,
> \qquad 4\pi iF_2=\frac5{64},
> \qquad
> F_1=\frac1{8\pi i},
> \quad F_2=\frac5{256\pi i}.
> $$
>
> The exact periods fix the prepotential data through $a_D=\partial_aF$.
>
> $\square$

These instantons change the special Kähler metric and period map.

[Seiberg–Witten](https://arxiv.org/abs/hep-th/9407087); [Matone](https://arxiv.org/abs/hep-th/9506102)

<a id="np-641ee493e77f2a26-1"></a>

### Monopole condensate induced by an adjoint mass

Add the $\mathcal N=1$ deformation $m\operatorname{Tr}\Phi^2$.

The light chiral fields on the monopole patch are $M,\widetilde M$, with local superpotential

$$
\widehat W=\sqrt2A_DM\widetilde M+mU(A_D).
$$

The scalar value of $U$ is $u$.

Here $M$ denotes the monopole field; the SQCD meson occurred in the preceding sections.

> **Proof.**
>
> The source magnetic period is
>
> $$
> a_D(u)=\frac{i\sqrt2}{\pi}(u-1)
> \int_0^1\sqrt{\frac{1-t}{t}}\frac{dt}{\sqrt{2+(u-1)t}}.
> $$
>
> Thus
>
> $$
> \left.\frac{da_D}{du}\right|_{u=1}
> =\frac i\pi\int_0^1\sqrt{\frac{1-t}{t}}\,dt
> =\frac i2,
> \qquad \left.\frac{du}{da_D}\right|_{a_D=0}=-2i.
> $$
>
> The F- and D-equations give
>
> $$
> a_D=0,
> \qquad M\widetilde M=i\sqrt2m,
> \qquad |M|^2=|\widetilde M|^2=\sqrt2|m|.
> $$
>
> Restoring the source scale gives
>
> $$
> u_{\mathrm{vac}}=\pm\Lambda^2,
> \qquad M\widetilde M=i\sqrt2m\Lambda,
> \qquad W_{\mathrm{vac},\pm}=\pm m\Lambda^2.
> $$
>
> The second vacuum uses a local dyon frame.
>
> $\square$

At $m=0$, the Coulomb branch survives.

The explicit deformation produces the two isolated vacua.

<a id="np-79a202a95b44b405-1"></a>

### Conformal $N_f=4$ and Argyres–Seiberg duality

In the microscopic Pauli–Villars coordinate of Dorey–Khoze–Mattis,

$$
\tau=\frac\theta\pi+\frac{8\pi i}{g_{\mathrm{PV}}^2},
\qquad q=e^{i\pi\tau},
$$

Here $q$ is the one-instanton factor.

The checked exact coefficients of the massless effective coupling are

$$
\lim_{q\to0}(\tau_{\mathrm{eff}}^{(0)}-\tau)
=\frac{4i}{\pi}\log2,
$$

$$
\lim_{q\to0}
\frac{\tau_{\mathrm{eff}}^{(0)}-\tau-(4i/\pi)\log2}{q^2}
=-\frac{7i}{15552\pi}.
$$

These coefficients depend on the microscopic coupling coordinate.

The massive geometry gives a more reparametrization-stable test.

For masses $(m,m,0,0)$, let $e_i$ be the three roots of the source elliptic cubic, with $e_1+e_2+e_3=0$:

$$
y^2=\prod_{\substack{i=1,2,3\\\mathrm{elliptic\ roots}}}
(x-e_i\widetilde u-e_i^2m^2),
$$

$$
\operatorname{disc}_x
=\prod_{\substack{1\le i<j\le3\\\mathrm{elliptic\ root\ pairs}}}(e_i-e_j)^2
\prod_{\substack{k=1,2,3\\\mathrm{elliptic\ roots}}}(\widetilde u-e_km^2)^2.
$$

At the three distinct double degenerations, the monodromy in each vanishing-cycle period basis is $\begin{pmatrix}1&2\\0&1\end{pmatrix}$.

The period transformation and flavor triality of S-duality preserve the central charge.

The scalar Coulomb generator remains $u$, and the Higgs branch has complex dimension $16-6=10$.

The Argyres–Seiberg curve for $SU(3)$ with six hypermultiplets uses $u,v$ of dimensions $2,3$, together with the dimensionless coupling $f$:

$$
y^2=(x^3-ux-v)^2-fx^6.
$$

On the patch $f=1,u=0$ with $\widetilde v\ne0$, take the source coordinates

$$
x=-\frac{i\widetilde x}{\widetilde v},
\qquad y=\frac{2\widetilde y}{\widetilde v},
\qquad v=2i\widetilde v,
$$

This gives the isolated rank-one $E_6$ curve

$$
\widetilde y^2=\widetilde x^3-\widetilde v^4.
$$

The resulting polynomial curve extends to $\widetilde v=0$.

The other Coulomb coordinate belongs to the gauged $SU(2)$.

The source indices $T(\mathbf3)=4,T(\mathbf2)=1$ give

$$
-2T(\mathbf3)+2T(\mathbf2)+k_{E_6}=0,
\qquad k_{E_6}=6.
$$

The complex dimensions of the two Higgs-branch descriptions are

$$
2\cdot3\cdot6-8-8=20,
\qquad
2\cdot2+22-3-3=20.
$$

The source also matches invariant generators and relations.

Its remaining twistor-line step concerns identification of the full hyperkähler metric.

[Argyres–Seiberg](https://arxiv.org/abs/0711.0054); [Gaiotto–Neitzke–Tachikawa](https://arxiv.org/abs/0810.4541)

<a id="bulk-translation"></a>

### Translation protects two infinite towers of bulk operators

Consider the flat-space holomorphic-topological twist of pure four-dimensional $\mathcal N=2$ Yang–Mills theory. Its supercharge is Khan's $Q_{\rm HT}=\overline Q^+_{\dot1}+Q^+_2$ and the holomorphic translation is $\partial=\partial_z$. In its tree complex the odd letters are adjoint $b,\partial b,\ldots$ and $\partial c,\partial^2c,\ldots$, restricted to constant gauge invariants. Each letter has ghost number one. The spins of $\partial^k b$ and $\partial^\ell c$ are $k+1$ and $\ell$, respectively, with $k\ge0,\ell\ge1$. Every monomial therefore has spin at least its ghost number. Use the Taylor coefficients of $b(z)$ and of $\tilde c(z)=c(z)-c(0)$, with the factorials of the source. [Khan, §§2–3](https://arxiv.org/html/2604.03500v1#S3)

**Theorem 21.1.** For a complex simple gauge algebra the complete relative tree cohomology at ghost number two has dimension one at each spin $s\ge2$ and zero at lower spins. For $\mathfrak{sl}_2$ its tower is generated by the translations of Khan's $X=-[\delta_{ij}b^i\partial c^j]$. The second tower consists of the translations of $Y=[\epsilon_{ijk}b^ib^jb^k]$, with $\delta_{ij}$ and $\epsilon_{ijk}$ the source color form and invariant three-form in its color basis. For a simple algebra the corresponding lines use its invariant form and its contraction with the Lie bracket. Both towers survive every order of a filtered bulk deformation that preserves holomorphic translation, raises ghost number by one, preserves spin, and has the indicated basic complex on its tree page. The deformation parameters have spin and ghost number zero.

> **Proof.** The tree differential on Taylor coefficients is the relative Chevalley differential of the current algebra in (3.65)–(3.66). Equivalently, the modes paired with $\tilde c$ have orders $p\ge1$ and bracket by adding orders. Their action on the modes paired with $b$, of orders $q\ge0$, is the adjoint action with order $p+q$. The latter modes commute with one another. This follows by expanding $Qb(z)=[\tilde c(z),b(z)]$ in the source color action and $Q\tilde c(z)=[\tilde c(z),\tilde c(z)]/2$. The factorial Taylor coefficients give ordinary convolution of mode orders. A scalar invariant bilinear in two adjoints of a simple algebra is a multiple of the invariant symmetric form. There are no invariant linear cochains, hence no ghost-number-two boundaries.
>
> An invariant alternating two-cochain is determined by its scalar coefficient on each ordered pair of modes. Its cycle condition is the cyclic sum of its evaluations on a bracket and a third mode. The color factor in each nonzero equation is the same nonzero invariant three-form. There are three possible field contents.
>
> 1. On two $b$ modes let $u_{q,r}=-u_{r,q}$ be the coefficient. Inserting a positive $c$ mode in the cycle equation gives $u_{p+q,r}=u_{q,p+r}$ for $p\ge1$. At a fixed sum of orders one can move every unit to one end. Antisymmetry then gives $u_{q,r}=0$.
> 2. On two $c$ modes at fixed total $s$, let $f(p)$ be the coefficient with orders $(p,s-p)$. Then $f(s-p)=-f(p)$. On three positive modes with orders adding to $s$ the cycle condition reads $f(p+r)=f(p)+f(r)$. Thus $f(p)=p f(1)$ for $1\le p<s$, and $s f(1)=0$. At $s=2$ antisymmetry already gives zero.
> 3. On a $c$ mode of order $p\ge1$ and a $b$ mode of order $q\ge0$, denote the coefficient by $a_{p,q}$. Two positive $c$ modes give
>
> $$
> a_{p+r,q}=a_{p,r+q}+a_{r,p+q},\qquad
> a_{p,q}=p\,a_{1,p+q-1}.
> $$
>
> Thus there is exactly one solution at $p+q=s-1$. The coefficient of the Taylor mode of order $s-2$ of the invariant pairing $b(z)\partial c(z)$ is precisely $p$ times the corresponding mixed pair. It is nonzero. This proves the complete ghost-number-two result.
>
> The invariant three-form on three $b(z)$ is closed by invariance. To detect a positive Taylor coefficient, put two independent color directions in the constant term and a third direction with nonzero bracket pairing in the desired jet. At order zero use three independent odd color components. Every tree boundary has a $c$ factor, since each nonzero tree variation introduces one. This proves all these pure-$b$ classes are nonzero at tree level.
>
> On any successive page let $d$ be the induced differential. A filtered translation induces a translation on the page with $[d,\partial]=0$. For example, at first order the coefficient of $[Q,\partial]=0$ relates the commutator of the correction with the original translation to a tree boundary, even when the representative for translation receives a correction. A filtered chain homotopy for this Ward identity gives the same induced identity. The generators $X,Y$ lie on $s=\mathrm{gh}$, so their outgoing targets at $\mathrm{gh}+1$ are empty. The translation identity gives zero outgoing differential for every derivative. The $X$ tower has no incoming source in ghost number one. All possible incoming sources for the $Y$ tower lie in the ghost-number-two space just computed and have zero outgoing differential. This argument repeats on every page. At fixed spin the tree space is finite dimensional. The successive solutions define representatives in the completed formal complex. $\square$

The direct $\mathfrak{sl}_2$ check expands the tree differential on every invariant bilinear at spins $2$ through $16$ and checks the pure-$b$ representatives through spin $16$. The recurrence proves the unbounded statement. [Exact exterior-algebra calculation](../checks/translation_jets.py)

<a id="khan-presentation-test"></a>

**Corollary 21.2.** In Khan's proposed presentation $A$, the even generator $X$ has $(s,\mathrm{gh})=(2,2)$ and the odd generator $Y$ has $(3,3)$. All other elements are generated from these by product and translation. Every translation-commuting product derivation of bidegree $(0,1)$ on $A$ is zero.

> **Proof.** The spaces $A^{2,3}$ and $A^{3,4}$ vanish by the spin bound. A derivation must annihilate $X,Y$, then their derivatives, then every product. In a multiplicative spectral sequence with this algebra as its full initial page, the first potentially nonzero differential would be such a derivation. Induction therefore gives zero on all pages under this presentation hypothesis. $\square$

Khan's explicitly proposed algebraic differential satisfies $Q_{\rm inst}X=0$, $Q_{\rm inst}(\partial X)=Y$ and $Q_{\rm inst}Y=0$. The source points out its translation commutator in §5. Theorem 21.1 excludes this particular bulk pairing already in the basic complex, independently of whether $A$ is the complete perturbative algebra. The corollary requires the full presentation to extend the conclusion to every operator. The algebraic cohomology computed for Khan's differential and the instanton corrections to the periods earlier in this chapter are compatible with these statements about the translation-preserving bulk supercharge. [Khan, equations (4.9)–(4.11), (5.1)–(5.2)](https://arxiv.org/html/2604.03500v1#S4)

Instanton charge gives an independent check on the proposed coefficient. The source assigns $r(\Lambda)=2$, $r(Q_{\rm HT})=1$, $r(\partial X)=0$ and $r(Y)=-3$. Thus its ansatz uses $\Lambda^2Y$ of charge $1$. In the ordinary integer-instanton expansion on flat space, the pure $SU(2)$ one-instanton scale is $\Lambda^4$, of spurion charge $8$. A scalar coefficient made from its $k$th power and neutral constants would require $8k-3=1$, or $k=1/2$. This excludes such a coefficient at integer instanton number for these operators. Choosing a sheet or inserting a further charged scale changes the hypotheses. The translation test applies independently of this scale-counting test. [Khan, (4.6)–(4.12) and (4.41)–(4.43)](https://arxiv.org/html/2604.03500v1#S4) · [Seiberg–Witten](https://arxiv.org/abs/hep-th/9407087)
