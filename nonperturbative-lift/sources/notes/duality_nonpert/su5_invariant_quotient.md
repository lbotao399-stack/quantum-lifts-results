# One traceless SU(5) matrix: exact invariant quotient

令 $X\in\mathfrak{sl}_5(\mathbb C)$，$\operatorname{Tr}$ 为 defining five-dimensional representation 的 ordinary trace，$I$ 为 identity matrix；按题目固定

$$
F(X)=X^2-\frac{\operatorname{Tr}X^2}{5}I.
$$

这里 $(F)$ 表示 $F$ 的 matrix entries 生成的 polynomial ideal，未取 radical。结论是

$$
\left(\mathbb C[\mathfrak{sl}_5]/(F)\right)^{SU(5)}
=\mathbb C[\operatorname{Tr}X^2]\big/\big((\operatorname{Tr}X^2)^3\big),
\qquad
[(\operatorname{Tr}X^2)^2]\ne0.
$$

## Completeness of the invariant ideal

$SU(5)$ 的 polynomial invariants/covariants 等于其 complexification $SL_5(\mathbb C)$ 的 invariants/covariants。Conjugation 下 scalar matrices 作用平凡，因此也等于 $GL_5(\mathbb C)$ 的结果。

先 restrict 到 diagonal matrices；dense diagonalizable locus 上 invariant polynomial 是 eigenvalues 的 symmetric polynomial。Newton identities 因而给出

$$
\mathbb C[\mathfrak{sl}_5]^{SU(5)}
=\mathbb C[\operatorname{Tr}X^2,\operatorname{Tr}X^3,
\operatorname{Tr}X^4,\operatorname{Tr}X^5].
$$

这四个 generators algebraically independent：它们与 characteristic polynomial 的四个非平凡 coefficients 之间是 invertible triangular polynomial transformation。

**Polynomial covariant completeness.** 任何 conjugation-equivariant polynomial matrix $B(X)$ 都与 $X$ commute。对 cyclic matrices，即存在一个 vector 使其前五个 iterates 构成 basis 的 matrices，centralizer 是 $I,X,\ldots,X^4$ 的 span。在 characteristic polynomial 的 standard companion matrix 上，这五个 powers 作用于 first standard basis vector 得到全部 standard basis vectors。因此 expansion coefficients 直接等于 $B$ 的 first-column entries；它们是 characteristic coefficients 的 polynomials。所有 cyclic matrices 都 conjugate 到该 companion matrix。Cyclic locus dense，故 polynomial identity 延拓到所有 $X$：

$$
B(X)=\sum_{\substack{j=0,\ldots,4\\\mathrm{matrix\ powers}}}
b_j(X)X^j,
\qquad
b_j\in\mathbb C[\operatorname{Tr}X^2,\operatorname{Tr}X^3,
\operatorname{Tr}X^4,\operatorname{Tr}X^5].
$$

$b_j$ 只是这一展开的 polynomial coefficients，没有改变 $X,F$ 的 normalization。

**Reynolds step.** 对 invariant polynomial $f\in(F)$，先写 $f(X)=\operatorname{Tr}(B(X)F(X))$，此时 $B$ 未必 equivariant。$dg$ 为 total mass $1$ 的 Haar measure；取 coefficientwise finite-dimensional average：

$$
\int_{SU(5)}g^{-1}B(gXg^{-1})g\,dg.
$$

它是 polynomial covariant，并且其与 $F$ 的 trace pairing 为

$$
\int_{SU(5)}
\operatorname{Tr}\big(B(gXg^{-1})F(gXg^{-1})\big)\,dg
=\int_{SU(5)}f(gXg^{-1})\,dg=f(X).
$$

代入 covariant completeness，且 $\operatorname{Tr}F=0$，得到完整 invariant ideal：

$$
(F)\cap\mathbb C[\mathfrak{sl}_5]^{SU(5)}
=\big(\operatorname{Tr}(XF),\operatorname{Tr}(X^2F),
\operatorname{Tr}(X^3F),\operatorname{Tr}(X^4F)\big).
$$

Reynolds average 同时保证 invariants functor 对本 quotient exact。因此没有遗漏由非 invariant coefficients 产生的新 invariant relation。

## Exact elimination

Cayley–Hamilton 与 Newton identities：

$$
X^5-\frac{\operatorname{Tr}X^2}{2}X^3
-\frac{\operatorname{Tr}X^3}{3}X^2
+\left(\frac{(\operatorname{Tr}X^2)^2}{8}
-\frac{\operatorname{Tr}X^4}{4}\right)X
-\left(\frac{\operatorname{Tr}X^5}{5}
-\frac{\operatorname{Tr}X^2\operatorname{Tr}X^3}{6}\right)I=0.
$$

乘 $X$ 后取 trace：

$$
\operatorname{Tr}X^6
=\frac34\operatorname{Tr}X^2\operatorname{Tr}X^4
+\frac13(\operatorname{Tr}X^3)^2
-\frac18(\operatorname{Tr}X^2)^3.
$$

四个 invariant generators 因而是

$$
\begin{aligned}
\operatorname{Tr}(XF)&=\operatorname{Tr}X^3,\\
\operatorname{Tr}(X^2F)&=\operatorname{Tr}X^4-\frac15(\operatorname{Tr}X^2)^2,\\
\operatorname{Tr}(X^3F)&=\operatorname{Tr}X^5-\frac15\operatorname{Tr}X^2\operatorname{Tr}X^3,\\
\operatorname{Tr}(X^4F)&=\frac{11}{20}\operatorname{Tr}X^2\operatorname{Tr}X^4
+\frac13(\operatorname{Tr}X^3)^2-\frac18(\operatorname{Tr}X^2)^3.
\end{aligned}
$$

直接 polynomial combination：

$$
\operatorname{Tr}(X^4F)
-\frac{11}{20}\operatorname{Tr}X^2\operatorname{Tr}(X^2F)
-\frac13\operatorname{Tr}X^3\operatorname{Tr}(XF)
=-\frac3{200}(\operatorname{Tr}X^2)^3.
$$

故完整 ideal 等于

$$
\left(
\operatorname{Tr}X^5,
\operatorname{Tr}X^4-\frac15(\operatorname{Tr}X^2)^2,
\operatorname{Tr}X^3,
(\operatorname{Tr}X^2)^3
\right).
$$

这也是 lex order $\operatorname{Tr}X^5>\operatorname{Tr}X^4>\operatorname{Tr}X^3>\operatorname{Tr}X^2$ 的 reduced Gröbner basis。Exact rational computation 独立核验其 normal forms：

$$
\mathrm{NF}\big((\operatorname{Tr}X^2)^2\big)=(\operatorname{Tr}X^2)^2,
\qquad
\mathrm{NF}\big((\operatorname{Tr}X^2)^3\big)=0.
$$

$\mathrm{NF}$ 表示除以上述 Gröbner basis 得到的 unique remainder。故 invariant quotient 的 vector-space basis 是

$$
1,\qquad [\operatorname{Tr}X^2],\qquad[(\operatorname{Tr}X^2)^2].
$$

## Reduced geometry is a different statement

若 $\operatorname{Tr}X^2\ne0$，$F=0$ 使 $X$ 的 eigenvalues 只有

$$
+\sqrt{\operatorname{Tr}X^2/5},\qquad
-\sqrt{\operatorname{Tr}X^2/5}.
$$

两 multiplicities 总和为 $5$，不可能相等，故 trace 不可能为零。于是 every complex solution 满足

$$
\operatorname{Tr}X^2=0,\qquad X^2=0,\qquad
\operatorname{Tr}X^j=0\quad(j\ge1).
$$

因此 **reduced invariant quotient** 是一个点。**Reduced matrix variety** 则是 square-zero matrices，含非零 matrices，例如只有 $(1,2)$ entry 等于 $1$ 的 matrix；不能称 entire matrix variety 是 origin。

该计算是 exact invariant-theory/F-term algebra。它本身不证明 quantum chiral ring 不再接收其他 relations，更不证明当前 nonperturbative differential。

Verification: `check_su5_invariant_quotient.py`; exact readback: `su5_invariant_quotient_readback.json`. All arithmetic is over the rationals; extension to complex coefficients preserves this Gröbner basis and its standard monomials.
