# Independent global SO11 return audit

## Result and unchanged conventions

The coefficient-functional identity is

\[
P^6G_C=25200R_C.
\]

Here the names refer to the accompanying SO11 calculation: `R_C` is precisely
the ordered-trace coefficient oracle `PfaffianCoefficient(3, 0)`;
`G_C` uses the same universal Pfaffian series, `PfaffianCoefficient(5, 0)`
scalar differential weights, exterior mask `0`, and no external minus sign.
`P^6` executes the six original odd scalar-to-ghost substitutions in the order
`A-, A+, B-, B+, C-, C+`. Execution order is stated explicitly; a transposed
coefficient-probe calculation executes them in reverse order.

This identity holds for every ordered trace coefficient in the indicated
degree, including products with repeated traces. The certificate checks the
complete auxiliary-marker functional, not only the two previously selected
trace probes. No differential, cochain, or Pfaffian normalization is changed.

## Reproduction

```sh
python -m pip install -r requirements.txt
python global_return_audit.py
```

Only Python and SymPy are required. The audit first looks for the pinned file
in the sibling `positive_sector/inputs/algorithms/dual_singlets_20261118/`
directory, then uses the unchanged copy bundled in `inputs/`. An explicit
source file and a separate receipt destination may be supplied:

```sh
python global_return_audit.py \
  --pfaffian-source inputs/orthogonal_glueball_first_pair.py \
  --out-dir replay
```

The source file must have SHA256
`eeb6d6d3acb94df435bbac3424c695abcd82156c569e2c9aa37be5676f553033`.
It comes from commit `9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9` of
`lbotao399-stack/quantum-lifts-results`, path
`algorithms/dual_singlets_20261118/orthogonal_glueball_first_pair.py`.
The audit compiles only its unchanged `scalar_code` and
`PfaffianCoefficient` AST definitions. It imports none of the candidate's
return routines, determinant routines, trace canonicalizers, or native
evaluation runtime. SymPy independently expands the determinant powers.
Python assertions must be enabled.

The released result contains:

| Exact check | Count or result |
|---|---:|
| Scalar-marker monomials in the complete degree | 1000 |
| Monomials with nonzero returned coefficient | 54 |
| Admissible plus/minus auxiliary-label assignments | 36 |
| Original factorial weights, determinant powers 3 and 5 | 54 and 222 |
| Original coefficient-method checks | 2376 |
| Trace multiplicity patterns through seven factors | 44 |
| Rejected exterior masks or double-ghost grade | 64 |
| Splittings of six ghosts across trace factors | 32 |
| Physical return multiplier | 25200 |

The complete 1000-entry coefficient table has SHA256
`4b17601b362e4a365fa0d5049aabc3f212211bc69f8e29f9e2e2db4e0df53586`.

## Global operator argument

### 1. The existing auxiliary variables and homogeneous degree

`X` denotes the 3-by-3 array of even scalar polarizations encoded by the
source's nine three-bit scalar exponents. Its rows label `A,B,C`; its columns
label the auxiliary matrices `H,E,F`. The six odd polarization markers have
the source's order

\[
\xi_{+0}\xi_{+1}\xi_{+2}\xi_{-0}\xi_{-1}\xi_{-2}.
\]

These describe the variables already used by the coefficient engine; they
do not introduce a new normalization. Each marker multiplies its odd colour
matrix. The universal Pfaffian series is obtained from
`exp(1/2 Tr log(1 + A))`, with `A` linear in these scalar and ghost
polarizations. The `doubles=0` sector contains no double-ghost tiles.

The scalar ancestor extracts scalar flavour degrees `(5,5,5)`, hence total
letter degree 15. The six returns remove two scalars of each flavour and
insert three plus and three minus ghosts. The result therefore has scalar
degrees `(3,3,3)` and ghost degrees `(3,3)`, still total degree 15. The
`PfaffianCoefficient(3,0)` extraction has exactly this same degree.
Neither side contains a coefficient in a different degree. In the SO11
application this is the degree-15 part of the same universal Pfaffian trace
identity; the argument does not transfer an identity between different
grades or between different `doubles` sectors.
With the source's two-dimensional auxiliary representation, the SO11
Pfaffian matrix has size 22. A Pfaffian of that linear matrix has degree at
most 11, so its homogeneous degree-15 coefficient is a finite-rank trace
identity. The formal degree-15 expression remains available before imposing
that identity, which is the setting of the coefficient calculation here.

### 2. Physical substitution and its transpose

On the linear auxiliary-colour matrix, the physical odd substitution
`phi_i -> lambda_sigma` agrees with

\[
\sum_{a=0}^{2}X_{ia}\frac{\partial}{\partial\xi_{\sigma a}},
\]

where the odd derivative is a left derivative. Both derivations vanish on
the ghost colour matrices and have the same value on the linear matrix.
The graded Leibniz rule therefore extends this equality to every matrix
power, trace, logarithm, and exponential in the universal series. It also
extends to their successive derivatives. Thus it is an identity before any
ordered colour-trace coefficient is selected.

For an even scalar differential polynomial, multiplication by `X_ia`
transposes to differentiation of that polynomial with respect to `X_ia`:

\[
\left.p(\partial_X)\bigl(X_{ia}f(X)\bigr)\right|_{X=0}
=
\left.(\partial_{X_{ia}}p)(\partial_X)f(X)\right|_{X=0}.
\]

This equality includes the ordinary monomial factorials. Applying it six
times transfers the six physical returns onto `(det X)^5`.

### 3. All odd assignments and the sign

An assignment survives only if the three plus labels use each of `0,1,2`
once, and the three minus labels do the same. There are `3! * 3! = 36`
assignments. With both label assignments the identity permutation, the
left derivatives remove the exterior markers in the order

\[
(3,0,4,1,5,2).
\]

Its inversion count is `3 + 2 + 1 = 6`, so its sign is `+1`. Permuting the
plus labels multiplies this sign by the plus permutation sign; permuting
the minus labels multiplies it by the minus permutation sign. Consequently
the six transferred scalar derivatives are exactly

\[
\det(\partial_X)^2(\det X)^5
=25200(\det X)^3,
\qquad
25200=(5\cdot6\cdot7)(4\cdot5\cdot6).
\]

The determinant expansion used by the independent check is

\[
\begin{aligned}
\det X={}&X_{00}X_{11}X_{22}+X_{01}X_{12}X_{20}
+X_{02}X_{10}X_{21}\\
&-X_{00}X_{12}X_{21}-X_{01}X_{10}X_{22}
-X_{02}X_{11}X_{20}.
\end{aligned}
\]

Its squared differential has precisely the 36 signed assignments above.
The code does not assume the displayed differential identity: it verifies
every coefficient-functional component. Each row of the remaining scalar
exponent matrix has degree 3, with `binomial(5,2)=10` possible exponent
triples. The complete component space therefore has `10^3=1000`
monomials. For each monomial, the program sums all 36 signed terms using
the independently expanded `(det X)^5` coefficient and its full monomial
factorial, then compares the sum with 25200 times the corresponding
factorial-weighted `(det X)^3` coefficient. All 1000 comparisons are exact
integer equalities. Outside these row degrees both functionals vanish.
Every universal ordered-trace coefficient is a linear combination of these
marker monomials, proving the identity globally in this degree.

### 4. Why the marker-insertion calculation has a different sign

Six successive left marker insertions reverse the order of the six markers
relative to the physical left-derivative removal calculation. Their extra
sign is

\[
(-1)^{6\cdot5/2}=(-1)^{15}=-1.
\]

Thus a separate calculation with left insertion operators gives
`-25200`, whereas the physical Fischer transfer above gives `+25200`.
These are different operator calculations; the insertion sign must not be
applied a second time to the physical return.

In the actual universal series, each odd marker is initially written
before its odd colour matrix. Moving all six markers to the left produces
another, already prescribed, sign `(-1)^15=-1`. The original upstream
`coefficient` method implements exactly this sign:

```python
value = sum(c*poly.get((63,x,self.doubles),0)
            for x,c in self.weights.items())
return -Fraction(value,multiplicity)
```

Accordingly `R_C` already includes the marker-colour sign. The physical
transfer multiplies this unchanged `R_C` by **positive** 25200. The scalar
ancestor has no odd markers and no such external minus sign.

### 5. Trace multiplicities, cyclicity, and reversal

The logarithm contributes the source's exact single-trace factor
`(-1)**(length+1)/(2*length)`. Expanding its exponential contributes

\[
\frac{1}{\prod_{\text{distinct traces}}
             (\text{multiplicity of that trace})!}.
\]

These factors are already part of the universal series on which the
derivation identity was proved. Repeated factors receive their full
Leibniz multiplicities when the physical returns act; no division by the
number of returns, or by `6!`, occurs. The determinant derivatives use the
same scalar monomial factorials as the original `weights` dictionary,
checked for all 54 and 222 nonzero entries. The original `coefficient`
method is additionally exercised on all 44 integer multiplicity patterns
through seven trace factors, for every one of its 54 weights. Seven is an
upper bound because a nonzero orthogonal trace has at least two letters
and the total degree is 15. This verifies 2376 exact rational factors.

Cyclic rotations and orthogonal reversal are imposed only after the
universal identity; their signs and the factor `1/(2*length)` are unchanged
on both sides. Thus coefficient extraction descends to the same canonical
graded trace algebra. Repeated odd colour-trace factors vanish there, as
usual; they do not create an additional factorial convention.

The external marker-colour sign is also independent of how the six ghosts
are distributed among trace factors. If the factors contain respectively
`k_1,k_2,...` odd letters, collecting markers within factors contributes
`sum k_t(k_t-1)/2` exchanges, while collecting them across factors
contributes `sum_{t<u} k_t k_u`. Their sum is always

\[
\frac{(\sum_t k_t)(\sum_t k_t-1)}{2}=15.
\]

All 32 compositions of six positive ghost counts are checked. Factors
with zero ghosts contribute zero to both exchange counts. This supplies
the same global minus sign for every trace factorization.

## Scope in the fortuitous calculation

This package proves the global return identity only. Combining it with
the separate complete scalar-quotient certificate gives the arithmetic

\[
9\left(\frac16\right)
+\frac{25600}{25200}\left(-\frac{189}{128}\right)
=\frac32-\frac32=0.
\]

The sibling `so11` package additionally supplies an explicit 16-term
`Q0(I)` preimage of `9 G_D + (64/63) G_C`. That completeness and preimage
are separate certificates; they are not inferred from the vanishing of
one selected pairing. Together with the original return action, these
certificates can establish that the aligned closed combination has zero
classical fortuitous obstruction. They do not establish that its finite
cohomology class vanishes, that a monotone component is nonzero, or that a
full quantum stable lift has been constructed.
