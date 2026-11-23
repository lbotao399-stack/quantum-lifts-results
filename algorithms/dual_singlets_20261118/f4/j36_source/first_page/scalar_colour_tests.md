# Finite-colour identities for the sixth glueball of $F_4$

The second supercharge differential into $S^6$ at $j=36$ is open. The five-dimensional space under test consists of **marked scalar-commutator images**. Each tested vector includes its complete trace polynomial. Its five coordinates and the corrected row $\eta Q_1$ are specified in the [preceding calculation](trace_obstruction.md).

## Removing common glueball factors

Let $U,V$ be the two first dotted gauginos and use the original defining trace $S=\operatorname{Tr}_{26}(UV)$ for counting spectator factors. The 19 trace sources have at least two factors $S$ after reducing the inactive classical splittings. These splittings impose the three graded ghost brackets and the six ghost–scalar brackets. Each has a primitive with zero $Q_1$. Test a native identity with a common factor $S^k$ **before** these splittings. The two odd colour spaces together have dimension 104. Their trace pairing gives a nondegenerate exterior two-form. Exterior Lefschetz injectivity for a symplectic space of dimension 104 gives injectivity of multiplication by $S^k$ on odd degree $d$ whenever $d+k\le52$. This applies coefficient by coefficient in the commuting scalar variables. A raw identity divisible by $S^k$ therefore comes from a raw identity on the shorter words.

Enumerate stable spin and flavour singlets on the shorter words. Evaluate their full ordered traces on integer $F_4$ matrices with exterior coefficients. Reduce each stable word by the inactive splittings, restore its factors of $S$, and impose the polarized quartic identities of the original calculation. Every genuine raw identity is in the kernel of the integer evaluation. The reduction vanishes on **the entire computed kernel** in each row.

| Common factors removed | Raw length | Stable trace singlets | Full native character | Integer evaluation rank | Evaluation nullity | Rank after restoring factors and reducing |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 6 | 5 | 4 | 4 | 1 | 0 |
| 3 | 8 | 41 | 22 | 22 | 19 | 0 |
| 2 | 10 | 370 | 162 | 153 | 217 | 0 |

The last evaluation can leave up to nine native tensor directions undetected. Its kernel is an upper bound on the space of true trace identities. Vanishing on this upper bound proves that **none of the raw defining-trace identities with two common glueball factors can supply the missing relation**. The quotient by inactive splittings can have torsion under multiplication by $S$. Raw primitives or relations with fewer common factors, and general invariant colour contractions, still enter the full problem. [Rational reduction](native_trace_relations.py) · [Frames and ranks](scalar_colour_certificates.json)

With one common factor removed the shorter word has twelve letters. The full colour character is 1,278. Restrict the lengths of the individual traces as indicated. The last column is the rank of the image of the entire integer-evaluation kernel under the same rational reduction.

| Maximum trace length | Stable singlets | Terms per sampled letter | Evaluation rank | Nullity | Reduced rank |
|---:|---:|---:|---:|---:|---:|
| 5 | 275 | 2 | 123 | 152 | 6 |
| 5 | 275 | 6 | 161 | 114 | 0 |
| 5 | 275 | 8 | 161 | 114 | 0 |
| 6 | 630 | 6 | 269 | 361 | 0 |

The first evaluation admitted two apparent nonzero obstruction pairings. Every genuine trace identity satisfies each of the two denser evaluations independently. Both give zero reduced rank. The last row excludes every raw identity in the one-factor sector whose remaining traces have length at most six. Longer traces and identities coupling sectors with different numbers of factors remain open. [Exact kernel data and frame seeds](scalar_colour_certificates.json)

## All ordinary scalar triples give zero

Fix the three scalar colours at arbitrary complex elements of $\mathfrak f_4$. The linear ghost–scalar brackets restrict both gauginos to their common Lie centralizer. The trace of a product can be evaluated on the successive simple quotients of the 26 under the commuting actions of the scalar-generated Lie algebra and its algebraic centralizer.

If the scalar-generated Lie algebra is solvable, it acts by scalars on each such quotient. Indeed, its image is a solvable ideal of the joint Lie algebra, and a solvable radical acts by scalars in an irreducible finite-dimensional complex representation. Every marked scalar commutator then has trace zero with any spectators. This includes nilpotent scalar triples generating a solvable algebra.

Otherwise its Levi factor contains an $\mathfrak{sl}_2$. The unipotent radical of the scalar centralizer acts by zero on each simple quotient, since it acts nilpotently and its Lie algebra acts there by scalars. Thus the trace polynomial factors through the reductive quotient of the centralizer. The projection to this quotient is a Lie map, so a primitive for a pure ghost bracket there pulls back to one on the actual fibre. The trace at the fixed triple is invariant under the connected reductive quotient.

A Levi subgroup of the centralizer commutes with the scalar $\mathfrak{sl}_2$, hence has rank at most three. For rank at most two the maximal number of each dotted ghost in an invariant class is at most three. This follows from the classical chiral-ring theorem, which gives the maximal glueball power $h^\vee-1$ for each simple factor, and from the exterior algebra on the central torus. At rank three the Cartans of the Levi and of the scalar $\mathfrak{sl}_2$ together form a full $F_4$ Cartan. Their roots are then ambient roots. The possible excessive types are $B_3$ and a $G_2$ factor with a rank-one complement. The latter requires root-length-squared ratio three, whereas the ambient ratio is two. The centralizer of a long-root $\mathfrak{sl}_2$ has twelve short and six long roots, of type $C_3$. For a short root it has twelve long roots, of type $A_3$. Neither admits $B_3$ with its eighteen roots. The native bracket test checks all 48 roots. The remaining rank-three types have maximal equal-dot invariant degree at most three, including central torus contributions. The desired degree is $(4,4)$. [Chiral-ring theorem](https://doi.org/10.1002/prop.202400027) · [Native root checks](scalar_fibre_controls.py)

Consequently **every fixed ordinary scalar triple kills the whole marked trace space in the inactive quotient**. A nonzero class of the corresponding polynomial module can vanish on every such fibre. The same test in $G_2$ has maximal degree one in the centralizer of a nonabelian scalar Levi. Yet for the unrestricted three scalars the earlier rational calculation gives ranks $10$ and $11$ before and after adjoining $S^2\operatorname{Tr}(X_1X_2X_3)$ in the defining representation. The ordered trilinear is half the commutator trace. The current replay recovers both ranks. This is a control for the inactive quotient. [Full $G_2$ calculation](../../README.md#restricting-the-scalar-correction-to-simple-factors) · [Rational separating covector](scalar_colour_certificates.json)

There is a polynomial strengthening on **every fixed three-dimensional colour subspace**. The scalar singlet at weight $(2,2,2)$ is the square of the determinant in the coefficients of a scalar triple in a chosen basis. Its value is this square times its value on the basis triple. A single determinant times each basis ghost–scalar bracket is a polynomial combination of the three varying ghost–scalar brackets, by the adjugate identity. Multiplying the fixed-triple primitive by the determinant squared therefore gives a polynomial primitive on the whole subspace. This also holds after base change to a ring with nilpotents.

## Polynomial tests on larger colour subspaces

On a fixed colour subspace keep the coefficients of all three scalars as independent variables. The scalar singlets are products of two $3$-by-$3$ minors. In a relation replace one row of a minor by the three components of the ghost–scalar bracket; multiply by another scalar minor at degree two. Straighten the resulting products of minors. Every row constructed in this way is an actual inactive relation. In the seven-dimensional calculation their right kernel has dimension 260 at both test primes. For two identical ghosts a dual tensor must be alternating. Contracting either ghost slot gives a member of that one-ghost kernel. Imposing antisymmetry on twelve of the 490 scalar coordinates already leaves zero solutions. The full two-ghost quotient on this colour subspace is therefore zero over the rationals, and multiplication gives the same vanishing at four identical ghosts. Full rank modulo one prime suffices for this conclusion. The same full-rank condition holds on a Zariski-open family of seven-dimensional subspaces.

| Colour dimension | Degree in scalar minors | Relation rows | One-ghost columns | Certified rank lower bound |
|---:|---:|---:|---:|---:|
| 4 | 1 | 312 | 208 | 208 |
| 5 | 1 | 520 | 520 | 518 |
| 5 | 2 | 3,900 | 2,600 | 2,600 |
| 6 | 2 | 10,920 | 9,100 | 9,100 |
| 7 | 2 | 25,480 | 25,480 | 25,220 |

The defining-matrix commutators reproduce all 37,232 selected Ward rows at a separate integer assignment. A Cartan control gives $4!\,(5184)(1728)=214990848$ for $S^4$ times the scalar Gram determinant. This detects a commuting scalar singlet on a special four-plane. A8 jobs 5704264 and 5704343 replayed the six- and seven-colour minors and the seven-colour alternating quotient. [Ward matrix](scalar_plane_ward.py) · [Exterior replay](scalar_plane_exterior.py) · [Pivots and receipts](scalar_colour_certificates.json)

Two special restrictions also admit exact tests over polynomial scalar variables. Take the component of both dotted species of top degree in the central torus of the Levi. For $B_2$ keep all three scalars in the semisimple factor. For $B_3$ fix one scalar along the central torus and keep two scalars in the semisimple factor. The defining traces span the orthogonal invariant tensors here. Compute the native relations on the full ordered traces, project the relations into the stable ghost-commuting quotient, and compare with imposing scalar commutation as well.

| Factor and remaining weight | Stable singlets | Finite singlets | Finite relations | Ghost quotient | Quotient after scalar commutation |
|---|---:|---:|---:|---:|---:|
| $B_2$, $(2,2\mid2,2,2)$ | 370 | 104 | 266 | 0 | 0 |
| $B_3$, $(3,3\mid2,2,0)$ | 370 | 328 | 42 | 1 | 1 |

Both marked images vanish in these components. The orthogonal calculation uses rational kernels and integer evaluation frames attaining the full cochain characters. [Calculation](scalar_orthogonal_walls.py) · [Reduced relation matrices](scalar_colour_certificates.json)

## Replay

Use Python with NumPy, SciPy, SymPy and python-flint from this directory. The short native runs construct the integer frames and rational null spaces from scratch.

~~~sh
for q in 0 1 2; do python3 native_trace_relations.py --dotted "$q" --out "/tmp/f4-common-$q.json"; done
python3 native_trace_relations.py --dotted 3 --max-trace 5 --sample-terms 2 --out /tmp/f4-sparse.json
python3 native_trace_relations.py --dotted 3 --max-trace 5 --sample-terms 6 --seed 369112 --out /tmp/f4-dense.json
python3 native_trace_relations.py --dotted 3 --max-trace 5 --sample-terms 8 --seed 850774 --out /tmp/f4-second.json
python3 native_trace_relations.py --dotted 3 --max-trace 6 --sample-terms 6 --seed 369112 --out /tmp/f4-six-letter-traces.json
python3 scalar_fibre_controls.py
python3 scalar_plane_ward.py --size 6 --power 2 --grading 7 --out /tmp/f4-six.json
python3 scalar_plane_ward.py --size 7 --power 2 --grading 7 --out /tmp/f4-seven.json
python3 scalar_plane_exterior.py /tmp/f4-seven.json --prime 1009 --out /tmp/f4-exterior.json
python3 scalar_plane_exterior.py /tmp/f4-seven.json --prime 1013 --out /tmp/f4-exterior-second.json
python3 scalar_orthogonal_walls.py --colours 5 --out /tmp/b2-wall.json
python3 scalar_orthogonal_walls.py --colours 7 --out /tmp/b3-wall.json
python3 ../../../f4_subgroup_equation.py --group G2 --power 2 --workers 1 --out /tmp/g2-control.json
~~~
