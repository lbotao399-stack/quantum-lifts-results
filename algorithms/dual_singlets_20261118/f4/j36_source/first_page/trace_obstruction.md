# The corrected second differential on $F_4$ trace sources at level 36

For the sixth power of the native glueball, the corrected obstruction on a classically closed source $X$ equals $-\eta Q_1X$. Here $\eta$ is the rational covector satisfying $\lambda Q_1=\eta Q_0$ for the target covector $\lambda$. The intermediate classical map is surjective. A classical cycle with nonzero corrected pairing would therefore establish $S^6\in\operatorname{im}d_2$. [Full relative complex and the target on the second page](README.md#one-row-decides-the-second-differential)

## Classical rank on the selected Taylor diagram

The source has 3,720 cochains and its classical split takes values in a 704-dimensional space. Alternating the four columns of heights $(4,4,3,3)$ evaluates the corresponding scalar polynomial on four elements of $\mathfrak f_4$. A nonzero rational polarization identifies this copy with the dotted/flavour singlet. With diagonal scalar parameters $(1,1,1)$ and $(2,3,5)$ its physical coefficients are $50400$ and $45360000$. For at most three traces, each four-slot column has two entries in a common trace. After commuting the scalars its alternation vanishes. The adjacent swap in a trace is a marked Lie bracket, which supplies a classical primitive. The selected short four-trace polynomials have a marked commutator from their construction. [Physical polarization](physical_flag.py) · [Certificate](physical_flag.json)

| Classical family included | Independent directions | Selected determinant modulo 1009 | Same determinant modulo 1013 |
|---|---:|---:|---:|
| Short traces and single traces in the 26 and 52 | 433 | 582 | |
| Products of two traces added | 532 | 869 | 671 |
| Products of three traces added | 542 | 373 | 809 |

The last increase consists of nine length patterns $(10,2,2)$ and one $(9,3,2)$, counting independent rows. Partitioning the oriented trace products before evaluation gives an independent integer check on the unreduced Young alternation. The two primes certify the same rational lower bound. The Cartan direction gives the upper bound 703. [Single traces](single_trace_flag.py) · [Double-trace minor and A8 receipts](double_a8.json) · [Triple-trace minor](triple_selected.json) · [A8 receipt](triple_a8.json)

## Moving the gauginos on colour graphs

Write $U,V$ for the two dotted first gauginos in the trace algorithms. Classical splitting of second ghosts and differentiated scalars imposes their graded commutators with each other and their commutators with first scalars. Each such split has a primitive with zero $Q_1$, since no pair of its letters contains three flavour slots. In a colour graph the Ward identity at each invariant vertex gives the signed incidence relation on placing a gaugino on an incident edge. Moving past scalar spectators uses the above splits. For each dotted species, placements anticommute and their squares vanish. Eliminating a spanning forest therefore leaves one independent placement per graph cycle. The fourth exterior power is zero on graphs with fewer than four cycles. A six-vertex $K_{3,3}$ graph has four cycles. In a fundamental cycle basis its four-edge minors have multiplicities $39,45,42$ at values $-1,0,1$. Both dotted groups together have $81^2=6561$ nonzero placement pairs. An invariant cubic in the 26 has 913952 verified Ward entries. Its $K_{3,3}$ Cartan evaluation divided by the undecorated graph is $972/2275$. [Graph calculation](colour_graph_ghosts.py) · [Native values](colour_graph_ghosts.json)

The quantum colour tensor on a matter fermion and its complementary scalar joins the two slots by the inverse metric on the adjoint. The two emitted ghosts act by adjoint commutators on this edge. In the intermediate quotient their adjoint actions square to zero and anticommute. Moving another ghost across either active slot gives a repeated action on the edge and hence zero. The intermediate covector $\eta$ annihilates the discarded ghost splittings by $\eta Q_0=\lambda Q_1$. Thus both the reduced classical split and the corrected score permit moving ghosts past the matter fermion. A graph of at most three cycles gains at most one cycle from the exchange. Its five identical outgoing ghosts have zero projection. [Native pair tensor](colour.py) · [Full-length forward kernel](../first_action.py)

The 542 selected classical columns have at most four traces. After commuting the ghosts, four traces have four independent formal scalar-commutator directions. The order-four minor on the short rows is $291259281309696000$. The character count gives the matching upper bound of four before finite-colour specialization. Thus at least $542-4=538$ independent classical directions admit a primitive with vanishing one-loop action. A candidate extension of the dual row has at most $704-538=166$ free coordinates before the Cartan freedom is removed. [Character and trace counts](ghost_trace_quotient.json) · [Four-row determinant and coordinates](flag_ghost_rows.json)

## Trace identities and corrected scores

After commuting $U,V$, a trace contains at most one of each. The source has one matter fermion and four first scalars. A trace with a single ghost needs a further letter. A trace with both ghosts and just one other letter vanishes by orthogonal reversal. Counting eight ghosts against five other letters forces at least two factors $\operatorname{Tr}_{26}(UV)$ in every reduced defining-trace source monomial.

At the terminal weight there are 14 formal scalar-commutator singlets. Polarize the $F_4$ quartic trace identity on four Lie words and multiply by all spectator traces of the required weight. Exact elimination leaves five commutator directions. The dimensions using at most four, five, six, and seven traces are $4,5,5,5$. This calculation uses all three scalar Lie alphabets and the two ghost arguments with their Koszul signs. Finite-colour identities involving other invariant tensors can reduce this quotient further. [Quartic identity on the Cartan](quartic_ghost_relations.py) · [Exact ideal and Weyl multiplicities](trace_ideal.py) · [Rank at every weight](trace_ideal.json)

The same five-coordinate presentation applies to trace polynomials made using other colour representations. Indeed, symmetrize a word one adjacent swap at a time. Each commutator term has one fewer Lie argument. After moving the ghosts the longest terminal trace has eight arguments. The $F_4$ Weyl invariant degrees are $2,6,8,12$. The defining traces at degrees $2,6,8$ supply the primitive generators through degree eight by exact Cartan minors. Polarization expresses every symmetrized trace in these generators. Induction on the word length then expresses each ordered trace in the defining-trace presentation. For example,

$$
\operatorname{Tr}_{52}(X^6)=-7\operatorname{Tr}_{26}(X^6)+\frac5{36}\big(\operatorname{Tr}_{26}(X^2)\big)^3 .
$$

The same induction with a marked scalar commutator preserves its vanishing when the scalars commute. The reduction applies to products of traces. General networks with independent invariant vertices also occur in the full complex. [Cartan minors, degree-eight formulas and two-matrix checks](trace_representations.json) · [Replay](trace_representations.py)

At weight $(4,4\mid2,2,2)$, 162 source monomials remain. Exactly 33, 99 and 30 contain two, three and four pure defining traces of $UV$, respectively. The three simple raising maps give 19 singlets. To lift them to actual trace polynomials, alternate the two ghosts in each dotted pair and average the position of that block around the other cyclic letters. The resulting 19 polynomials satisfy all 57 raising identities before any ghost splitting. The classical split replaces the matter fermion by its two-scalar commutator. The rank on these singlets is 14 before the quartic and five afterwards. Each of the 14 resulting kernel vectors can be completed to a classical cycle using the ghost-splitting primitives. [Basis and kernel coordinates](trace_source_cycles.json.gz) · [Full trace singlet test](trace_source_section.py)

For the exact pairing, contract the intermediate covector with two copies of $\operatorname{Tr}_{26}(UV)$. The direct two-by-two metric minor and two successive single contractions agree. Push the two-input one-loop colour tensor backwards and evaluate the remaining traces using $B=\operatorname{Tr}_{26}/6$. There are 2747440, 770934, and 35809 distinct colour assignments for the three glueball multiplicities. Eight singlets have nonzero individual scores. On all 14 constructed classical cycles the corrected score is zero. The classical and augmented ranks are both five. [Corrected covector](trace_score_frame.py) · [Native contraction](trace_source_contract.cpp) · [All scores](trace_source_score.json)

In the five quotient coordinates listed in the certificate, the exact row for $\eta Q_1$ is

$$
4961442649363200\,(0,16,36,0,3).
$$

The selected five-column minor is $-17414258688$. Thus a further native relation in these coordinates with nonzero pairing would give a $d_2$ lift. All five coordinates refer to the displayed trace quotient. The calculation leaves the full $F_4$ differential open. An independent forward action on the original thirteen letters gives $-564,-13608,-1119744$ on three probes, one for each removed trace power. The shortened computation reproduces each coefficient. [Five trace coordinates, determinant and factorization](trace_score_descent.json) · [Forward and transposed controls](trace_score_controls.json)

## Replay

Run in this directory with Python, NumPy, SymPy and a C++17 compiler. The published intermediate frame can also be rebuilt from Dot's pinned archive.

~~~sh
python3 trace_ideal.py | diff -u trace_ideal.json -
python3 trace_representations.py | diff -u trace_representations.json -
python3 trace_source_section.py | diff -u trace_source_section.json -
python3 trace_source_score.py | diff -u trace_source_score.json -
python3 trace_score_descent.py | diff -u trace_score_descent.json -
python3 trace_score_frame.py --archive /path/to/F4-j36-d2-verification.zip --out /tmp/f4-frame.gz
cmp trace_score_frame.json.gz /tmp/f4-frame.gz
python3 trace_score_controls.py /path/to/F4-j36-d2-verification.zip | diff -u trace_score_controls.json -
OPENBLAS_NUM_THREADS=1 python3 triple_trace_flag.py --prime 1013 --replay triple_selected.json --workers 8 --out /tmp/f4-triple.json
~~~
