# Two odd BMN generators in $G_2$

At level $j=2(E+J_L)=36$, an odd-letter-length generator has flavour partition $(6,6,6)$. At least two odd generators occur by level 44. By level 62 there are odd generators in at least two distinct flavour partitions. These are rational classical BMN bounds. [Index and product calculation](rank_two_odd_generators.json)

The primary index at $(6,6,6)$, signed by letter length, is $34961-34962=-1$. A primary is a cohomology class modulo the images of the three raising supercharges. All 293 odd classical blocks through level 28 vanish. Moving raising operations off an odd generator onto its dressing by the graded Leibniz rule leaves three partitions that can supply this primary. The entire odd column at $(6,5,4)$ was previously excluded. At $(6,6,4)$ the last two classical cohomology dimensions are

$$
21650-8255-13395=0\quad (11\text{ letters}),\qquad
18935-12339-6596=0\quad (13\text{ letters}).
$$

The nonzero minors of the incoming and outgoing classical supercharge matrices at prime $2147483647$ bound their rational ranks from below. Their sum attains the rational cochain dimension obtained independently from the colour and flavour characters. Nilpotence makes both bounds equalities. This leaves the singlet $(6,6,6)$ itself, at length 11, 13 or 15. [Index and lower levels](README.md) · [A8 length eleven](5691485/receipt.json) · [A8 length thirteen](5691572/receipt.json)

At $(8,7,7)$ the primary index is $-2$. Tensoring with $(6,6,6)$ shifts all three parts of a flavour partition by six, since that representation is the sixth power of the determinant. The unique complementary dressing is $(2,1,1)$, with even cohomology dimension one. The singlet generator can therefore account for at most one of the required two primary classes. At $(19,6,6)$ the primary index is $-1$. The unique dressing of $(6,6,6)$ would be $(13,0,0)$. Its cohomology is the degree-thirteen part of the one-scalar invariant ring, whose generators have degrees two and six. This part has dimension zero. Any number of copies of the singlet gives zero in this primary through multiplication and raising. [Ordinary and primary indices](7658750/receipt.json) · [Tensor-product bounds](generator_bounds.py)

If exactly two odd generators suffice for all the archived indices, the second one must contribute at $(8,7,7)$. There are 33 candidate partitions after the low-degree rank exclusions. Testing the ordinary and primary indices through levels 162 and 192, with first part at most 50, leaves 22. These counts concern the two-generator hypothesis. The exclusions at levels 44 and 62 hold independently of that hypothesis. [Partitions and first failed inequalities](rank_two_odd_generators.json)

## Four odd generators in USp(4) by level 48

The primary indices at $(12,5,5)$ and $(8,8,8)$ are $-3$ and $-2$ at levels 44 and 48. The two independent odd generators at $(10,4,4)$ can supply at most two classes in the first sector and zero in the second. If one additional generator suffices, it must supply at least one and two respectively. Every odd generator has degree at least 15 by the 297 rational rank equalities through level 28. Its diagram must fit in the three-row rectangle $(8,8,8)$. Every complementary dressing has degree at most nine, where the even cohomology dimension equals the ordinary index. Moving the raising operators onto the dressing modulo total descendants and applying the Littlewood–Richardson rule gives all the candidates below. [Two existing generators](../b2_level36/README.md#length-and-lower-products-without-the-published-parent) · [Ordinary and primary indices](rank_two_indices/) · [All inequalities and rank checks](rank_two_odd_generators.json)

| Proposed third generator | Maximum contribution at $(12,5,5)$ | Maximum contribution at $(8,8,8)$ |
|---|---:|---:|
| $(6,5,4)$ | 3 | 5 |
| $(6,5,5)$ | 2 | 2 |
| $(7,4,4)$ | 6 | 3 |
| $(7,5,3)$ | 5 | 7 |
| $(7,5,4)$ | 3 | 5 |
| $(8,5,2)$ | 6 | 2 |

All six odd columns have zero rational cohomology. The last one to close was $(7,5,4)$. At lengths $7,9,11,13,15$ the independent character dimensions equal the sum of the two adjacent modular rank lower bounds,

$$
52=0+52,\quad 5732=934+4798,\quad 29157=12070+17087,
\quad 21748=14416+7332,\quad 2256=2024+232 .
$$

Nilpotence makes these equalities exact over the rationals. This excludes the entire odd source column and proves the fourth-generator bound. Within this index test $(8,5,3)$ passes all earlier inequalities and first fails at level 48, where its capacity is one. The bound counts classical highest-weight generators and supplies no new operator representative. [A6 outer lengths](7658856/receipt.json) · [A8 length 11](5700336/receipt.json) · [A8 length 13](5700337/receipt.json)

### The earlier four level-30 columns

At level 30 each of the following odd columns vanishes over the rationals. The table gives the central calculation at length eleven. The replay checks all 19 allowed odd lengths in these four columns against their independent colour characters. The $(6,5,5)$ calculation gives the fifth zero column. [Symplectic index argument](symplectic_655/README.md) · [All lengths and ranks](rank_two_odd_generators.json)

| Flavour | A6 job | Cochains | Incoming rank | Outgoing rank | Cohomology |
|---|---:|---:|---:|---:|---:|
| $(6,5,4)$ | [7658843](7658843/receipt.json) | 12711 | 6652 | 6059 | 0 |
| $(7,4,4)$ | [7658842](7658842/receipt.json) | 9053 | 4418 | 4635 | 0 |
| $(7,5,3)$ | [7658855](7658855/receipt.json) | 14783 | 7148 | 7635 | 0 |
| $(8,5,2)$ | [7658854](7658854/receipt.json) | 6713 | 2747 | 3966 | 0 |

From the repository root, the rank, character, product and index checks use the Python standard library.

```sh
d=algorithms/dual_singlets_20261118/bmn_primary
python3 "$d/verify.py"
python3 "$d/rank_two_odd_generators.py" --out /tmp/n4-rank-two-odd-generators.json
cmp "$d/rank_two_odd_generators.json" /tmp/n4-rank-two-odd-generators.json
```
