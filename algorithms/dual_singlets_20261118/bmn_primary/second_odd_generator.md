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

## Odd columns in $USp(4)$

At level 30 each of the following odd columns vanishes over the rationals. The table gives the central calculation at length eleven. The replay checks all 19 allowed odd lengths in these four columns against their independent colour characters. The previously excluded $(6,5,5)$ gives a fifth zero column. Among the six single-additional-generator explanations of the symplectic index, only $(7,5,4)$ remains after these completed jobs. [Symplectic index argument](symplectic_655/README.md) · [All lengths and ranks](rank_two_odd_generators.json)

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
