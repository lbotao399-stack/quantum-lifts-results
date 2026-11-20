# A first differential at level 24 in SO(7)

At doubled source charges $(2,2,12,6,6)$ and target charges $(1,1,12,6,6)$, apply two spin-to-first-flavour substitutions and then the first odd translation to the known source $F$ for the fourth gaugino power. Denote this ordered operation by its action on the two archived polynomials. It produces a 23,220-term highest-weight source and a six-term target from the fourth power. The classical variation of the source is the negative of the same operation on $Q_0F$. All 81 independent colour coordinates of $Q_0F$ vanish. The six-term target is classically closed term by term after collecting the variation. [Source and substitutions](../dual633_pairing.py)

The first action has 502 terms with common denominator 126. In 250 independent target coordinates, the classical images of all 476 formal highest-weight predecessors have rational rank 209. Adding the six-term target or the first action gives rank 210. Adding their sum gives rank 209. Thus the first supercharge differential on the displayed source is minus the six-term class. A8 job `5679131` and a second integer frame computed locally give the same five rational ranks. [A8 coordinates](a8_coordinates.json.gz) · [Independent coordinates](local_coordinates.json.gz) · [Determinants and frame definitions](a8_pairing.json)

The full filtered maps on both sides of the same doubled source charges give the following dimensions over the field with 1,000,003 elements. The four derivative bidegrees are $(t,t)$ with $0\le t\le3$. The incoming and outgoing ranks include every order allowed in these intervals.

| Group | Cochains | Incoming rank | Outgoing rank | Surviving dimension |
|---|---:|---:|---:|---:|
| $SO(7)$ | 14,302 | 8,349 | 5,803 | 150 |
| $USp(6)$ | 14,280 | 8,343 | 5,787 | 150 |

The orthogonal outgoing rank consists of 5,802 classical pairs and the first-order pair above. The other three maps have only classical pairs at this prime. Equality of the two surviving dimensions leaves the instanton question open. [Orthogonal outgoing](so1.json) · [Orthogonal incoming](so2.json) · [Symplectic outgoing](sp1.json) · [Symplectic incoming](sp2.json) · [A8 receipts](jobs.json)

From the repository root, run

```sh
python3 algorithms/dual_singlets_20261118/so7_level24_633/verify.py
python3 algorithms/dual_singlets_20261118/dual633_pairing.py --out /tmp/so7-level24-pairing
```

To regenerate the filtered maps, run `dual_charge_bands.py` with `--flavours 6 3 3`, the desired group and defining dimension, and `--deficit 1` or `--deficit 2`. Run each grade from 0 through 3 in the same output directory, then omit `--grade` to assemble the intervals. [Implementation](../dual_charge_bands.py)
