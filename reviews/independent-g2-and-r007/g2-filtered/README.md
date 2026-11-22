# A level-24 second-differential test in the relative G₂ jet complex

Let $Y$ denote the number of letters. We order weights by the three flavour degrees followed by the two dotted Taylor degrees. At $(5,2,2;1,1)$ and $Y=9$, the cohomology of the classical supercharge $Q_0$ over the rationals is one-dimensional. Its class reaches $E_2$, the cohomology of the first induced differential $d_1$. Two independent prime-field calculations then give a complete second primitive. [Li, Liu and Chang, quantum supercharge and marked product](https://arxiv.org/abs/2611.24160) · [Public classical character, trace and ordered-pair algorithms at ee9e064](https://github.com/lbotao399-stack/quantum-lifts-results/tree/ee9e06494a926fb88baeb6fd875cee36762f63a7/algorithms/dual_singlets_20261118)

## Characteristic-zero source

At lengths $8,9,10,11$ the source-weight cochain dimensions are $3215,2872,1222,176$. The modular ranks out of lengths $8,9,10$ are $1825,1046,176$. The last rank fills its target. Nilpotence therefore bounds the rational rank out of length nine by $1222-176=1046$.

At length eight, adjoining rational products of the single-graviton factors of degrees two and six to the incoming tree boundaries gives a minor of order $1390$. All its columns are closed. Thus the rational rank out of length eight is at most $3215-1390=1825$. Both bounds are attained by the modular minors. Consequently

$$
\dim_{\mathbb Q}H(Q_0)_{(5,2,2;1,1),9}=2872-1825-1046=1.
$$

The incoming first source at $(6,3,3;0,0),Y=9$ has dimension $733-307-426=0$. The outgoing first target at $(4,1,1;2,2),Y=9$ has dimension $944-796-148=0$. Both equalities are rational because the displayed lower rank bounds saturate the exact cochain dimension. The only incoming higher maps would require a negative dotted degree.

A first correction has weight $(4,1,1;2,2)$ and length eight. At that length the joined columns of tree boundaries and rational graviton products have rank at least $1328$. The outgoing tree rank is at least $796$, and $1328+796=2124$ is the cochain dimension. Every closed freedom in the first correction is therefore a sum of a tree boundary and leading terms of quantum products. The marked-product theorem gives a closed perturbative lift of each product, so this freedom has zero induced first differential over the rationals.

## Two full second-correction systems

The source has total Taylor degree eleven and nine positive-degree letters. Three contracted letters have total Taylor degree at most five. The three-input kernel $Q_2$ requires six flavour slots, so $Q_2 O_0=0$ on every source cochain. For the native classical cycle $O_0$, the corrected second obstruction is represented by $Q_1O_1$, where $Q_0O_1=-Q_1O_0$. The target weight is $(3,0,0;3,3)$ at length eight.

| Prime | First minor of order 796 | Nonzero entries of $Q_1O_0$ | Nonzero entries of $Q_1O_1$ | First correction terms | Second correction terms |
|---:|---:|---:|---:|---:|---:|
| 2147483647 | 48031520 | 944 | 115 | 796 | 101 |
| 1000033 | 439099 | 944 | 115 | 796 | 101 |

In each field, the packet verifies $Q_0O_2=-Q_1O_1$ on all 115 target coordinates. One target dual detects the single modular target class and kills the residual. It also kills the first images of all 46 modular homogeneous corrections. The stored first and second primitives use the generator and evaluation bases at their respective prime. The characteristic-zero rank of the corrected second differential requires independent evidence beyond these two zero reductions.

At $(5,2,2;2,0)$ and $(3,3,2;3,0)$ the additional length-nine source dimensions are $2986-1926-1060=0$ and $4345-3098-1247=0$ over the rationals.

## Replay

With Python 3.13, NumPy and python-flint, check the coefficient packets using

~~~sh
python3 replay_g2_filtered.py
~~~

For a complete regeneration from the pinned public trace runtime, install SciPy and SymPy as well, and use a C++ compiler. Set REPO to a checkout of the linked commit. The program builds the exceptional colour contractions using the defining seven-dimensional trace and its dual, independently checks them by summing the 14 colour generators, regenerates the kernels and ranks, checks closure of all product factors, and compares every saved quantum coefficient.

~~~sh
python3 replay_g2_filtered.py --rebuild --repo "$REPO" --workers 2
~~~
