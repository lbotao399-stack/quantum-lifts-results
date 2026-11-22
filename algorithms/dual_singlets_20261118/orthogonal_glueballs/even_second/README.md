# The seventh and eighth glueball powers of $SO(12)$ and $SO(14)$ survive to the third page

In the full relative holomorphic jet complex, take the glueball $S$ with the defining orthogonal trace. Its powers $S^7$ for $SO(12)$ at level 42 and $S^8$ for $SO(14)$ at level 48 have zero outgoing supercharge on the displayed representatives. Both remain nonzero after quotienting by incoming $d_1$ and $d_2$ images. The target is colour-reflection even. Averaging over the two components of the orthogonal group commutes with the supercharge.

In the stable tree quotient at the pure gaugino weight, each nonempty graded-symmetric even trace is the bilinear $S$. The pure target is consequently one-dimensional. There are zero finite-colour trace relations at the pure target weight and at the first dual weight, even before flavour projection. The stable first dual therefore descends. The stable second dual obeys the two transposed descent identities tested on 33 first and 1,266 second source columns. Its finite-colour descent requires its value on every relation at the second dual weight. The calculation uses the full balanced dotted weight and projects the three scalar flavours to the singlet. This contains the dotted singlet of the target. [Stable dual construction](../../orthogonal_glueball_second_pair.py) · [Character and polarization](../../orthogonal_glueball_relation_span.py)

| Group | Power | Character dimension of relations | Independent ordered-word rows | Cartan rank | Rank with second pairing |
|---|---:|---:|---:|---:|---:|
| $SO(12)$ | 7 | 114 | 114 | 3 | 3 |
| $SO(14)$ | 8 | 4 | 4 | 1 | 1 |

Each row is an explicitly polarized Pfaffian of degree exceeding the defining dimension. In twelve colours, 263 polarizations produce a basis with a nonzero 114-by-114 minor in ordered traces. In fourteen colours the first four suffice. Their ordered-word minor is

$$
\det\begin{pmatrix}
0&0&144&240\\
-384&10464&0&0\\
0&0&0&-576\\
0&128&0&0
\end{pmatrix}=-4076863488.
$$

In the order of the 257 distinct-plane Cartan coordinates in the data, the second pairing on every relation is the linear functional with first three coefficients $(-1/320,-3/64,-5/64)$ for twelve colours and with first coefficient $1/1280$ for fourteen colours. All its remaining coefficients vanish. Subtracting this functional preserves the adjoint equation since the commuting coordinates annihilate tree boundaries. The corrected dual descends and annihilates all earlier images. [Twelve-colour basis and minor](span12.json.gz) · [Fourteen-colour basis and minor](span14.json.gz)

A8 job `5693393` performed the initial second-pairing search. A8 job `5693483` constructed both full minors. The following reproduces the physical dual, the characters, the selected Pfaffians, their ordered-word coefficients and their Cartan restrictions, using the repository trace runtime and Python with NumPy, SciPy and python-flint.

~~~sh
python3 verify.py --rebuild --workers 4
~~~
