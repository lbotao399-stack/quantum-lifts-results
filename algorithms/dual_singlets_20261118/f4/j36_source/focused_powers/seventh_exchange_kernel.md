# Two marked endpoints of the F₄ seventh-power second source

At $j=42$, select the Taylor diagram $(5,5,3,3)$ and its ordinary five-scalar polarization. On each endpoint of the Albert contraction alternate a column of five and then a column of three. The first column contains either a marked second jet $Y_{ij}$ or its classical split $[X_i,X_j]$. Between the endpoints use the identity, one of the six pair exchanges, or two ordered exchanges. There are $1+6+36=43$ middle operators. Marking either end gives 86 source recipes. Every closed linear combination is a classical boundary. [Network and point generator](seventh_crossings.py) · [Exact identities and replay](seventh_exchange_relations.py)

For a fixed middle operator the two marks give the same outgoing classical column. Their difference is the classical derivative of the corresponding two-mark cochain by the graded product rule. Among the 43 common outgoing columns, the three disjoint pairs of exchanges commute already on free marked tensors. The fourth identity comes from the right endpoint.

Use the zero-based tensor slots 1 and 3 on that endpoint. Before either alternation the tensor in these slots is a raised Albert cubic with its third argument fixed by the other slots. Identify the pair with an endomorphism using the native invariant form. It is symmetric, and its endomorphism is a scalar multiple of the corresponding Jordan multiplication matrix $L$. For any two native Lie generators $U,V$ one has $\operatorname{Tr}_{26}(LUV)=0$ by all $26\cdot52^2=70304$ integer tests. All Lie generators are skew for this form.

In the first alternation there is one scalar on each of slots 1 and 3. The other scalar and the mark act on the complementary slots. Upon exchanging the two displayed scalar flavours the contribution, up to its common sign, is $ULV-VLU$. It is skew. Testing against any Lie generator $D$ gives a trace of $L$ against

$$
VDU-UDV=[V,D]U+D[V,U]+[D,U]V .
$$

Every trace on the right vanishes by the preceding integer identity. Thus the intermediate skew endomorphism is trace-orthogonal to the Lie algebra. The final alternating triple again acts once on each of slots 1 and 3. It gives a symmetric endomorphism whose trace pairs that intermediate endomorphism with the commutator of two Lie generators, hence is zero. This argument works before or after replacing the right mark by its bracket.

Write $E$ for the native integer exchange on the pair, $P$ for swap, and $h$ for the supplied native metric of the 26. Flattening the two slots in the native order gives the exact $676$-dimensional matrix identity

$$
(E-6I)(E+36I)(I+P)=240\,\operatorname{vec}(3h^{-1})\operatorname{vec}(h)^{\mathsf T}.
$$

The right endpoint is symmetric and trace-free in this pair. Consequently $E^2+30E-216I$ annihilates it. This yields the independent fourth relation on both choices of free mark.

The fixed order-39 outgoing minor equals $721\bmod1009$ and $842\bmod1013$. These values and the four rational identities give rational rank 39. In the 86-dimensional space of formal recipes the classical kernel has dimension $86-39=47$. It consists of the 43 changes of mark, all boundaries, and four relations that already vanish as cochains on a single choice of mark. It therefore supplies zero classes to the second differential. Arbitrary repetitions of the same rightmost exchange on slots 1 and 3 reduce to the identity and that exchange by the same polynomial. The full selected source has 62866 cochains. [Fixed rows and native identities](seventh-exchange-result.json) · [First minor](seventh-crossings-1009.json) · [Second minor](seventh-crossings-1013.json)

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -B seventh_exchange_relations.py --regenerate
~~~
